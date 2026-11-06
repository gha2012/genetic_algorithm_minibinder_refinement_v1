#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0

# version: 2025-07-03 (refactored input handling)
import os
import random
import shutil
import numpy as np
import uuid
import concurrent.futures
import multiprocessing
from pathlib import Path
import argparse
import hashlib

from Bio.PDB import PDBParser, PPBuilder, MMCIFParser, is_aa
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO
from typing import Optional, Dict, Any

# ----------------------------
# CONFIGURATION & PARAMETERS
# ----------------------------
INPUT_FOLDER = "input"
MASK_FOLDER = "masks"
POOL_FOLDER = "pool"
CHAI1_OUTPUT_FOLDER = "chai1_outputs"
POOL_SIZE = 500
GENERATIONS = 50
TOURNAMENT_SIZE = 15
MUTATION_RATE = 0.05 # base mutation rate for '+' residues in the mask
HIGH_MUTATION_RATE = 2  # for '*' residues in the mask
CROSSOVER_PROBABILITY = 0
PRUNE_FRACTION = 0.9
JOBS_PER_GPU = 2
NUM_DIFFN_SAMPLES = 2
SEQUENCE_HASH_LENGTH = 10

# PTM penalty settings
MIN_BINDER_PTM_THRESHOLD = 0.9
PTM_PENALTY_WEIGHT = 1.0

# InDel rates
INSERTION_RATE = 0
DELETION_RATE = 0

# Epitope penalty settings
#Numbering starts at 1 in chain B (target)!
HOTSPOT_RESIDUES = [81]
EPITOPE_PENALTY_WEIGHT = 0 # chai-1 
EPITOPE_DISTANCE_THRESHOLD = 5.0  # Ångstroms, not used anymore
# NEW FLAG: if False, ignore the mask and include all binder residues
USE_BINDER_MASK = True
# NEW: van der Waals radii for common elements (in Å)
VDW_RADII = {'C':1.70, 'N':1.55, 'O':1.52, 'S':1.80}
GAUSSIAN_SIGMA        = 3.0   # Ångstrom for your Gaussian fall-off, controls steepness of contact penalty: smaller = steeper, i.e. more sensitive to small distances

TARGET_PI = None
PI_PENALTY_WEIGHT = 0.0

TARGET_GRAVY = None
GRAVY_PENALTY_WEIGHT = 0.0

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

# Global chain identifiers (can be overridden via CLI)
BINDER_CHAIN = "B"
TARGET_CHAIN = "A"

# ----------------------------
# WORKER INITIALIZER
# ----------------------------
def worker_initializer():
    gpu_str = os.environ.get("CUDA_VISIBLE_DEVICES", "0")
    gpu_ids = [int(g.strip()) for g in gpu_str.split(',') if g.strip().isdigit()]
    gpu_ids = gpu_ids or [0]
    try:
        worker_id = multiprocessing.current_process()._identity[0] - 1
    except Exception:
        worker_id = 0
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_ids[worker_id % len(gpu_ids)])
    print(f"Worker {worker_id} -> GPU {os.environ['CUDA_VISIBLE_DEVICES']}", flush=True)

# ----------------------------
# MASK READING & SEQUENCE EXTRACTION
# ----------------------------
def read_mask(pdb_file: str) -> Optional[str]:
    base = Path(pdb_file).stem
    mask_path = Path(MASK_FOLDER) / f"{base}_mask.fasta"
    if not mask_path.exists():
        return None
    try:
        record = next(SeqIO.parse(str(mask_path), "fasta"), None)
    except Exception as e:
        print(f"Warning: failed to parse mask {mask_path}: {e}", flush=True)
        return None
    seq = str(record.seq) if record and record.seq else None
    if not seq or any(c not in '+-*' for c in seq):
        print(f"Invalid mask at {mask_path}, using default '+'.", flush=True)
        return None
    return seq


def extract_sequences(path: str) -> (str, str):
    p = Path(path)
    suffix = p.suffix.lower()
    if suffix in ('.fa', '.fasta', '.faa'):
        records = list(SeqIO.parse(str(p), 'fasta'))
        if len(records) != 2:
            raise ValueError(f"Expected 2 sequences in {p}, found {len(records)}")
        ids = {rec.id.lower(): str(rec.seq) for rec in records}
        return ids.get('binder', str(records[0].seq)), ids.get('target', str(records[1].seq))
    parser = PDBParser(QUIET=True) if suffix == '.pdb' else MMCIFParser(QUIET=True)
    struct = parser.get_structure('complex', str(p))
    model = struct[0]
    ppb = PPBuilder()
    binder_peps = ppb.build_peptides(model[BINDER_CHAIN])
    target_peps = ppb.build_peptides(model[TARGET_CHAIN])
    if not binder_peps or not target_peps:
        raise ValueError(f"Unable to extract both chains from {p}")
    return str(binder_peps[0].get_sequence()), str(target_peps[0].get_sequence())

# ----------------------------
# STRUCTURE LOADING & PENALTIES
# ----------------------------
def load_structure(path: str):
    print(path, flush=True)
    parser = MMCIFParser(QUIET=True) if path.lower().endswith('.cif') else PDBParser(QUIET=True)
    return parser.get_structure('model', path)


def compute_interface_centroid(coord_file: str, mask: str) -> np.ndarray:
    struct = load_structure(coord_file)
    coords = []
    seq_idx = 0
    for res in struct[0][BINDER_CHAIN]:
        if not is_aa(res):
            continue
        if seq_idx < len(mask) and mask[seq_idx] in ('+', '*') and 'CA' in res:
            coords.append(res['CA'].get_coord())
        seq_idx += 1
    if not coords:
        raise ValueError("No variable-interface residues found for centroid.")
    return np.mean(np.vstack(coords), axis=0)


def compute_epitope_penalty(coord_file: str, mask: str) -> float:
    struct = load_structure(coord_file)
    tgt_chain = struct[0][TARGET_CHAIN]
    hotspots = [tgt_chain[r]['CA'].get_coord() for r in HOTSPOT_RESIDUES if r in tgt_chain]
    if not hotspots:
        return 0.0
    cog = compute_interface_centroid(coord_file, mask)
    deviations = [max(0.0, np.linalg.norm(cog - hp) - EPITOPE_DISTANCE_THRESHOLD) for hp in hotspots]
    return EPITOPE_PENALTY_WEIGHT * (sum(deviations) / len(deviations))

def compute_surface_contact_penalty(coord_file: str, mask: str) -> float:
    """
    Penalty ∈ [0, EPITOPE_PENALTY_WEIGHT].
    Debug prints included.
    """
    print(f"[DEBUG] Loading structure from {coord_file}")
    struct = load_structure(coord_file)

    # 1) binder interface atoms
    bind_atoms = []
    seq_idx = 0
    for res in struct[0][BINDER_CHAIN]:
        if not is_aa(res):
            continue
        include = (not USE_BINDER_MASK) or (seq_idx < len(mask) and mask[seq_idx] in ('+','*'))
        if include:
            for atom in res:
                el = atom.element
                if el in VDW_RADII:
                    bind_atoms.append((atom.get_coord(), VDW_RADII[el]))
        seq_idx += 1

    print(f"[DEBUG] Collected {len(bind_atoms)} binder atoms (mask {'ON' if USE_BINDER_MASK else 'OFF'})")
    if not bind_atoms:
        print("[DEBUG] No binder atoms → returning max penalty")
        return EPITOPE_PENALTY_WEIGHT

    # 2) hotspot atoms on target
    tgt = struct[0][TARGET_CHAIN]
    hot_atoms = []
    for pos in HOTSPOT_RESIDUES:
        if pos not in tgt:
            print(f"[DEBUG] Hotspot residue {pos} not found in chain {TARGET_CHAIN}")
            continue
        for atom in tgt[pos]:
            el = atom.element
            if el in VDW_RADII:
                hot_atoms.append((atom.get_coord(), VDW_RADII[el]))
    print(f"[DEBUG] Collected {len(hot_atoms)} hotspot atoms at residues {HOTSPOT_RESIDUES}")
    if not hot_atoms:
        print("[DEBUG] No hotspot atoms → returning zero penalty")
        return 0.0

    # 3) compute minimal surface-to-surface distances
    import numpy as np
    min_dists = []
    for hi, (h_coord, h_rad) in enumerate(hot_atoms, start=1):
        dists = []
        for bi, (b_coord, b_rad) in enumerate(bind_atoms, start=1):
            center_dist = np.linalg.norm(h_coord - b_coord)
            surface_dist = max(center_dist - (h_rad + b_rad), 0.0)
            dists.append(surface_dist)
        min_d = min(dists)
        min_dists.append(min_d)
        print(f"[DEBUG] Hotspot atom #{hi}: min surface dist = {min_d:.3f} Å")

    dmins = np.array(min_dists)
    print(f"[DEBUG] Distance stats — min: {dmins.min():.3f}, mean: {dmins.mean():.3f}, max: {dmins.max():.3f}")

    # 4) gaussian contact scores
    contact_scores = np.exp(-0.5 * (dmins / GAUSSIAN_SIGMA)**2)
    mean_contact = float(contact_scores.mean())
    print(f"[DEBUG] Contact scores (Gaussian σ={GAUSSIAN_SIGMA}) — mean = {mean_contact:.3f}")

    # 5) penalty
    penalty = EPITOPE_PENALTY_WEIGHT * (1.0 - mean_contact)
    print(f"[DEBUG] Returning penalty = {penalty:.3f} (weight = {EPITOPE_PENALTY_WEIGHT})")
    return penalty

# ----------------------------
# METRICS & SCORING
# ----------------------------
from typing import Dict, Any
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def compute_custom_iptm(pae: np.ndarray,
                        idxA: list[int],
                        idxB: list[int],
                        d0: float = 10.0) -> float:
    """
    pae:  N×N array of predicted aligned errors from your Chai-1 output
    idxA: list of residue indices (0-based) in chain A you care about
    idxB: list of residue indices (0-based) in chain B you care about
    d0:   length constant, default ~10 Å
    """
    weights = []
    for i in idxA:
        for j in idxB:
            err = pae[i, j]
            w   = 1.0 / (1.0 + (err / d0)**2)
            weights.append(w)
    return float(np.mean(weights))

def compute_candidate_metrics(
    chainA_seq: str,
    candidates,
    mask: str
) -> Dict[str, Any]:
    cif_paths = candidates.cif_paths
    all_chai1_ptm_scores = [rd.ptm_scores for rd in candidates.ranking_data]
    print(all_chai1_ptm_scores, flush=True)

    # 1) Find best interface pTM
    best_iptm = 0.0
    best_idx  = -1
    for idx, scores in enumerate(all_chai1_ptm_scores):
        this_iptm = float(scores.interface_ptm)
        if this_iptm > best_iptm:
            best_iptm, best_idx = this_iptm, idx

    # 2) Pull out the winning model
    cif_file     = str(cif_paths[best_idx])
    iptm         = best_iptm
    chai1_scores = all_chai1_ptm_scores[best_idx]

    # 2.1) custom ipTM calculation
    #print(candidates, flush=True)
    #print(candidates.pae, flush=True)
    # --- usage example --------------------------------
    # Suppose 'rd' is one of your Chai-1 RankingData objects:
    #pae_matrix = candidates.pae[best_idx]  # shape (L,L)

    # define your interface residues, e.g. mask or manual list:
    # here we pick residues 10–20 on A and 100–110 on B,
    # but you can define any lists based on your mask or structure.
    #interface_A = list(range(10, 21))
    #interface_B = list(range(100, 111))

    #my_iptm = compute_custom_iptm(pae_matrix, interface_A, interface_B)
    #print("Custom ipTM over my interface:", my_iptm)

    # 3) Get per-chain pTM (chain A = index 0, chain B = index 1)
    chain_A_ptm = float(chai1_scores.per_chain_ptm[0][0])
    chain_B_ptm = float(chai1_scores.per_chain_ptm[0][1])

    # 4) PTM penalty
    delta_ptm    = max(0.0, MIN_BINDER_PTM_THRESHOLD - chain_A_ptm) \
                   if MIN_BINDER_PTM_THRESHOLD is not None else 0.0
    ptm_penalty  = PTM_PENALTY_WEIGHT * delta_ptm

    # 5) pI penalty
    pa           = ProteinAnalysis(chainA_seq)
    pI           = pa.isoelectric_point()
    delta_pI     = abs(pI - TARGET_PI) if TARGET_PI is not None else 0.0
    pI_penalty   = PI_PENALTY_WEIGHT * delta_pI

    # 6) GRAVY penalty
    gravy        = pa.gravy()
    delta_gravy  = abs(gravy - TARGET_GRAVY) if TARGET_GRAVY is not None else 0.0
    gravy_penalty= GRAVY_PENALTY_WEIGHT * delta_gravy

    # 7) Epitope penalty
    print(cif_file, mask, flush=True)
    epi_pen = compute_surface_contact_penalty(cif_file, mask)

    # 8) Final score (clamped at zero)
    score = iptm - (ptm_penalty + pI_penalty + gravy_penalty + epi_pen)
    if score < 0:
        print(
            f"Warning: negative score {score:.4f} for {cif_file}, "
            f"iptm={iptm:.4f}, penalties: ptm={ptm_penalty:.4f}, "
            f"pI={pI_penalty:.4f}, gravy={gravy_penalty:.4f}, epitope={epi_pen:.4f}",
            flush=True
        )
        #score = 0.0

    return {
        'cif_file':        cif_file,
        'chai1_scores':    chai1_scores,
        'iptm':            iptm,
        'chain_A_ptm':     chain_A_ptm,
        'chain_B_ptm':     chain_B_ptm,
        'pI':              pI,
        'delta_pI':        delta_pI,
        'gravy':           gravy,
        'delta_gravy':     delta_gravy,
        'epitope_penalty': epi_pen,
        'score':           score
    }

# ----------------------------
# CHAI-1 INTERFACE
# ----------------------------
def score_chai1(chainA_seq: str, chainB_seq: str, design_name: str,
                device: Optional[str]=None) -> str:
    device = device or os.environ.get('CHAI1_DEVICE', 'cuda:0')
    design_folder = Path(CHAI1_OUTPUT_FOLDER) / design_name
    design_folder.mkdir(parents=True, exist_ok=True)
    fasta_path = design_folder / f"{design_name}.fasta"
    with open(fasta_path, 'w') as f:
        f.write(f">protein|name=binder\n{chainA_seq}\n")
        f.write(f">protein|name=target\n{chainB_seq}\n")
    outputs = design_folder / 'outputs'
    shutil.rmtree(outputs, ignore_errors=True)
    outputs.mkdir()
    from chai_lab.chai1 import run_inference
    candidates = run_inference(
        fasta_file=fasta_path,
        output_dir=outputs,
        num_trunk_recycles=3,
        num_diffn_timesteps=200,
        num_diffn_samples=NUM_DIFFN_SAMPLES,
        seed=42,
        device=device,
        use_esm_embeddings=True
    )
    #npz_files = list(outputs.glob('scores.model_idx_*.npz'))
    #cif_files = list(outputs.glob('pred.model_idx_*.cif'))
    return candidates

# ----------------------------
# MUTATION & RECOMBINATION
# ----------------------------
def mutate_with_mask(seq: str, mask: str, base_rate: float=MUTATION_RATE,
                     high_factor: float=HIGH_MUTATION_RATE) -> (str, str):
    n = len(seq)
    allowed_insert = [False] * (n+1)
    for i in range(n+1):
        if i == 0:
            allowed_insert[i] = (mask[0] != '-')
        elif i == n:
            allowed_insert[i] = (mask[-1] != '-')
        else:
            allowed_insert[i] = not (mask[i-1] == '-' and mask[i] == '-')
    new_seq, new_mask = [], []
    for i, aa in enumerate(seq):
        if allowed_insert[i] and random.random() < INSERTION_RATE:
            new_seq.append(random.choice(AMINO_ACIDS)); new_mask.append('+')
        if mask[i] in ('+', '*') and random.random() < DELETION_RATE:
            continue
        if mask[i] != '-' and random.random() < base_rate * (high_factor if mask[i]=='*' else 1.0):
            aa = random.choice(AMINO_ACIDS)
        new_seq.append(aa); new_mask.append(mask[i])
    if allowed_insert[n] and random.random() < INSERTION_RATE:
        new_seq.append(random.choice(AMINO_ACIDS)); new_mask.append('+')
    return ''.join(new_seq), ''.join(new_mask)

def crossover_with_mask(seq1: str, seq2: str, mask1: str, mask2: str) -> (str, str):
    if len(seq1)<2 or len(seq2)<2: return seq1, mask1
    pt = random.randint(1, min(len(seq1), len(seq2))-1)
    return seq1[:pt]+seq2[pt:], mask1[:pt]+mask2[pt:]

# ----------------------------
# POOL MANAGEMENT
# ----------------------------
def get_score(cand: Dict[str, Any]) -> float: return cand.get('score', -np.inf)
def tournament_selection(pool: list, k: int=TOURNAMENT_SIZE) -> Dict[str, Any]: return max(random.sample(pool, k), key=get_score)

def create_offspring_candidate(pool: list, gen: int, var: int) -> Dict[str, Any]:
    p1, p2 = tournament_selection(pool), tournament_selection(pool)
    if random.random()<CROSSOVER_PROBABILITY:
        seq,mask = crossover_with_mask(p1['chainA_seq'],p2['chainA_seq'],p1['mask'],p2['mask'])
        parents=[p1['design_name'],p2['design_name']]
    else:
        parent = random.choice([p1,p2]); seq,mask = parent['chainA_seq'],parent['mask']; parents=[parent['design_name']]
    seq,mask=mutate_with_mask(seq,mask)
    name=f"gen{gen}_{var}_variant_{uuid.uuid4().hex[:8]}"
    return {'chainA_seq':seq,'chainB_seq':p1['chainB_seq'],'mask':mask,'generation':gen,
            'progenitors':p1.get('progenitors',[])+p2.get('progenitors',[])+[parents],
            'design_name':name}

def process_input_file(input_file: str)->Optional[Dict[str,Any]]:
    try: binder, targ = extract_sequences(input_file)
    except Exception as e: print(f"Skipping {input_file}: {e}",flush=True); return None
    mask=read_mask(input_file) or '+'*len(binder)
    if len(mask)!=len(binder): raise ValueError(f"Mask length mismatch for {input_file}")
    return {'chainA_seq':binder,'chainB_seq':targ,'mask':mask,
            'generation':0,'progenitors':[],'design_name':Path(input_file).stem}

def score_candidate(candidate: Dict[str,Any])->Dict[str,Any]:
    candidates = score_chai1(candidate['chainA_seq'],candidate['chainB_seq'],candidate['design_name'])
    candidate.update(compute_candidate_metrics(candidate['chainA_seq'],candidates,candidate['mask']))
    return candidate

def process_input_wrapper(input_file: str)->Optional[Dict[str,Any]]:
    cand=process_input_file(input_file);
    return None if cand is None else score_candidate(cand)

def prune_pool(pool:list)->list: return sorted(pool,key=get_score,reverse=True)[:max(1,int(len(pool)*PRUNE_FRACTION))]

def save_candidate(cand:Dict[str,Any],fn:str):
    seq = cand['chainA_seq']
    h = hashlib.sha1(seq.encode()).hexdigest()[:SEQUENCE_HASH_LENGTH]
    hdr = (
        f"# REMARK length: {len(seq)}\n"
        f"# REMARK hash: {h}\n"
        f"# REMARK iptm: {cand['iptm']:.4f}\n"
        f"# REMARK chainA_ptm: {cand['chain_A_ptm']:.4f}\n"
        f"# REMARK chainB_ptm: {cand['chain_B_ptm']:.4f}\n"
        f"# REMARK pI:   {cand['pI']:.2f}\n"
        f"# REMARK GRAVY:{cand['gravy']:.3f}\n"
        f"# REMARK epitope:  {cand['epitope_penalty']:.4f}\n"
        f"# REMARK score:{cand['score']:.4f}\n"
    )
    try:
        cont = open(cand['cif_file']).read()
    except Exception:
        cont = ''
    with open(fn, 'w') as out:
        out.write(hdr + cont)

def update_pool_folder(pool:list):
    shutil.rmtree(POOL_FOLDER,ignore_errors=True); os.makedirs(POOL_FOLDER,exist_ok=True)
    for c in pool: save_candidate(c, str(Path(POOL_FOLDER)/f"{c['design_name']}.cif"))

def write_generation_stats(gen:int,pool:list,fn:str):
    scores=[get_score(c) for c in pool]; hist,bins=np.histogram(scores,bins=20,range=(-1,1))
    maxc=hist.max() if hist.size>0 else 0; scale=50/maxc if maxc>0 else 1; lines=[f"Gen {gen} stats"]
    for i,count in enumerate(hist): lines.append(f"[{bins[i]:.2f},{bins[i+1]:.2f}]: {'*'*int(count*scale)} ({count})")
    with open(fn,'w') as f: f.write("\n".join(lines))

def write_top_designs(gen:int, pool:list, fn:str):
    top = sorted(pool, key=get_score, reverse=True)[:20]
    name_w = max(len(c['design_name']) for c in top)
    # Add chainA_ptm, chainB_ptm, and epitope_penalty columns after iPTM
    fmt = f"{{:<5}}{{:<{name_w}}}  {{:<6}}  {{:<8}}  {{:<8}}  {{:<8}}  {{:<8}}  {{:<8}}  {{:<6}}  {{:<6}}  {{:<6}}\n"
    with open(fn, 'a') as f:
        f.write("\nTop 20 Designs:\n")
        header = fmt.format("Rank", "Name", "Len", "Hash", "Score", "iPTM", "A-PTM", "B-PTM", "Epitope", "pI", "GRAVY")
        f.write(header)
        f.write("-" * (len(header) - 1) + "\n")
        for i, c in enumerate(top, 1):
            h = hashlib.sha1(c['chainA_seq'].encode()).hexdigest()[:8]
            f.write(fmt.format(
                i,
                c['design_name'],
                len(c['chainA_seq']),
                h,
                f"{c['score']:.4f}",
                f"{c['iptm']:.4f}",
                f"{c['chain_A_ptm']:.4f}",
                f"{c['chain_B_ptm']:.4f}",
                f"{c['epitope_penalty']:.4f}",
                f"{c['pI']:.4f}",
                f"{c['gravy']:.4f}"
            ))

# ----------------------------
# MAIN EVOLUTIONARY LOOP
# ----------------------------


def main():
    # Count GPUs visible inside the Slurm job environment
    visible_gpus = [x.strip() for x in os.environ.get("CUDA_VISIBLE_DEVICES", "0").split(",") if x.strip()]
    ngpus = max(1, len(visible_gpus))
    workers = max(1, ngpus * JOBS_PER_GPU)

    print(f"CUDA_VISIBLE_DEVICES={os.environ.get('CUDA_VISIBLE_DEVICES', '')}", flush=True)
    print(f"Visible GPUs: {visible_gpus}", flush=True)
    print(f"Workers: {workers}", flush=True)

    inputs = [
        str(p) for p in Path(INPUT_FOLDER).iterdir()
        if p.suffix.lower() in (".pdb", ".cif", ".fa", ".fasta")
    ]
    if not inputs:
        raise RuntimeError("No valid files in input folder.")

    # CUDA-safe multiprocessing context
    ctx = multiprocessing.get_context("forkserver")

    pool = []
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=workers,
        mp_context=ctx,
        initializer=worker_initializer,
    ) as ex:
        for c in ex.map(process_input_wrapper, inputs):
            if c:
                pool.append(c)

    if len(pool) < POOL_SIZE:
        founders = pool.copy()
        offspring = []
        idx = 1

        while len(pool) + len(offspring) < POOL_SIZE:
            parent = random.choice(founders)
            seq, mask = mutate_with_mask(parent["chainA_seq"], parent["mask"])
            name = f"gen0_{idx}_variant_{uuid.uuid4().hex[:8]}"
            offspring.append({
                "chainA_seq": seq,
                "mask": mask,
                "chainB_seq": parent["chainB_seq"],
                "generation": 0,
                "progenitors": [parent["design_name"]],
                "design_name": name,
            })
            idx += 1

        with concurrent.futures.ProcessPoolExecutor(
            max_workers=workers,
            mp_context=ctx,
            initializer=worker_initializer,
        ) as ex:
            offspring = list(ex.map(score_candidate, offspring))

        pool += offspring

    pool = pool[:POOL_SIZE]
    update_pool_folder(pool)
    write_generation_stats(0, pool, "gen_0_stats.txt")
    write_top_designs(0, pool, "gen_0_stats.txt")

    for gen in range(1, GENERATIONS + 1):
        pool = prune_pool(pool)
        best = max(pool, key=get_score)
        print(f"Gen {gen}: trimmed to {len(pool)}; best {best['score']:.4f}", flush=True)

        offspring = []
        cnt = 1
        while len(pool) + len(offspring) < POOL_SIZE:
            offspring.append(create_offspring_candidate(pool, gen, cnt))
            cnt += 1

        with concurrent.futures.ProcessPoolExecutor(
            max_workers=workers,
            mp_context=ctx,
            initializer=worker_initializer,
        ) as ex:
            offspring = list(ex.map(score_candidate, offspring))

        pool += offspring
        update_pool_folder(pool)
        write_generation_stats(gen, pool, f"gen_{gen}_stats.txt")
        write_top_designs(gen, pool, f"gen_{gen}_stats.txt")

    best = max(pool, key=get_score)
    print("\nEvolution complete.")
    print(f"Best: {best['design_name']} (score {best['score']:.4f})")
    print("Binder seq:", best["chainA_seq"])
    print("Target seq:", best["chainB_seq"])
    print("CIF:", best["cif_file"])


if __name__ == "__main__":
    multiprocessing.set_start_method("forkserver", force=True)

    parser = argparse.ArgumentParser(description="Evolutionary protein binder design")
    parser.add_argument("--binder", type=str, default="A")
    parser.add_argument("--target", type=str, default="B")
    parser.add_argument("--target-pi", type=float, default=None)
    parser.add_argument("--pi-penalty-weight", type=float, default=0.0)
    parser.add_argument("--target-gravy", type=float, default=None)
    parser.add_argument("--gravy-penalty-weight", type=float, default=0.0)
    parser.add_argument("--insertion-rate", type=float, default=INSERTION_RATE)
    parser.add_argument("--deletion-rate", type=float, default=DELETION_RATE)
    parser.add_argument("--epitope-penalty-weight", type=float, default=EPITOPE_PENALTY_WEIGHT)
    parser.add_argument("--epitope-threshold", type=float, default=EPITOPE_DISTANCE_THRESHOLD)
    args = parser.parse_args()

    BINDER_CHAIN = args.binder.upper()
    TARGET_CHAIN = args.target.upper()
    TARGET_PI = args.target_pi
    PI_PENALTY_WEIGHT = args.pi_penalty_weight
    TARGET_GRAVY = args.target_gravy
    GRAVY_PENALTY_WEIGHT = args.gravy_penalty_weight
    INSERTION_RATE = args.insertion_rate
    DELETION_RATE = args.deletion_rate
    EPITOPE_PENALTY_WEIGHT = args.epitope_penalty_weight
    EPITOPE_DISTANCE_THRESHOLD = args.epitope_threshold

    main()
