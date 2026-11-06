# Genetic algorithm for minibinder refinement

GPU-oriented tools for diversifying protein minibinders and prioritizing candidates predicted to bind a target.

## Two alternative workflows

This repository contains **two separate, standalone approaches**. Choose one for a given run:

- **Option I — [`design_genetic.sh`](design_genetic.sh):** runs a SLURM-oriented structural diversification loop using [RFdiffusion](https://github.com/RosettaCommons/RFdiffusion), [ProteinMPNN/Rosetta and AlphaFold2 initial guess](https://github.com/nrbennet/dl_binder_design).
- **Option II — [`chai1_genetic_ptm_multiGPU.py`](chai1_genetic_ptm_multiGPU.py):** evolves binder sequences, predicts each binder-target complex with [Chai-1](https://github.com/chaidiscovery/chai-lab), and ranks candidates using interface pTM plus optional physicochemical and epitope penalties.

These scripts are not stages of one pipeline, do not exchange results, and are not intended to run in conjunction. They have different dependencies, algorithms, input assumptions, and outputs. Use separate working directories and software environments. Both happen to use a directory named `pool`; the Python option deletes and rebuilds that directory during each generation.

`MiSeqAnalysisClean_v4.R` is a third, unrelated, interactive MiSeq analysis script. Neither refinement option calls it.

## Repository layout

```text
.
├── design_genetic.sh              # Option I: RFdiffusion/MPNN/AF2 loop
├── chai1_genetic_ptm_multiGPU.py  # Option II: Chai-1 sequence evolution
└── MiSeqAnalysisClean_v4.R        # separate sequencing analysis
```

## Option I: RFdiffusion/ProteinMPNN/AlphaFold2 shell workflow

Option I is independent of Option II. It does not use Chai-1, Option II's masks, its CIF population, or its scoring files.

### What it does

For 50 iterations, `design_genetic.sh`:

1. fills a structural pool to 500 members by partially diffusing randomly selected PDB complexes;
2. packs generated PDBs into a Rosetta silent file;
3. designs interface sequences with `dl_interface_design.py`;
4. predicts and scores complexes with the AlphaFold2 initial-guess pipeline;
5. extracts the predicted structures; and
6. deletes the half with the worst `pae_interaction` values.

This script is a site-specific SLURM job template, not a portable installer or general command-line program. It accepts no command-line options; paths, environment names, iteration count, and pool size must be edited in the file.

### Dependencies

Option I uses the following external software and runtime tools:

| Dependency | How it is used |
|---|---|
| SLURM | Schedules the GPU job and provides `SLURM_JOB_ID`. The current script requests one `a5000` GPU. |
| Bash | Runs the iteration, random-selection, file-management, and pruning logic. |
| Conda | Activates the three hard-coded environments: `SE3nv`, `proteinmpnn_binder_design`, and `af2_binder_design`. |
| [RFdiffusion](https://github.com/RosettaCommons/RFdiffusion) and its model weights | Partially diffuses selected binder-target PDB structures to generate new backbones. |
| [`dl_binder_design`](https://github.com/nrbennet/dl_binder_design) | Supplies `mpnn_fr/dl_interface_design.py` for ProteinMPNN/FastRelax interface sequence design and `af2_initial_guess/predict.py` for AlphaFold2 screening. |
| ProteinMPNN and Rosetta/PyRosetta | Perform interface sequence design and structural relaxation through `dl_binder_design`. |
| AlphaFold2 model parameters and required databases | Predict and score candidate complexes through the initial-guess workflow. |
| Rosetta silent-file utilities | `silentfrompdbs` packages structures and `silentextract` restores predicted PDB files. Both must be available on `PATH`. |
| Standard Unix tools | `uuidgen`, `grep`, `awk`, `sort`, `find`, `cp`, `mv`, `rm`, and `wc` are called directly by the script. |

System and hardware dependencies are:

- A Linux SLURM cluster with an NVIDIA GPU partition.
- A compatible NVIDIA driver and the CUDA/PyTorch stacks required by RFdiffusion, ProteinMPNN, and AlphaFold2.
- Storage for model weights, AlphaFold databases, Rosetta silent files, PDB structures, and intermediate predictions.

Install and test RFdiffusion and `dl_binder_design` from their upstream instructions before adapting this driver. Their CUDA, PyTorch, Rosetta/PyRosetta, AlphaFold, database, and model-weight requirements are substantial and cluster-specific.

### Site-specific configuration

At minimum, review these lines in `design_genetic.sh`:

```bash
#SBATCH --gpus=a5000:1

poolFolder="pool"
iterations=50
poolSize=500

source /home/mhoe/miniconda3/etc/profile.d/conda.sh

/home/${USER}/RFdiffusion/scripts/run_inference.py ...
/home/${USER}/dl_binder_design/mpnn_fr/dl_interface_design.py ...
/home/${USER}/dl_binder_design/af2_initial_guess/predict.py ...
```

The hard-coded `/home/mhoe/...` path will not work for other users. Change it to the conda initialization script used on your cluster. Confirm that the environments are named exactly:

```text
SE3nv
proteinmpnn_binder_design
af2_binder_design
```

You will usually also need cluster-specific SLURM directives for wall time, CPUs, memory, partition/account, and log location.

### Inputs

Prepare an isolated Option I run directory:

```text
run_option_i/
├── design_genetic.sh
├── input/
│   ├── seed_01.pdb
│   └── seed_02.pdb
└── pool/
```

Run the job from `run_option_i/`. The script assumes:

- at least one `input/*.pdb` file;
- chain `A` is the binder to diversify;
- chain `B` is the fixed target;
- ATOM records use conventional whitespace and numeric residue IDs that the script can read with `awk '{print $6}'`; and
- the target chain's first-to-last residue range is valid as one RFdiffusion contig.

The script copies all input PDBs into `pool/` and chooses parents from the PDBs currently there.

### Usage

After adapting the paths and scheduler settings:

```bash
sbatch design_genetic.sh
```

For an initial cluster test, set `iterations=1` and a small `poolSize` such as `4`. This still runs the complete Option I chain.

### Outputs

The Option I `pool/` directory accumulates iteration-specific artifacts, including:

```text
inference_output_iteration_<N>.silent
interfaceDesignOutput_iteration_<N>.silent
afCheck_iteration_<N>.silent
afCheck_iteration_<N>.sc
*.pdb                                  # retained extracted predictions/seeds
slurm-<job_id>.out                     # moved here at completion
```

RFdiffusion `.trb` files and intermediate poly-glycine PDBs are deleted. Approximately the worst half of extracted PDBs, ranked by `pae_interaction`, are also deleted during every iteration.

### Current assumptions and caveats

- The shell script has no strict error mode, so a failed external command may not stop later destructive cleanup commands.
- The input-glob test does not reliably detect an empty `input/` directory under default Bash glob behavior.
- The selection loop expects every `*.pdb` considered for sorting to contain a `pae_interaction` record. Copied founder PDBs may not have one, so inspect how blank values sort on your system.
- The `files_array` selection array is declared inside the iteration loop but is not explicitly cleared, which can carry stale entries into later iterations.
- `check.point` is deleted at startup, but no checkpoint or resume mechanism is implemented.
- The final `mv ../slurm-${job_id}.out .` assumes SLURM's default log name and that the log is created one directory above `pool/`.

Use a disposable run directory, retain backups of founder structures, and validate a one-iteration run before scaling up.

## Option II: Chai-1 sequence evolution

### What it does

For each founder complex in `input/`, the Python program:

1. extracts a binder and target sequence;
2. applies an optional per-residue mutation mask;
3. predicts the complex with Chai-1 and keeps the sample with the best interface pTM (`ipTM`);
4. fills a population of 500 candidates;
5. repeats 50 generations of tournament selection, mutation, optional crossover, Chai-1 prediction, and pruning; and
6. writes the current population as CIF files plus a text summary for every generation.

The default run is a large production-scale job: approximately 3,000 candidate evaluations and, with two Chai-1 diffusion samples per evaluation, up to approximately 6,000 predicted structures. Start with the smoke-test settings below and estimate cluster time and storage before using the defaults.

### Dependencies

Direct runtime dependencies used by Option II are:

| Dependency | How it is used |
|---|---|
| Python 3.10 or later | Runs the evolutionary algorithm and multiprocessing code. |
| [`chai_lab`](https://github.com/chaidiscovery/chai-lab) | Predicts binder-target structures and supplies ipTM and per-chain pTM scores. Its installation brings its own dependencies, including PyTorch. |
| [`numpy`](https://numpy.org/) | Coordinate calculations, Gaussian contact scoring, score histograms, and numerical utilities. |
| [Biopython](https://biopython.org/) | Reads PDB, mmCIF, FASTA, and mask files; extracts peptide sequences; calculates pI and GRAVY. |

System and hardware dependencies are:

- Linux.
- An NVIDIA GPU, compatible NVIDIA driver, and CUDA runtime with bfloat16 support. Chai-1's upstream documentation recommends high-memory data-center GPUs; smaller complexes may fit on lower-memory cards.
- Sufficient local storage for Chai-1 weights, per-candidate FASTA/CIF outputs, and generation pools.
- Network access on the first run so Chai-1 can obtain its model weights, unless those weights are already cached.
- Conda or another Python environment manager is recommended for isolation, but the script does not import it directly.

The repository does not provide a lock file or state which dependency versions were tested. The pinned Chai-1 release below is therefore a reproducible starting point, not a claim of author-tested compatibility.

### Installation

```bash
git clone https://github.com/gha2012/genetic_algorithm_minibinder_refinement_v1.git
cd genetic_algorithm_minibinder_refinement_v1

conda create -n minibinder-chai python=3.10 -y
conda activate minibinder-chai

python -m pip install --upgrade pip
python -m pip install "chai_lab==0.6.1" biopython numpy

python -c "from chai_lab.chai1 import run_inference; import Bio, numpy; print('imports OK')"
```

Chai-1 documents `chai_lab==0.6.1` as its published package version at the time this README was prepared. If the script fails because the Chai-1 return structure has changed, test against the current upstream release and pin the exact working version in your environment.

By default, Chai-1 manages downloaded weights itself. To place them on a larger or shared filesystem:

```bash
export CHAI_DOWNLOADS_DIR=/path/to/chai-downloads
```

### Inputs

Create an `input/` directory in the Option II working directory. Each file may be a PDB, mmCIF, or two-record FASTA file:

```text
input/
├── seed_01.pdb
├── seed_02.cif
└── seed_03.fasta
```

#### PDB and mmCIF

With the defaults, chain `A` is the binder and chain `B` is the target. Override those chain identifiers with `--binder` and `--target`.

Only the first continuous peptide assembled by Biopython's `PPBuilder` is used from each chain. Structures with chain breaks, unusual residue numbering, or nonstandard amino acids should be checked carefully.

#### FASTA

A FASTA input must contain exactly two records. Prefer the exact record IDs `binder` and `target`:

```fasta
>binder
MKKL...
>target
EVQL...
```

If those IDs are not present, the first record is treated as the binder and the second as the target. The `--binder` and `--target` options affect structural inputs; they do not reorder FASTA records.

#### Mutation masks

Masks are optional. A mask for `input/seed_01.pdb` must be named `masks/seed_01_mask.fasta` and have exactly the same length as the extracted binder sequence:

```text
masks/
└── seed_01_mask.fasta
```

```fasta
>seed_01_mask
------++++++++****++++------
```

| Symbol | Meaning |
|---|---|
| `-` | Keep the residue fixed; do not substitute or delete it. |
| `+` | Mutate with the base substitution probability, 0.05 by default. |
| `*` | Mutate at twice the base probability, 0.10 with the current defaults. |

If a mask is absent or invalid, all binder positions receive `+`. Insertions and deletions are disabled by default. When enabled, inserted residues receive `+`; insertions are suppressed at boundaries surrounded by fixed residues.

### Usage

Single GPU, default scoring:

```bash
mkdir -p input masks
export CUDA_VISIBLE_DEVICES=0
python chai1_genetic_ptm_multiGPU.py
```

Two GPUs with pI, hydropathy, and indel penalties enabled:

```bash
export CUDA_VISIBLE_DEVICES=0,1
python chai1_genetic_ptm_multiGPU.py \
  --binder A \
  --target B \
  --target-pi 6.5 \
  --pi-penalty-weight 0.02 \
  --target-gravy -0.25 \
  --gravy-penalty-weight 0.10 \
  --insertion-rate 0.002 \
  --deletion-rate 0.002
```

Use numeric, comma-separated CUDA device identifiers. The program starts two worker processes per visible GPU by default, and each worker performs its own Chai-1 inference.

#### Command-line options

| Option | Default | Effect |
|---|---:|---|
| `--binder ID` | `A` | Binder chain identifier used when parsing PDB/mmCIF founders. |
| `--target ID` | `B` | Target chain identifier used when parsing PDB/mmCIF founders. |
| `--target-pi FLOAT` | unset | Desired binder isoelectric point. No pI penalty is applied when unset. |
| `--pi-penalty-weight FLOAT` | `0.0` | Multiplier for `abs(observed_pI - target_pI)`. |
| `--target-gravy FLOAT` | unset | Desired binder GRAVY hydropathy value. |
| `--gravy-penalty-weight FLOAT` | `0.0` | Multiplier for `abs(observed_GRAVY - target_GRAVY)`. |
| `--insertion-rate FLOAT` | `0.0` | Per-eligible-boundary insertion probability. Use a value in `[0, 1]`. |
| `--deletion-rate FLOAT` | `0.0` | Per-variable-residue deletion probability. Use a value in `[0, 1]`. |
| `--epitope-penalty-weight FLOAT` | `0.0` | Scales the active Gaussian surface-contact penalty. A zero value disables it. |
| `--epitope-threshold FLOAT` | `5.0` | Parsed and stored, but **does not affect the active scoring path in the current code**. The active contact penalty uses the hard-coded `GAUSSIAN_SIGMA` instead. |

Several important settings are source constants rather than command-line options:

| Constant | Default | Purpose |
|---|---:|---|
| `POOL_SIZE` | `500` | Candidates retained per generation. |
| `GENERATIONS` | `50` | Evolutionary generations after generation 0. |
| `TOURNAMENT_SIZE` | `15` | Candidates sampled during parent selection. |
| `MUTATION_RATE` | `0.05` | Base substitution probability at `+` positions. |
| `CROSSOVER_PROBABILITY` | `0` | Probability of single-point crossover; disabled by default. |
| `PRUNE_FRACTION` | `0.9` | Fraction retained before refilling each generation. |
| `JOBS_PER_GPU` | `2` | Concurrent worker processes assigned per visible GPU. |
| `NUM_DIFFN_SAMPLES` | `2` | Chai-1 diffusion samples per candidate. |
| `HOTSPOT_RESIDUES` | `[81]` | One-based target residue IDs used by the contact penalty. |
| `GAUSSIAN_SIGMA` | `3.0 Å` | Width of the active surface-contact penalty. |

For a first end-to-end smoke test, temporarily use values such as:

```python
POOL_SIZE = 8
GENERATIONS = 1
TOURNAMENT_SIZE = 4
JOBS_PER_GPU = 1
NUM_DIFFN_SAMPLES = 1
```

`TOURNAMENT_SIZE` must not exceed the population remaining after pruning.

### Scoring

For each candidate, the program selects the Chai-1 sample with the highest interface pTM and calculates:

```text
score = ipTM
        - binder-pTM penalty
        - pI penalty
        - GRAVY penalty
        - epitope-contact penalty
```

The binder-pTM penalty starts below a hard-coded threshold of `0.9` and has a default weight of `1.0`. The pI, GRAVY, and epitope penalties are disabled by default. Negative total scores are retained rather than clamped to zero.

The epitope term evaluates surface contact between masked binder atoms and `HOTSPOT_RESIDUES` on the target. If you enable it, verify the target residue IDs and chain IDs in Chai-1's output CIFs before relying on the ranking. Non-default input chain identifiers may not be preserved by Chai-1 in the same way the scoring code assumes.

### Outputs

```text
chai1_outputs/
└── <design_name>/
    ├── <design_name>.fasta
    └── outputs/                  # Chai-1 samples and score data

pool/
└── <design_name>.cif             # current population with REMARK metrics

gen_0_stats.txt
gen_1_stats.txt
...
gen_50_stats.txt
```

Each population CIF is prefixed with the sequence length, sequence hash, ipTM, per-chain pTM, pI, GRAVY, epitope penalty, and final score. Each generation summary contains a score histogram and the top 20 designs.

The program does not checkpoint its in-memory population and does not resume a partial run. `chai1_outputs/` is retained, but restarting begins a new evolution and rewrites `pool/` and generation summaries.

## Reproducibility notes

- Option II's genetic operations use Python's unseeded random generator; reruns are not deterministic. Chai-1 inference itself is called with `seed=42`.
- Option I uses Bash's unseeded `$RANDOM` for parent and diffusion-depth selection.
- Neither option records dependency versions, source revisions, GPU/software details, or a complete run configuration. For a real campaign, archive those alongside the outputs.
- These scores are computational prioritization signals, not experimental evidence of binding, specificity, stability, safety, or biological activity. Validate selected designs experimentally.

## Citation

If you use this code, cite:

> Tan, E., Frechen, S. C., Broske, B., et al. **Evolutionary algorithms accelerate de novo design of potent Nectin-4-specific cancer biologics.** *bioRxiv* (2026). <https://doi.org/10.64898/2026.03.04.709551>

This manuscript is currently under revision. The citation will be updated when the peer-reviewed version is published.

## License

The original code in this repository is free to use, modify, and redistribute under the [Apache License 2.0](LICENSE). See `LICENSE` for the complete terms, including attribution, redistribution, patent, and warranty provisions.

This license applies only to code authored in this repository. It does not license or relicense external software, model weights, parameters, databases, or other assets used by either option. Those components remain subject to their own licenses and terms.

Key third-party terms include:

- **Chai-1:** code and model weights are distributed under Apache-2.0.
- **RFdiffusion:** source code and referenced model weights are distributed under BSD-3-Clause.
- **AlphaFold2:** source code is distributed under Apache-2.0; model parameters are distributed under CC BY 4.0.
- **`dl_binder_design`:** distributed under the MIT License.
- **ProteinMPNN:** distributed under the MIT License.
- **Rosetta/PyRosetta:** governed by the Rosetta license; commercial use requires a separate commercial license.

See [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md) for links and additional details. Users are responsible for complying with the terms that apply to the exact software versions, model assets, and databases they use.

## Upstream projects

- [Chai-1 / `chai_lab`](https://github.com/chaidiscovery/chai-lab)
- [RFdiffusion](https://github.com/RosettaCommons/RFdiffusion)
- [`dl_binder_design`](https://github.com/nrbennet/dl_binder_design)
- [Biopython](https://biopython.org/)
