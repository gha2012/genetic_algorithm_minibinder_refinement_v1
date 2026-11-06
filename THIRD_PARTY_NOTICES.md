# Third-party software notices

The Apache License 2.0 in this repository covers only the original code distributed in this repository. No third-party source code, model weights, model parameters, or biological databases are included here. The scripts import or invoke software installed separately by the user.

Use of each external component is governed by the license and terms supplied with the exact version and assets installed by the user. Upstream terms control if they differ from this summary.

## Option I: `design_genetic.sh`

### RFdiffusion

- Project: <https://github.com/RosettaCommons/RFdiffusion>
- License: BSD 3-Clause
- License text: <https://github.com/RosettaCommons/RFdiffusion/blob/main/LICENSE>
- Scope stated upstream: source code and model weights referenced for download by the RFdiffusion README.

### `dl_binder_design`

- Project: <https://github.com/nrbennet/dl_binder_design>
- License: MIT
- License text: <https://github.com/nrbennet/dl_binder_design/blob/main/LICENSE>

The shell workflow calls `mpnn_fr/dl_interface_design.py` and `af2_initial_guess/predict.py` from this project.

### ProteinMPNN

- Project: <https://github.com/dauparas/ProteinMPNN>
- License: MIT
- License text: <https://github.com/dauparas/ProteinMPNN/blob/main/LICENSE>

### AlphaFold2

- Project: <https://github.com/google-deepmind/alphafold>
- Source-code license: Apache License 2.0
- Model-parameter license: Creative Commons Attribution 4.0 International (CC BY 4.0)
- License and disclaimer: <https://github.com/google-deepmind/alphafold#license-and-disclaimer>

AlphaFold code, model parameters, prediction data, and external sequence or structure databases can have distinct terms.

### Rosetta and PyRosetta

- Project: <https://github.com/RosettaCommons/rosetta>
- License information: <https://github.com/RosettaCommons/rosetta/blob/main/LICENSE.md>

Rosetta and PyRosetta are not distributed under this repository's Apache-2.0 license. Rosetta is available under its own University of Washington license. Its upstream documentation states that commercial use requires a separate commercial license.

### Additional Option I components

The CUDA toolkit, NVIDIA drivers, PyTorch, Conda, SLURM, Unix utilities, downloaded databases, and any other transitive dependencies remain subject to their respective licenses and terms.

## Option II: `chai1_genetic_ptm_multiGPU.py`

### Chai-1 / `chai_lab`

- Project: <https://github.com/chaidiscovery/chai-lab>
- License: Apache License 2.0
- License text: <https://github.com/chaidiscovery/chai-lab/blob/main/LICENSE>
- Scope stated upstream: Chai-1 code and model weights.

### NumPy

- Project: <https://github.com/numpy/numpy>
- License information: <https://github.com/numpy/numpy/blob/main/LICENSE.txt>

### Biopython

- Project: <https://github.com/biopython/biopython>
- License information: <https://github.com/biopython/biopython/blob/master/LICENSE.rst>

### Additional Option II components

PyTorch and the other transitive dependencies installed with `chai_lab`, as well as the CUDA toolkit, NVIDIA drivers, Conda, and downloaded assets, remain subject to their respective licenses and terms.
