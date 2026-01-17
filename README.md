# psi4me

A Python repository for notebook-driven quantum chemistry workflows combining [Psi4](https://psicode.org/) and [RDKit](https://www.rdkit.org/).

## Overview

This repository provides a clean, modular structure for conducting quantum chemistry calculations and molecular structure analysis through Jupyter notebooks. Reusable helper functions are organized in the `src/` module, making it easy to maintain consistent workflows across multiple notebooks.

## Repository Structure

```
psi4me/
├── notebooks/          # Jupyter notebooks for workflows and analysis
│   └── 01_basic_workflow.ipynb
├── src/                # Reusable helper modules
│   ├── __init__.py
│   ├── molecular_helpers.py      # RDKit molecular structure utilities
│   ├── quantum_helpers.py        # Psi4 quantum chemistry calculations
│   └── visualization_helpers.py  # Plotting and visualization tools
├── data/               # Input data and molecular structures
├── figures/            # Output figures and plots
├── environment.yml     # Conda environment specification
└── README.md           # This file
```

## Installation

### Prerequisites

- [Conda](https://docs.conda.io/en/latest/miniconda.html) or [Mamba](https://mamba.readthedocs.io/) package manager

### Setup Environment

1. Clone the repository:
```bash
git clone https://github.com/juanjosecas/psi4me.git
cd psi4me
```

2. Create the conda environment:
```bash
conda env create -f environment.yml
```

3. Activate the environment:
```bash
conda activate psi4me
```

4. Launch JupyterLab:
```bash
jupyter lab
```

## Quick Start

Open the example notebook `notebooks/01_basic_workflow.ipynb` to see demonstrations of:

- Creating molecular structures from SMILES strings
- Generating 3D coordinates with RDKit
- Running quantum chemistry energy calculations with Psi4
- Optimizing molecular geometries
- Comparing energies across multiple molecules
- Visualizing molecular structures and results

## Available Modules

### `molecular_helpers.py`

Functions for working with molecular structures:
- `smiles_to_mol()` - Convert SMILES to RDKit molecule object
- `add_hydrogens()` - Add explicit hydrogens
- `generate_3d_coordinates()` - Generate 3D structure using ETKDG
- `mol_to_xyz()` - Convert molecule to XYZ format

### `quantum_helpers.py`

Functions for quantum chemistry calculations:
- `run_energy_calculation()` - Single point energy calculation
- `optimize_geometry()` - Geometry optimization
- `calculate_frequencies()` - Vibrational frequency analysis

### `visualization_helpers.py`

Functions for plotting and visualization:
- `plot_molecule_2d()` - Display 2D molecular structure
- `plot_energy_comparison()` - Bar chart for energy comparisons
- `plot_frequency_spectrum()` - Vibrational frequency spectrum

## Usage Example

```python
import sys
sys.path.insert(0, '../src')

from molecular_helpers import smiles_to_mol, add_hydrogens, generate_3d_coordinates, mol_to_xyz
from quantum_helpers import run_energy_calculation

# Create water molecule
mol = smiles_to_mol("O")
mol = add_hydrogens(mol)
mol = generate_3d_coordinates(mol)

# Calculate energy
xyz_string = mol_to_xyz(mol)
energy = run_energy_calculation(xyz_string, method="scf", basis="sto-3g")
print(f"Energy: {energy:.6f} Hartrees")
```

## Environment

The environment includes:
- **Python 3.11** - Base programming language
- **Psi4** - Quantum chemistry calculations
- **RDKit** - Molecular structure handling and cheminformatics
- **NumPy** - Numerical computations
- **Pandas** - Data manipulation and analysis
- **Matplotlib** - Plotting and visualization
- **JupyterLab** - Interactive notebook environment

## Contributing

This is a notebook-driven research repository. Feel free to:
- Add new notebooks with different workflows
- Extend helper functions in the `src/` module
- Add molecular data to the `data/` directory
- Save figures to the `figures/` directory

## License

See [LICENSE](LICENSE) file for details.

## Notes

- This repository is designed for research and educational purposes
- No PyPI packaging or distribution is intended
- Helper functions in `src/` should be general-purpose and reusable
- Keep notebooks focused on specific workflows or analyses