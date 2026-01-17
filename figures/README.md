# Figures Directory

This directory is for storing output figures and plots generated from notebooks.

## Usage

Save figures from your notebooks to this directory:

```python
import matplotlib.pyplot as plt
from visualization_helpers import plot_molecule_2d

# Create and save figure
fig = plot_molecule_2d(mol, save_path='../figures/water_molecule.png')
```

## Organization

Consider organizing figures by:
- Notebook name (e.g., `01_basic_workflow_energy_comparison.png`)
- Analysis type (e.g., `energy_comparison_small_molecules.png`)
- Date (e.g., `2024-01-15_optimization_results.png`)

This keeps your figures organized and easy to reference in reports or presentations.
