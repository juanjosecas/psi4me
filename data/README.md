# Data Directory

This directory is for storing:
- Input molecular structures (XYZ, SDF, MOL files)
- SMILES datasets
- Calculation results (energies, geometries)
- Any other input data for notebooks

## Example Data Files

Add your molecular data files here. They can be referenced from notebooks using relative paths:

```python
import pandas as pd

# Read molecular data
df = pd.read_csv('../data/molecules.csv')
```

## Recommended Format

For sharing molecular structures, common formats include:
- CSV files with SMILES strings and metadata
- XYZ files for 3D coordinates
- SDF/MOL files for more complex structures
