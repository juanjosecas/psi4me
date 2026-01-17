"""
Molecular structure helpers for RDKit operations
"""

def smiles_to_mol(smiles):
    """
    Convert SMILES string to RDKit molecule object
    
    Parameters
    ----------
    smiles : str
        SMILES string representation of molecule
    
    Returns
    -------
    rdkit.Chem.Mol
        RDKit molecule object
    """
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    return mol


def add_hydrogens(mol):
    """
    Add explicit hydrogens to molecule
    
    Parameters
    ----------
    mol : rdkit.Chem.Mol
        RDKit molecule object
    
    Returns
    -------
    rdkit.Chem.Mol
        Molecule with explicit hydrogens
    """
    from rdkit import Chem
    return Chem.AddHs(mol)


def generate_3d_coordinates(mol, random_seed=42):
    """
    Generate 3D coordinates for molecule using ETKDG method
    
    Parameters
    ----------
    mol : rdkit.Chem.Mol
        RDKit molecule object
    random_seed : int, optional
        Random seed for reproducibility (default: 42)
    
    Returns
    -------
    rdkit.Chem.Mol
        Molecule with 3D coordinates
    
    Raises
    ------
    RuntimeError
        If embedding or optimization fails
    """
    from rdkit.Chem import AllChem
    
    result = AllChem.EmbedMolecule(mol, randomSeed=random_seed)
    if result != 0:
        raise RuntimeError("Failed to generate 3D coordinates for molecule")
    
    opt_result = AllChem.MMFFOptimizeMolecule(mol)
    if opt_result != 0:
        raise RuntimeError("Failed to optimize molecule with MMFF")
    
    return mol


def mol_to_xyz(mol):
    """
    Convert RDKit molecule to XYZ format string
    
    Parameters
    ----------
    mol : rdkit.Chem.Mol
        RDKit molecule object with 3D coordinates
    
    Returns
    -------
    str
        XYZ format string
    
    Raises
    ------
    ValueError
        If molecule has no conformers (no 3D coordinates)
    """
    from rdkit import Chem
    
    if mol.GetNumConformers() == 0:
        raise ValueError("Molecule has no 3D coordinates. Call generate_3d_coordinates() first.")
    
    conf = mol.GetConformer()
    xyz_lines = [str(mol.GetNumAtoms()), ""]
    
    for i, atom in enumerate(mol.GetAtoms()):
        pos = conf.GetAtomPosition(i)
        symbol = atom.GetSymbol()
        xyz_lines.append(f"{symbol:2s} {pos.x:12.6f} {pos.y:12.6f} {pos.z:12.6f}")
    
    return "\n".join(xyz_lines)
