"""
Quantum chemistry calculation helpers for Psi4
"""

def run_energy_calculation(xyz_string, method="scf", basis="sto-3g"):
    """
    Run single point energy calculation with Psi4
    
    Parameters
    ----------
    xyz_string : str
        Molecular geometry in XYZ format
    method : str, optional
        Quantum chemistry method (default: "scf")
    basis : str, optional
        Basis set (default: "sto-3g")
    
    Returns
    -------
    float
        Calculated energy in Hartrees
    """
    import psi4
    
    psi4.core.clean()
    psi4.set_memory('500 MB')
    psi4.set_options({'reference': 'rhf'})
    
    # Parse XYZ string
    lines = xyz_string.strip().split('\n')
    num_atoms = int(lines[0])
    geom_lines = lines[2:2+num_atoms]
    
    # Create Psi4 geometry string
    geom_str = "0 1\n" + "\n".join(geom_lines)
    mol = psi4.geometry(geom_str)
    
    energy = psi4.energy(f"{method}/{basis}", molecule=mol)
    
    return energy


def optimize_geometry(xyz_string, method="scf", basis="sto-3g"):
    """
    Optimize molecular geometry with Psi4
    
    Parameters
    ----------
    xyz_string : str
        Initial molecular geometry in XYZ format
    method : str, optional
        Quantum chemistry method (default: "scf")
    basis : str, optional
        Basis set (default: "sto-3g")
    
    Returns
    -------
    tuple
        (optimized_energy, optimized_xyz_string)
    """
    import psi4
    
    psi4.core.clean()
    psi4.set_memory('500 MB')
    psi4.set_options({'reference': 'rhf'})
    
    # Parse XYZ string
    lines = xyz_string.strip().split('\n')
    num_atoms = int(lines[0])
    geom_lines = lines[2:2+num_atoms]
    
    # Create Psi4 geometry string
    geom_str = "0 1\n" + "\n".join(geom_lines)
    mol = psi4.geometry(geom_str)
    
    energy = psi4.optimize(f"{method}/{basis}", molecule=mol)
    
    # Get optimized geometry
    opt_geom = mol.save_string_xyz()
    
    return energy, opt_geom


def calculate_frequencies(xyz_string, method="scf", basis="sto-3g"):
    """
    Calculate vibrational frequencies with Psi4
    
    Parameters
    ----------
    xyz_string : str
        Molecular geometry in XYZ format
    method : str, optional
        Quantum chemistry method (default: "scf")
    basis : str, optional
        Basis set (default: "sto-3g")
    
    Returns
    -------
    dict
        Dictionary containing energy and frequency information
    """
    import psi4
    
    psi4.core.clean()
    psi4.set_memory('500 MB')
    psi4.set_options({'reference': 'rhf'})
    
    # Parse XYZ string
    lines = xyz_string.strip().split('\n')
    num_atoms = int(lines[0])
    geom_lines = lines[2:2+num_atoms]
    
    # Create Psi4 geometry string
    geom_str = "0 1\n" + "\n".join(geom_lines)
    mol = psi4.geometry(geom_str)
    
    energy, wfn = psi4.frequency(f"{method}/{basis}", molecule=mol, return_wfn=True)
    
    freqs = wfn.frequencies().to_array()
    
    return {
        'energy': energy,
        'frequencies': freqs,
        'num_frequencies': len(freqs)
    }
