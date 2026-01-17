"""
Visualization helpers for molecular structures and quantum chemistry results
"""

import matplotlib.pyplot as plt
import numpy as np


def plot_molecule_2d(mol, figsize=(6, 6), save_path=None):
    """
    Plot 2D molecular structure using RDKit
    
    Parameters
    ----------
    mol : rdkit.Chem.Mol
        RDKit molecule object
    figsize : tuple, optional
        Figure size (default: (6, 6))
    save_path : str, optional
        Path to save figure (default: None, shows plot)
    
    Returns
    -------
    matplotlib.figure.Figure
        Figure object
    """
    from rdkit.Chem import Draw
    
    img = Draw.MolToImage(mol, size=(600, 600))
    
    fig, ax = plt.subplots(figsize=figsize)
    ax.imshow(img)
    ax.axis('off')
    
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=150)
    
    return fig


def plot_energy_comparison(labels, energies, ylabel="Energy (Hartrees)", 
                          figsize=(10, 6), save_path=None):
    """
    Create bar plot comparing energies
    
    Parameters
    ----------
    labels : list
        Labels for each energy value
    energies : list
        Energy values to plot
    ylabel : str, optional
        Y-axis label (default: "Energy (Hartrees)")
    figsize : tuple, optional
        Figure size (default: (10, 6))
    save_path : str, optional
        Path to save figure (default: None, shows plot)
    
    Returns
    -------
    matplotlib.figure.Figure
        Figure object
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    x = np.arange(len(labels))
    ax.bar(x, energies, color='steelblue', alpha=0.7)
    ax.set_xlabel('Molecule', fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=150)
    
    return fig


def plot_frequency_spectrum(frequencies, figsize=(10, 6), save_path=None):
    """
    Plot vibrational frequency spectrum
    
    Parameters
    ----------
    frequencies : array-like
        Vibrational frequencies in cm^-1
    figsize : tuple, optional
        Figure size (default: (10, 6))
    save_path : str, optional
        Path to save figure (default: None, shows plot)
    
    Returns
    -------
    matplotlib.figure.Figure
        Figure object
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # Filter out imaginary frequencies (negative values)
    real_freqs = [f for f in frequencies if f > 0]
    
    ax.stem(real_freqs, np.ones_like(real_freqs), basefmt=' ')
    ax.set_xlabel('Frequency (cm⁻¹)', fontsize=12)
    ax.set_ylabel('Intensity (arb. units)', fontsize=12)
    ax.set_title('Vibrational Frequency Spectrum', fontsize=14)
    ax.grid(axis='x', alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=150)
    
    return fig
