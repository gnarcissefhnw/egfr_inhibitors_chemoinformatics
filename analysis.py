"""
Shared analysis script. Both students edit THIS file on their own branches.

Student A implements : compute_descriptors, plot_by_generation
Student B implements : fingerprints, similarity_heatmap
BOTH implement       : load_molecules, main, summary_sentence   <- expect a merge conflict here

Run: python analysis.py
"""
import yaml
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CONFIG = yaml.safe_load(open("config.yaml"))

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, AllChem

# ---------- BOTH ----------
def load_molecules(path):
    """Return the DataFrame with an extra column 'mol' holding RDKit Mol objects."""
    raise NotImplementedError


# ---------- Student A ----------
def compute_descriptors(df):
    """Add columns MW, logP, TPSA, HBD, HBA (rdkit.Chem.Descriptors). Return df."""
    raise NotImplementedError


def plot_by_generation(df, out="results/properties.png"):
    """One box plot per descriptor, grouped by generation (a 1x5 panel is fine)."""
    raise NotImplementedError


# ---------- Student B ----------
def fingerprints(df, radius):
    """Return a list of Morgan fingerprints (2048 bits) in the row order of df."""
    raise NotImplementedError


def similarity_heatmap(fps, names, out="results/similarity.png"):
    """Tanimoto similarity matrix as a heatmap with the names on both axes.
    Return the matrix (numpy array)."""
    raise NotImplementedError


# ---------- BOTH ----------
def summary_sentence(df, sim):
    """One sentence: which generation is heaviest / most polar, and the most similar pair."""
    raise NotImplementedError


def main():
    df = load_molecules(CONFIG["compounds"])
    # Student A: descriptors + box plots
    # Student B: fingerprints + heatmap
    # After the merge: both, then print(summary_sentence(df, sim))
    raise NotImplementedError


if __name__ == "__main__":
    main()
