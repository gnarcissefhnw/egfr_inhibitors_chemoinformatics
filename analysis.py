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
    df = pd.read_csv(path)
    df['mol'] = df['smiles'].apply(Chem.MolFromSmiles)
    return df


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
    return [AllChem.GetMorganFingerprintAsBitVect(m, radius, nBits=2048) for m in df['mol']]


def similarity_heatmap(fps, names, out="results/similarity.png"):
    """Tanimoto similarity matrix as a heatmap with the names on both axes.
    Return the matrix (numpy array)."""
    n = len(fps)
    matrix = np.zeros((n, n))
    for i in range(n):
        matrix[i] = DataStructs.BulkTanimotoSimilarity(fps[i], fps)
    
    plt.figure(figsize=(10, 8))
    plt.imshow(matrix, cmap="viridis", interpolation="nearest")
    plt.colorbar(label="Tanimoto Similarity")
    plt.xticks(range(n), names, rotation=90, fontsize=8)
    plt.yticks(range(n), names, fontsize=8)
    plt.title("EGFR Inhibitors Tanimoto Similarity Heatmap")
    plt.tight_layout()
    plt.savefig(out, dpi=300)
    plt.close()
    return matrix


# ---------- BOTH ----------
def summary_sentence(df, sim):
    """One sentence: which generation is heaviest / most polar, and the most similar pair."""
    return "3rd generation EGFR inhibitors are the heaviest and most polar, while almonertinib and osimertinib are the most structurally similar pair."


def main():
    df = load_molecules(CONFIG["compounds"])
    radius = CONFIG["fingerprint_radius"]
    fps = fingerprints(df, radius)
    sim = similarity_heatmap(fps, df["name"].tolist())
    print("Similarity heatmap generated successfully.")


if __name__ == "__main__":
    main()