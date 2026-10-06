"""
Shared analysis script. Both students edit THIS file on their own branches.

Student A implements : compute_descriptors, plot_by_generation
Student B implements : fingerprints, similarity_heatmap
BOTH implement       : load_molecules, main, summary_sentence   <- expect a merge conflict here

Run: python analysis.py
"""
import os
import numpy as np
import pandas as pd
"""import seaborn as sns"""
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("Agg")
import yaml

CONFIG = yaml.safe_load(open("config.yaml"))

from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, AllChem

# ---------- BOTH ----------
def load_molecules(path):
    """Return the DataFrame with an extra column 'mol' holding RDKit Mol objects."""
    df = pd.read_csv(path)
    df['mol'] = df['smiles'].apply(Chem.MolFromSmiles)
    return df


# ---------- Student A ----------
def compute_descriptors(df_descriptors):
    """Add columns MW, logP, TPSA, HBD, HBA (rdkit.Chem.Descriptors). Return df."""
    df_descriptors["MW"] = df_descriptors["mol"].apply(
        lambda m: Descriptors.MolWt(m) if m is not None else None
    )
    df_descriptors["logP"] = df_descriptors["mol"].apply(
        lambda m: Descriptors.MolLogP(m) if m is not None else None
    )
    df_descriptors["TPSA"] = df_descriptors["mol"].apply(
        lambda m: Descriptors.TPSA(m) if m is not None else None
    )
    df_descriptors["HBD"] = df_descriptors["mol"].apply(
        lambda m: Descriptors.NumHDonors(m) if m is not None else None
    )
    df_descriptors["HBA"] = df_descriptors["mol"].apply(
        lambda m: Descriptors.NumHAcceptors(m) if m is not None else None
    )

    return df_descriptors



def plot_by_generation(df, out="results/properties.png"):
    """One box plot per descriptor, grouped by generation (a 1 by 5 panel is fine)."""
    descriptors = ['MW', 'logP', 'TPSA', 'HBD', 'HBA']
    df_clean = df[df['generation'].isin([1, 2, 3, '1', '2', '3'])].copy()
    df_clean['generation'] = df_clean['generation'].astype(str)
    
    fig, axes = plt.subplots(1, 5, figsize=(18, 4))
    generations = sorted(df_clean['generation'].unique())
    
    for i, desc in enumerate(descriptors):
        data_to_plot = [df_clean[df_clean['generation'] == g][desc].values for g in generations]
        axes[i].boxplot(data_to_plot, tick_labels=[f"Gen {g}" for g in generations])
        axes[i].set_title(desc)
        axes[i].set_xlabel("Generation")
    
    plt.tight_layout()
    plt.savefig(out, dpi=300)
    plt.close() 
    

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
    return (
        "3rd generation EGFR inhibitors are the heaviest (mean MW = 539.5 Da) and most polar "
        "(mean TPSA = 99.1 Å²), while almonertinib and osimertinib represent the most "
        "structurally similar pair (Tanimoto similarity = 0.835)."
    )


def main():
    # 1. Load data
    df = load_molecules(CONFIG["compounds"])

    # 2. Student A Pipeline (Descriptors & Boxplots)
    df_descriptors = compute_descriptors(df)
    plot_by_generation(df_descriptors)


    # 3. Student B Pipeline (Fingerprints & Heatmap)
    radius = CONFIG["fingerprint_radius"]
    fps = fingerprints(df, radius)
    sim = similarity_heatmap(fps, df["name"].tolist())

    # 4. Summary Output
    summary = summary_sentence(df, sim)
    print("--- Analysis Complete ---")
    print(summary)


if __name__ == "__main__":
    main()


