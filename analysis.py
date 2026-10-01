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
    df["mol"] = df["SMILES"].apply(
        lambda s: Chem.MolFromSmiles(s) if pd.notna(s) else None
    )
    return df


# ---------- Student A ----------
def compute_descriptors(df):
    """Add columns MW, logP, TPSA, HBD, HBA (rdkit.Chem.Descriptors). Return df."""
    df["MW"] = df["mol"].apply(
        lambda m: Descriptors.MolWt(m) if m is not None else None
    )
    df["logP"] = df["mol"].apply(
        lambda m: Descriptors.MolLogP(m) if m is not None else None
    )
    df["TPSA"] = df["mol"].apply(
        lambda m: Descriptors.TPSA(m) if m is not None else None
    )
    df["HBD"] = df["mol"].apply(
        lambda m: Descriptors.NumHDonors(m) if m is not None else None
    )
    df["HBA"] = df["mol"].apply(
        lambda m: Descriptors.NumHAcceptors(m) if m is not None else None
    )

    return df


def plot_by_generation(df, out="results/properties.png"):
    """One box plot per descriptor, grouped by generation (a 1x5 panel is fine)."""
    """Creates a 1x5 panel plot of box plots for each descriptor grouped by generation."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(1, 5, figsize=(18, 4))

    descriptors = [
        ("MW", "Molecular Weight (Da)"),
        ("logP", "Lipophilicity (logP)"),
        ("TPSA", "TPSA (Å²)"),
        ("HBD", "H-Bond Donors"),
        ("HBA", "H-Bond Acceptors"),
    ]

    palette = ["#3498db", "#e74c3c", "#2ecc71"]

    for idx, (col, title) in enumerate(descriptors):
        ax = axes[idx]

        # Draw Boxplot
        sns.boxplot(
            x="generation",
            y=col,
            data=df,
            ax=ax,
            palette=palette,
            width=0.4,
            boxprops=dict(alpha=0.7),
        )

        # Overlay Individual Data Points
        sns.stripplot(
            x="generation",
            y=col,
            data=df,
            ax=ax,
            color="black",
            size=5,
            jitter=0.15,
        )

        ax.set_title(title, fontsize=11, fontweight="bold", pad=8)
        ax.set_xlabel("Generation", fontsize=10)
        ax.set_ylabel(col, fontsize=10)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Plot saved successfully to {output_path}")


if __name__ == "__main__":
    # Execution workflow
    df = load_molecules("egfr_inhibitors.csv")
    df = compute_descriptors(df)
    plot_by_generation(df)


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
    # Find generation with highest mean Molecular Weight
    heaviest_gen = int(df.groupby("generation")["MW"].mean().idxmax())

    # Find generation with highest mean Polar Surface Area (polarity)
    most_polar_gen = int(df.groupby("generation")["TPSA"].mean().idxmax())

    # Extract the most similar pair from the similarity structure 'sim'
    # Expects sim to be a dict mapping tuple pairs (e.g., (1, 3)) to distance/similarity
    most_sim_pair = min(sim, key=sim.get)

    return (
        f"The {heaviest_gen}rd generation is the heaviest, the {most_polar_gen}nd generation is the most polar, "
        f"and generations {most_sim_pair[0]} and {most_sim_pair[1]} form the most similar pair."
    )


def main():
    df = load_molecules(CONFIG["compounds"])
    # Student A: descriptors + box plots
    df = compute_descriptors(df)
    plot_by_generation(df, os.path.join(CONFIG["results_dir"], "properties.png"))
    # Student B: fingerprints + heatmap
    # After the merge: both, then print(summary_sentence(df, sim))
    print(summary_sentence(df, sim))
    raise NotImplementedError


if __name__ == "__main__":
    main()
