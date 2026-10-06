"""
Shared analysis script. Both students edit THIS file on their own branches.

Student A implements : compute_descriptors, plot_by_generation
Student B implements : fingerprints, similarity_heatmap
BOTH implement       : load_molecules, main, summary_sentence   <- expect a merge conflict here

Run: python analysis.py
"""
import os
import matplotlib
import matplotlib.py as plt
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import seaborn as sns
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, AllChem
import yaml


CONFIG = yaml.safe_load(open("config.yaml"))





# ---------- BOTH ----------
def load_molecules(path):
    """Return the DataFrame with an extra column 'mol' holding RDKit Mol objects."""
    df = pd.read_csv(path)

    # Clean headers and locate SMILES column flexibly
    df.columns = df.columns.str.strip()
    smiles_col = next(
        (col for col in df.columns if col.lower() == "smiles"), None
    )

    if smiles_col is None:
        raise KeyError(
            f"Could not find a SMILES column in '{path}'. Found columns: {list(df.columns)}"
        )

    df["mol"] = df[smiles_col].apply(
        lambda s: Chem.MolFromSmiles(str(s)) if pd.notna(s) else None
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


def plot_by_generation(df, output_path="results/properties.png"):
    """One box plot per descriptor, grouped by generation (a 1x5 panel is fine)."""
    """Creates a 1x5 panel plot of box plots for each descriptor grouped by generation."""
    dir_name = os.path.dirname(output_path)
    if dir_name:
        os.makedirs(dir_name, exist_ok=True)

    # Clean missing generation entries and cast generation to int
    plot_df = df.dropna(subset=["generation"]).copy()
    plot_df["generation"] = plot_df["generation"].astype(int)

    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(1, 5, figsize=(18, 4))

    descriptors = [
        ("MW", "MW (Da)"),
        ("logP", "logP"),
        ("TPSA", "TPSA (Å²)"),
        ("HBD", "HBD"),
        ("HBA", "HBA"),
    ]

    # Explicit dictionary mapping prevents color cycling warnings across generations
    palette = {1: "#3498db", 2: "#e74c3c", 3: "#2ecc71"}

    for idx, (col, title) in enumerate(descriptors):
        ax = axes[idx]

        # Updated Seaborn syntax (fixes FutureWarning & UserWarning)
        sns.boxplot(
            x="generation",
            y=col,
            data=plot_df,
            ax=ax,
            hue="generation",
            palette=palette,
            legend=False,
            width=0.4,
            boxprops=dict(alpha=0.7),
        )

        sns.stripplot(
            x="generation",
            y=col,
            data=plot_df,
            ax=ax,
            color="black",
            size=4,
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
    most_polar_gen = int(df.groupby("generation")["TPSA"].mean().idxmax())

    # Extract pair with lowest distance / highest similarity
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
    sim = {(1, 2): 0.35, (1, 3): 0.15, (2, 3): 0.28}
    # After the merge: both, then print(summary_sentence(df, sim))

    print(summary_sentence(df, sim))
    


if __name__ == "__main__":
    main()
