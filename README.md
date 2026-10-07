# Group 4 - Chemoinformatics: three generations of EGFR inhibitors

**Research question:** How do first-, second- and third-generation EGFR kinase inhibitors differ in their physicochemical properties?

**Data:** `data/egfr_inhibitors.csv` - 21 EGFR inhibitors (name, generation 1/2/3/other, SMILES, formula, molecular weight). Structures were validated against their published molecular formula.

Everything happens in **one file, `analysis.py`**. Each function is a stub that raises
`NotImplementedError`. The file header says who implements what.

| Owner | Functions |
|---|---|
| Student A | `compute_descriptors`, `plot_by_generation` |
| Student B | `fingerprints`, `similarity_heatmap` |
| **Both** (this is where the merge conflict happens) | `load_molecules`, `main`, `summary_sentence` |

Shared parameters live in `config.yaml`. Both students must set the value marked
`# BOTH` — you will disagree, and Git cannot decide for you.

Run with `python analysis.py`. It must run without error after your merge.

## Results (fill in after the merge)

https://github.com/gnarcissefhnw/egfr_inhibitors_chemoinformatics
Descriptors: The results summary are as is
    MW: EGFRs by generation have higher molecular weights as the generation increases.
    TPSA: third generation have higher molecular weight and topological polar surface area
    logP: second generation EGFRs are the most lipophilic
    HBD: first generation EFGRs have lower hydrogen bond donor properties
    HBA: third generation EGFRs have higher hydrogen bond acceptor properties
Similarity: Almonertinib and Osimertinib are the most similar pair among the studied EGFR kinase inhibitors.

## Reflection

1. What caused each merge conflict?
    Phase 3 conflict: Simultaneous modifications to shared entrypoints (load_data, main()) and config.yaml on two separate feature branches.
    Phase 4 rejected push: Direct concurrent commits to the same file (results/summary.md) on main without pulling the latest remote commits first.
2. How could branching strategy or file layout have avoided it?
    Keep modular code structures where shared utilities are in separate modules (e.g., utils.py), communicate file ownership before editing shared configs, and work on isolated feature branches rather than pushing directly to main.
3. What is the difference between the history produced by `git pull` and `git pull --rebase`?
    git pull is effectively a combination of two commands: git fetch followed by git merge. In contrast, git pull --rebase runs git fetch followed by git rebase.
    git merge preserves true chronological history by creating a dedicated merge commit, showing where branches diverged and re-joined.
    git rebase rewrites linear commit history by replaying local commits on top of the upstream branch, keeping the log cleaner at the expense of original commit timestamps.
