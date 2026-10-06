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

<!-- one sentence answering the research question, one figure -->

## Reflection

1. What caused each merge conflict?
2. How could branching strategy or file layout have avoided it?
3. What is the difference between the history produced by `git pull` and `git pull --rebase`?
