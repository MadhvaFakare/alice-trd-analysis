# TRD macros

Every macro here reads one of two inputs, and the distinction matters:

- **`AnalysisResults.root`** — the histograms the PIDStudy task writes. Entries are
  per *tracklet*, not per track, because the charge histograms are filled inside a
  loop over six layers. No per-track momentum is available.
- **`trdpid_tree.root`** — the flat tree. One row per track, carrying `fQ00..fQ25`,
  `fPt` and `fSpecies` together. Anything that needs momentum and per-layer charge
  at the same time has to use this.

All macros save their figures to `trd/plots/` by absolute path, so it does not
matter which directory ROOT is launched from.

## What is in here

| Macro | Input | Purpose |
|---|---|---|
| `trd_io.h` | — | Shared plumbing: opens the flat tree, chains every `DF_*` directory, hands back one `Track` at a time. Included by the `lq_*` macros; no physics in it. |
| `inspect_trd_output.C` | histograms | Lists everything the task wrote, with class, dimensionality and entries. Run this before assuming `Q0_El`/`Q0_Pi` are all that exists. |
| `count_trdpid.C` | flat tree | Sums entries across every `DF_*/trdpid`. Sample yields after quality cuts. |
| `explore_tree.C` | flat tree | Per-track walk through the tree: totals and electron/pion/untagged counts. |
| `trd_q_simple.C` | histograms | Overlays Q0, Q1 and Q2 for electrons and pions, and prints the means. |
| `trd_q0_comparison.C` | histograms | Q0 only, electrons versus pions. |
| `trd_q_slices.C` | histograms | Runs the same electron-efficiency scan on each slice and prints which one separates best. Q1 should win; it covers the TR peak. |
| `trd_fig55.C` | flat tree | Figure 5.5: Q0, Q1, Q2 in layer 0 for both species inside a momentum window. Needs the tree, since the histograms carry no momentum. |
| `sumQ_layer0.C` | flat tree | Per-track `Q0+Q1+Q2` in layer 0, with an optional pT window. Header records the expected yields as a sanity check. |
| `trd_q_separation.C` | flat tree | All three slices across all six layers: separation and pion-rejection tables, plus the per-track average. |
| `trd_likelihood.C` | flat tree | LQ1D and LQ3D as in Schlepper Sec. 5.3: reference distributions, LUTs, combined likelihood. |
| `trd_likelihood_felix.C` | flat tree | Same method with histogram PDFs plus a pseudo-count and 2-fold cross validation, chosen because the electron sample is only ~150 tracks. |
| `lq_step1.C` | flat tree | Step 1 of the LQ1D rewrite on top of `trd_io.h`: count electrons and pions in 1 < pT < 2 GeV/c. |
| `trd_lq1d.C` | flat tree | LQ1D rewrite, in progress. Currently includes only. |
| `trd_all_q.C` | flat tree | Opens the tree and reports the chain size. Output captured in `trd/all_q.txt`. |
| `FORM.C` | — | `TString::Form` syntax scratch. |

## Running

```bash
root -l -q trd/macros/trd_q_slices.C
root -l -q 'trd/macros/sumQ_layer0.C("trdpid_tree.root", 1.0, 2.0)'
```

## `archive/`

Numbered variants of the macros above — `trd_separation01.C`, `trd_separation02.C`,
`trd_separation03.C` and so on. These are practice iterations, where a higher number
is a later and better-understood version of the same exercise. They are kept because
the progression is the record of how the method was learned, but the unnumbered macro
of each name is the one to use.
