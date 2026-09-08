# Stage 2.11C — H3 Control Redesign Manifest V2

Status: **prospective — supersedes the final-control section of V1 before any H3 similarity analysis**

## Purpose
Test whether the representation pipeline reproduces Yamanaka-like temporal structure in biologically unrelated temporal processes. H3 is a negative control for generic temporal/representation structure, not a control for shared biology.

## Current design decision
The previous H3 candidate GSE3945 is **rejected for execution** because the technical audit found only ~0.78% probe-to-HGNC mapping, ~1% PROGENy overlap and ~0.86% DoRothEA overlap, plus mostly singleton timepoints. It therefore cannot provide a meaningful H3 test under the locked representation contract.

GSE129486 is **reclassified as H2 inflammatory nuisance control**, not a clean H3 control. Its temporal structure is strong, but the nuisance-axis audit shows marked JAK-STAT and NF-κB activity, making it non-orthogonal to the predeclared H2 hypothesis.

## New H3 candidate
**GSE263713 — human primary skin fibroblast circadian RNA-seq time course.**

NCBI GEO reports six healthy individuals measured at 13 time points over 48 hours (0–48 h, every 4 h), with RNA-seq expression data and repeated temporal trajectories. This gives an explicit repeated-measures unit (individual fibroblast line) and a process class that is orthogonal to cellular reprogramming and to the predeclared inflammatory H2 axes.

Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263713

Expected processed input:
`Data/GSE263713_raw_counts.tsv.gz`

Expected acquisition URL:
`https://ftp.ncbi.nlm.nih.gov/geo/series/GSE263nnn/GSE263713/suppl/GSE263713_raw_counts.tsv.gz`

The exact SHA-256 is intentionally **not locked yet**. It must be recorded from the locally acquired file before technical validation is accepted. No similarity statistic may be computed before that lock.

## H3 eligibility criteria
Before H3 similarity:

- 6 repeated-measures individuals must be recoverable from metadata/sample identifiers;
- 13 expected timepoints must be recoverable: 0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48 h;
- each individual must have a complete temporal course;
- matrix must be gene × sample and predominantly numeric;
- identifier namespace must be mapped reproducibly to HGNC;
- mapping/collision statistics must be persisted;
- representation/network overlap with the existing Yamanaka contract must be checked before H3 similarity;
- repeated-measures permutation/exchangeability unit must be the individual trajectory, never individual samples independently;
- H3 similarity must remain locked to the prospective `delta_H3=0.10`, one-sided `p<0.05`, BH `q<0.05` framework.

## H2 control
GSE129486 remains available as a nuisance-axis control. Its actual JAK-STAT/NF-κB activity is audited independently of H3 similarity. It must not be used as the sole orthogonal H3 control.

## Remaining requirement
The final H3 design still requires a second independent H3 process class. GSE263713 is the first strong candidate; a second orthogonal temporal process must be selected and technically audited before the H3 similarity run.

## Guardrails
- Candidate selection must not inspect Yamanaka similarity.
- Technical failure returns `H3_UNRESOLVED`; it does not trigger replacement after looking at similarity.
- H3 support does not prove Yamanaka specificity or causality.
- H3 failure does not identify a particular confounder by itself.
