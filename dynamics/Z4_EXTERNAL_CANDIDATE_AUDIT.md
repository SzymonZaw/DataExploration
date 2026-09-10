# Z4 External Validation Candidate Audit

## Purpose

This document records the **pre-analysis selection of external validation candidates** for Z4. It is deliberately based on dataset provenance and design metadata, not on results obtained by transferring the frozen GSE67462 representation.

The selection question is:

> Which independent public experiment can test whether the GSE67462-derived state/transition signal is biologically specific rather than a consequence of the original mouse secondary-reprogramming context?

No candidate may be promoted or rejected using its eventual transfer score.

## Discovery object

**GSE67462** is the frozen discovery system. NCBI describes it as mouse OSKM-mediated secondary reprogramming with bulk expression profiling, 18 samples, two replicates at each time point, and a matched iPSC endpoint. The experiment also has the regulatory companion GSE67520. This makes GSE67462 suitable for discovery, but not independent validation of itself.

## Candidate ranking

| Candidate | Biological relationship to target | Independence | Temporal design | Replication | Independent endpoint | Main nuisance/context | Z4 role | Pre-analysis status |
|---|---|---|---|---|---|---|---|---|
| **GSE297234** | Human fibroblast partial reprogramming with Sendai OSKM | Strong: different study, donors, species, protocol and technology | D0/D3/D7/D10 | Two donors; not technical replicates | Biological rejuvenation/mesenchymal-drift phenotype is defined independently of this project | species, donor age, scRNA-seq, Sendai, partial rather than full reprogramming | **Primary external target candidate** | **ELIGIBLE WITH GATES** |
| **GSE28688** | Human fibroblast OSKM induction | Strong: independent study, donor, protocol and platform | D0/24h/48h/72h plus iPSC/ESC endpoints | Two biological/technical branches in early time points | iPSC/ESC endpoint | viral transduction, strong early stress/redox/p53 response | **Specificity/context challenge** | **ELIGIBLE AS NEGATIVE/CHALLENGE CONTROL** |
| **GSE67462** | Mouse secondary OSKM reprogramming | None; discovery dataset | D0–D18 | Two replicates/timepoint | iPSC endpoint | original discovery context | Discovery only | **EXCLUDED** |

## Candidate 1 — GSE297234

NCBI describes GSE297234 as 10x Genomics scRNA-seq of human fibroblasts from a young donor (GM23815, 22 years) and an aged donor (GM00731, 96 years), treated with Sendai-virus OSKM for up to 10 days. Samples are available at days 0, 3, 7 and 10. The study explicitly focuses on reversal of mesenchymal drift during partial reprogramming before dedifferentiation/pluripotency acquisition.

Why this is valuable for Z4:

- independent experiment and investigators;
- different species from the discovery system;
- different cell source and donor structure;
- different reprogramming implementation (Sendai OSKM vs secondary inducible mouse system);
- single-cell rather than bulk measurement;
- biological target is not defined merely as "time after OSKM" but includes an independently characterized mesenchymal-drift/rejuvenation phenotype.

Why it is **not automatically sufficient**:

- it is partial reprogramming, not the same endpoint as GSE67462 iPSC acquisition;
- only four time points;
- species and technology differences create a large domain shift;
- donor age is biologically meaningful and must not be absorbed into the target score.

Therefore GSE297234 should be treated as a **homologous external target**, not as a literal replication of the GSE67462 trajectory.

## Candidate 2 — GSE28688

NCBI describes GSE28688 as human foreskin fibroblast HFF1 subjected to OSKM transduction, with measurements at 24 h, 48 h and 72 h, plus iPSC and hESC samples. The study reports that viral transduction itself induces redox perturbation, oxidative damage, p53 activation, senescence and apoptosis.

This makes GSE28688 scientifically useful even though it is not a clean positive validation of the same biological phase:

- it is independent;
- it contains the same broad OSKM intervention class;
- it contains early temporal structure;
- importantly, it exposes a known context/stress component that can mimic or dominate reprogramming-associated expression.

For Z4 it should therefore serve as a **specificity challenge**. A representation that simply recognizes "OSKM + time + stress" should perform positively here. A biologically specific reprogramming representation should show a distinguishable target pattern rather than treating the acute viral/stress response as equivalent to the discovery transition.

It must not be used as a post-hoc negative control selected because its result is unfavorable.

## Why these two candidates are complementary

The strongest public-data Z4 design currently available is therefore asymmetric:

```text
                         FROZEN GSE67462 REPRESENTATION
                                      │
                 ┌────────────────────┴────────────────────┐
                 │                                         │
        GSE297234                                   GSE28688
   homologous target transfer                  context/specificity challenge
   human partial OSKM                           human early OSKM + stress
                 │                                         │
          should score +                         should NOT be treated as
                 │                              equivalent target biology
                 └────────────────────┬────────────────────┘
                                      │
                           Z4 specificity decision
```

A positive result in GSE297234 alone would establish external signal transfer, but not full Z4 specificity. The stronger claim requires discrimination from the GSE28688 context challenge and robustness to donor/protocol effects.

## Commensurability gates before transfer

Before any validation score is inspected, the following gates must pass:

### Gate C1 — target construct

The external dataset must contain a biologically interpretable state/transition homologous to the frozen discovery target. For GSE297234 this is **partial OSKM reprogramming / mesenchymal-drift reversal**, not full iPSC acquisition.

### Gate C2 — molecular feature mapping

A deterministic human-to-mouse orthology mapping must be frozen before validation. One-to-many mappings must not be resolved using validation outcomes.

### Gate C3 — endpoint independence

The target endpoint must not be computed from the same frozen genes used to construct the transferred score. Prefer independently defined phenotype/marker information from the source study.

### Gate C4 — nuisance observability

Donor, age, protocol, time, platform and sequencing/technical metadata must be available sufficiently to test context dependence.

### Gate C5 — control comparability

The specificity challenge must share plausible nuisance structure with the target. A completely unrelated dataset is a weak negative control for this question.

### Gate C6 — no adaptive tuning

If a candidate fails a commensurability gate, it is excluded before looking at transfer results. No thresholds or mappings may be changed to rescue it.

## Frozen candidate decision

As of this audit:

- **GSE297234 = primary external target candidate**;
- **GSE28688 = predefined context/specificity challenge**;
- **GSE67462 = discovery only**.

The next implementation should therefore perform a **commensurability audit only**. It should not yet calculate whether the frozen Z6 score transfers successfully.

## Scientific interpretation boundary

If GSE297234 passes the commensurability gates and the frozen representation transfers, while GSE28688 does not receive an equivalent target score after accounting for its acute stress/redox context, that would provide meaningful evidence for Z4 specificity.

If both datasets score positively, the likely interpretation is generic OSKM/time/stress/context sensitivity rather than biological specificity.

If GSE297234 fails commensurability, Z4 remains unresolved rather than being declared negative.

## Provenance sources

- NCBI GEO GSE67462: mouse secondary OSKM reprogramming, bulk expression, two replicates per time point, matched iPSC endpoint.
- NCBI GEO GSE297234: human fibroblasts, young/aged donors, Sendai OSKM, D0/D3/D7/D10 scRNA-seq.
- NCBI GEO GSE28688: human HFF1 fibroblasts, OSKM transduction at 24/48/72 h plus iPSC/ESC endpoints; source study reports viral-transduction-associated redox/p53/stress effects.

These provenance facts are used only for candidate selection. No transfer result is incorporated into this document.
