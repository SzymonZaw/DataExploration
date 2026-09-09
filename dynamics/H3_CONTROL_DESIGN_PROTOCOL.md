# H3 Control Design Protocol (prospective redesign)

## Purpose

H3 asks whether the temporal structure observed in the Yamanaka representation could arise from a generic temporal process rather than reprogramming-specific biology.

This protocol is deliberately stricter than the original GSE263713 comparison. GSE263713 remains an exploratory control, not a decisive H3 falsification control, because its circadian trajectory is oscillatory and therefore geometrically mismatched to the approximately monotonic Yamanaka trajectory.

## 1. Eligibility criteria for a decisive H3 temporal control

A candidate control must satisfy all mandatory criteria before its outcome is inspected:

1. **Biological independence**: the process is not Yamanaka/OSKM reprogramming and does not share an obvious direct reprogramming mechanism.
2. **Temporal geometry**: the primary state transition is approximately monotonic or otherwise directional over the sampled interval; strongly periodic/cyclic processes are excluded from the decisive set.
3. **Time scale**: the experiment spans days, preferably with a duration comparable to the 10-day Yamanaka window. Hour-scale-only controls are exploratory.
4. **Transcriptional magnitude**: the process must produce a measurable transcriptomic state transition; a near-flat temporal process is not an adequate negative control.
5. **Independent trajectories**: biological replicates/donors/lines must be identifiable. Technical replicates must not be counted as independent trajectories.
6. **Temporal coverage**: at least 4 ordered timepoints spanning the transition, with no requirement to force the control onto the Yamanaka day labels.
7. **Representation compatibility**: sufficient gene overlap for the locked PROGENy and DoRothEA representation gates.
8. **H2 separation**: the control must not be selected primarily because it activates the same inflammatory/stress axes already under H2 investigation.
9. **Data integrity**: raw/processed count data and sample metadata must permit deterministic reconstruction of the time course.

## 2. Preferred biological class

The preferred class is an unrelated, directional differentiation or terminal state transition measured by bulk RNA-seq over several days. Examples include adipogenic or other lineage differentiation, provided the specific dataset passes the H2-separation and design audits.

A promising candidate identified during the redesign is **GSE249195**, primary human preadipocyte differentiation. GEO reports 24 bulk RNA-seq samples across six timepoints over 14 days, with four replicates per timepoint, and provides a raw-count supplementary file. This is a candidate only; it is **not pre-approved** until the repository-side representation/design/H2 audits pass.

## 3. Locked structural eligibility thresholds

For GSE249195 and subsequent candidates evaluated under the same audit implementation, the following thresholds are frozen before similarity testing:

- PROGENy overlap >= **0.90**.
- DoRothEA overlap >= **0.90**.
- PC1 absolute Spearman correlation with native time >= **0.80**.
- PC1 adjacent-step monotonic fraction >= **0.80** (at least 4 of 5 adjacent intervals for a six-timepoint design must have the same sign after accounting for arbitrary PC1 sign).
- Fraction of gene trajectories with absolute Spearman |rho| >= **0.80** >= **0.60**.
- Median absolute gene-level Spearman |rho| >= **0.60**.
- Median absolute day-final vs day-zero log2-CPM difference >= **0.50**.
- Total temporal span >= **7 days**.

These are structural eligibility criteria, not outcome thresholds for the later Yamanaka-similarity test. Failure means the candidate is not an eligible decisive H3 control; the similarity analysis is not run.

The H2 nuisance audit is reported separately and reviewed before eligibility is frozen. It is not converted into a post-hoc numerical exclusion rule in the GSE249195 audit script.

## 4. Temporal alignment rule

Do not equate biological hours/days across experiments. For an eligible directional control:

- retain the control's native ordered timepoints;
- normalize each trajectory to relative time [0,1];
- compare temporal shape using a predeclared statistic;
- if interpolation to the four Yamanaka positions is required, use one locked interpolation rule and do not optimize it after seeing the result.

The control must not be made more similar to Yamanaka by selecting favorable phases, windows, or interpolation points.

## 5. Primary statistic

For each feature family (PROGENy and DoRothEA), calculate the median across common features of the feature-wise Pearson correlation between the canonical Yamanaka four-point trajectory and the aligned control trajectory.

The feature-matching permutation null remains a diagnostic null: feature rows are permuted between target and control while preserving each feature's temporal trajectory.

This null does **not** test every possible generic temporal artifact. Therefore a negative result is interpreted only within the predeclared control design.

## 6. Aggregation rule across PROGENy and DoRothEA

PROGENy and DoRothEA are treated as two representation families, **not as independent biological experiments**.

For a decisive control, both families must satisfy the same predeclared directional decision criterion. A disagreement is classified as **INCONCLUSIVE**, not averaged away and not resolved post hoc.

Suggested locked labels:

- `NEGATIVE_CONCORDANT`: both families fail to show evidence of conserved Yamanaka-like temporal shape.
- `POSITIVE_CONCORDANT`: both families show evidence of conserved shape.
- `INCONCLUSIVE`: one family is positive/equivocal while the other is negative.

The exact alpha/empirical-tail threshold must be fixed in the analysis configuration before the candidate's result is inspected.

## 7. Aggregation rule across independent controls

A final H3 closure requires at least **two eligible, biologically orthogonal controls**.

- Two `NEGATIVE_CONCORDANT` controls: H3 receives strong falsification support under this operational definition.
- Any `POSITIVE_CONCORDANT` eligible control: H3 remains open and requires investigation of the shared temporal structure.
- Any `INCONCLUSIVE` control: H3 remains unresolved.
- A control failing eligibility after inspection is not silently converted into a negative control; it is excluded with the reason documented.

This rule is locked before the second decisive control is analyzed.

## 8. GSE263713 disposition

GSE263713 is retained because it is technically valid and useful as an exploratory contrast, but it is not counted toward the two-control decisive H3 gate.

Its observed results remain unchanged:

- PROGENy median trajectory correlation ≈ -0.337, empirical one-sided p ≈ 0.816.
- DoRothEA median trajectory correlation ≈ 0.158, empirical one-sided p ≈ 0.157.

These values should not be described as proof that H3 is false. The current script records the decision as `H3_EXPLORATORY_ONLY`.

## 9. Historical-data protection

This redesign does not modify the historical Stage 2.11C null files and does not rerun `run_yamanaka_control_null.py`.

The canonical Yamanaka trajectory used by the prospective H3 controls remains the pointwise arithmetic mean of the historical `a_values` and `b_values`, after verifying that both contain exactly four points.

## 10. Required audit sequence for the next decisive control

1. Dataset design audit.
2. Sample/trajectory independence audit.
3. Gene-ID/HGNC mapping audit.
4. PROGENy/DoRothEA representation overlap gate.
5. H2/nuisance-axis audit.
6. Directional temporal-geometry audit using the locked thresholds above.
7. Freeze eligibility decision before similarity testing.
8. Run the identical locked similarity/null procedure.
9. Apply the predeclared within-control aggregation rule.
10. Apply the predeclared cross-control H3 decision rule.

No threshold, time window, phase, feature family, or control inclusion decision may be changed after observing the similarity result.
