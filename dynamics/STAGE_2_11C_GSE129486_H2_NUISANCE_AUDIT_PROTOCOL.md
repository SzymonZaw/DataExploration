# Stage 2.11C — GSE129486 H2 Nuisance-Axis Audit

## Status

Prospective diagnostic gate. This audit does **not** perform Yamanaka comparison, H3 similarity, trajectory-agreement testing, or H3 decision-making.

## Purpose

GSE129486 was selected as an inflammatory temporal control because it has strong temporal replication and representation compatibility. The structural audit identified limited orthogonality to the predeclared H2 nuisance hypothesis. This audit quantifies whether the dataset actually engages the nuisance axes before any H3 similarity is considered.

## Predeclared H2 axes

1. JAK-STAT;
2. NF-kB;
3. interferon/STAT-like inflammatory response.

These axes were specified from the historical signal-attribution analysis. The present audit must not redefine them after inspecting GSE129486.

## Questions

- Are JAK-STAT and NF-kB activities detectable in GSE129486?
- How do those activities vary with time, stimulation and cell line?
- Is the inflammatory signal broadly shared across blocks or concentrated in particular conditions?
- Does the dataset provide evidence of nuisance-axis engagement strong enough that it should not be treated as an orthogonal H3 control?

## Analysis lock

The audit may describe pathway activity distributions and stratified means. It must not use Yamanaka similarity to choose thresholds, exclude samples, or redefine the control role.

No biological conclusion is attached to a pathway activity score alone. Activity is treated as evidence of axis engagement, not proof that the same axis caused the Yamanaka signal.

## Mapping

Use the same Ensembl-to-HGNC namespace normalization already locked for GSE129486. Duplicate canonical HGNC symbols are averaged before network scoring. Persist input hashes and mapping counts.

## Output

Persist `results/Dynamics/stage2_11c_h3_control_redesign/H3_GSE129486_H2_NUISANCE_AUDIT.json`.

## Interpretation

- `RISK_CONFIRMED`: nuisance-axis engagement is sufficiently clear that GSE129486 cannot serve as an orthogonal H3 control without an explicit H2-aware design.
- `RISK_LIMITED`: some engagement is present, but the control may retain a restricted diagnostic role; this still does not make it a clean H3 control.
- `RISK_NOT_EVIDENT`: the predeclared axes show little engagement; orthogonality risk is reduced, but H3 eligibility still requires the independent temporal-control requirements.
- `NOT_AVAILABLE`: pathway scoring could not be executed reproducibly; no H3 decision is permitted.

The final category must be based on a threshold committed before interpreting the resulting activity values. Until such a threshold is separately locked, the machine-readable result remains descriptive and carries no H3 decision authority.
