# Z6 GSE67462 replicate-level directional falsification

## Purpose

Diagnostic follow-up to the aggregate directional mechanistic audit. The goal is to test whether the regulatory-to-expression temporal association survives the two independent GSE67462 expression replicates rather than appearing only after averaging them.

GSE67462 contains two expression replicates for each time point from day 0 through day 18, with a separate iPSC endpoint. The paired GSE67520 regulatory profiles provide one temporal regulatory trajectory per modality/time point rather than two independent regulatory replicates.

Therefore this audit tests **replicate stability of the expression response**, not full replicate stability of the regulatory measurement itself.

## Frozen boundary

This audit does not modify:

- Z6 predictive support;
- LODO or within-system validation;
- temporal module assignments;
- multimodal support decisions;
- mechanistic falsification decisions already recorded.

## Primary test

For each module M1, M4 and M6 and each active modality:

1. retain the common regulatory trajectory;
2. use expression replicate 1 as the response trajectory;
3. use expression replicate 2 as the response trajectory;
4. compute concurrent Spearman rho;
5. compute one-step lead rho, `rho(regulator[t], expression[t+1])`;
6. compute `lead_gain = lead_rho - concurrent_rho`;
7. compare the observed median gene-level lead gain against a circular time-shift null separately for each expression replicate.

A directional effect is considered **replicate-concordant** only when both expression replicates have the same sign of median lead gain. Stronger evidence additionally requires both replicates to exceed their own null q95; this is deliberately conservative.

## Chain analysis

The data do not contain independent regulatory replicates for GSE67520. Consequently the audit must not claim replication of a complete chain such as:

`OCT4 -> H3K27ac -> expression`

Instead it can report the aggregate regulatory ordering separately and test whether its **expression endpoint** is stable across the two GSE67462 replicates.

## Interpretation

- Same-sign replicate effects strengthen reproducibility of the temporal association.
- Opposite signs falsify a stable replicate-level directional interpretation.
- Same sign without null exceedance is weak/ambiguous evidence.
- Even concordant replicate effects remain observational and do not establish causality.
- A negative result is especially informative because the aggregate audit already found no active modality above its circular null.
