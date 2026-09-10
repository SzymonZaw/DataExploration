# Z6 GSE67462 temporal-module analysis protocol

## Purpose

Resolve the multimodal mechanistic core into reproducible temporal modules before pathway-level interpretation.

## Scientific question

Does the multimodal core contain distinct temporal programs corresponding to different phases of the GSE67462 transition, rather than one undifferentiated expression program?

## Procedure

1. Restrict analysis to genes classified as `MULTIMODAL_CORE`.
2. Use the eight timed GSE67462 observations (0, 24, 72, 120, 168, 264, 360, 432 h), averaging replicate branches at each time point.
3. Standardize each gene trajectory across time.
4. Cluster trajectories using deterministic hierarchical clustering with correlation distance.
5. Report module size, centroid trajectory, within-module coherence and dominant active regulatory modalities.
6. Repeat clustering over a small predefined range of module counts and report adjusted-Rand agreement between solutions.
7. Compare module centroids with H3K27ac, H3K4me3, RNAPII, OCT4 and H3K27me3 trajectories.
8. Treat temporal lead as hypothesis-generating, not causal.

## Interpretation boundary

A stable temporal module is a representation-level structure, not a biological mechanism. Pathway interpretation should follow module stability and use an explicitly versioned gene-set file.

## Frozen constraints

- No Z6 predictive-support threshold is changed.
- No model is retrained.
- No cross-dataset result is pooled into GSE67462 modules.
- No causal claim is made from temporal ordering alone.
