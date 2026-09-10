# Z6 gene-set inputs

These files are external reference resources and should be downloaded locally rather than committed to the repository unless licensing permits.

Recommended source for the first enrichment pass: **Mouse MSigDB 2026.1.Mm**, using MGI gene-symbol GMT files. The collection page provides separate GMT downloads for GO Biological Process and the Reactome subset of canonical pathways, plus Hallmark sets useful for nuisance/context auditing.

Expected local filenames:

- `GO_BP_mouse.gmt`
- `Reactome_mouse.gmt`
- `Hallmark_mouse.gmt`

Record the exact source release and SHA-256 in the generated enrichment manifest. Do not mix human and mouse collections.

The current GSE67462 analysis is mouse and the validated expression space uses mouse gene symbols. The upstream experiment is the mouse OSKM secondary-reprogramming time series with two replicates per time point and integrated expression/Oct4/histone profiling. See NCBI GEO GSE67462 and the associated BioProject PRJNA279997.

For the first run, `Hallmark_mouse.gmt` may be supplied as `--nuisance`; this is a diagnostic context screen, not a declaration that every Hallmark pathway is a nuisance.
