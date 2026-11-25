# Public data sources

`config/public_sources.json` contains the exact URLs, byte counts, and SHA-256 values. Downloads reject HTML responses and checksum mismatches. Public matrices and author-processed data are the starting points; this workflow does not realign sequencing reads where the paper used deposited matrices.

| Cohort/reference | Public source | Role |
|---|---|---|
| GSE235046 | NCBI GEO count table and sample annotations | Source scores, DE and GSEA |
| TCGA-KIRC | UCSC Xena HiSeqV2, GDC primary clinical records, PanCanAtlas ABSOLUTE | Clinical and comparator models |
| TCGA immune composition | GDC PanImmune publication files | CIBERSORT and leukocyte fractions |
| TCGA STAR TPM | GDC Xena release metadata and GENCODE v36 probe map | EPIC linear-TPM input |
| CPTAC ccRCC | LinkedOmics RNA/protein matrices; Clark 2019 publisher Tables S1/S7 | Expression, clinical data, published immune subtypes |
| CPTAC phosphoproteome/TPM | Public Zenodo record 8394329, UMich and Broad exports | Phosphosites and EPIC |
| E-MTAB-1980 | BioStudies matrix and Sato clinical supplement | Independent array cohort |
| TRACERx Renal | Version-pinned `sanroman-24/tx100_rna_2024` deposit | Primary regions, clinical and evolution annotations |
| GSE159115 | GEO 10x matrices and author ccRCC annotations | Patient pseudobulks, UMAP and communication |
| GSE210041 | GEO GSM6415705 and GSM6415706 | Spatial counts and coordinates |
| JAVELIN Renal 101 | Motzer 2020 Nature Medicine publisher supplementary workbook | Clinical and molecular trial data |
| CheckMate 025 | Braun public supplementary workbook | Clinical and molecular trial data |
| ClearCode34 | Commit-pinned `aedin/KidneyCancerSubtypes` | Training matrix and scaling factors; publisher Table S3 supplies benchmark labels |
| EPIC | Commit-pinned `GfellerLab/EPIC` source archive | TRef reference and deconvolution implementation |
| UniProt/OmniPath | Checksum-locked sequence and curated phosphorylation downloads | Canonical peptide/site mapping and kinase scores |
| Pathway sets | MSigDB v2024.1 human GMTs and explicit custom/mouse definitions | Comparator and enrichment tests |

The workflow downloads these files automatically when `--download` is specified. It does not require the user's historical ZIP, an HPC scratch directory, or a CPTAC API token. Original notebook path aliases are created by matching public-source hashes in `config/original_public_source_manifest.csv`.

GDC clinical data are a dynamic API source: the request and full response are saved under `data/raw/tcga_kirc_public/`, and the observed checksum is incorporated into the runtime provenance manifest. Preserve that response to repeat a particular run exactly. Updated GDC records or intentional source updates require rerunning downstream results.

Downloaded data and third-party reference resources retain their originating terms. Cite the primary studies and software as well as the manuscript.
