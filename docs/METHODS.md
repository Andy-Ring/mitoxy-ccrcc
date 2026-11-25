# Computational methods

## General conventions

The patient is the inferential unit for human cohort analyses. Raw UMI counts, published expression matrices, study clinical workbooks, and author annotations are acquired from the public sources in `config/public_sources.json`. Every fixed file has a SHA-256 checksum. The GDC clinical API request and response are preserved; only explicitly primary diagnoses are used. Missing values are not filled with manuscript estimates. Models use complete cases and report sample sizes.

The four genes are BAK1, BAX, BID, and SLC7A11. Gene standardization generally uses population SD (`ddof=0`). Equal4 is an equal-weight mean of gene z scores, with a further patient-score standardization where a model reports effects per score SD. Core3 excludes SLC7A11. The main author-label cellular analysis uses sample SD (`ddof=1`) across eligible pseudobulks; the population-SD sensitivity is exported separately. Do not equate a mean gene-z unit with a subsequently standardized score SD.

Continuous/ordinal outcomes use OLS with HC3 covariance and normal-reference inference unless explicitly labelled otherwise. TCGA immune-proxy models use HC3 covariance with residual-degrees-of-freedom t inference. Binary outcomes use logistic regression; survival uses Cox models with Breslow ties. Confidence limits for Spearman correlations use a descriptive Fisher-z approximation. BH corrections are applied within the families indicated by each output table and implementation.

## Source experiment

GSE235046 has 15 libraries, three per condition: media, LPS, carbon starvation, LPS plus starvation, and the combined condition plus Torin. Low-count-filtered counts undergo median-ratio normalization and log transformation for scores/PCA. The LPS-by-starvation interaction uses the four factorial conditions. There is no Torin-only condition. Matched-null analyses draw 10,000 modules using the submitted matching algorithm and seed.

PyDESeq2 fits the reconstructed differential-expression models. Six contrasts and the fixed mouse gene sets in `config/` are used for GSEApy preranking with 1,000 permutations. Full gene statistics are exported. `Torin_vs_Combined` means Torin minus combined; negative NES denotes a program reduced with Torin. Figure 1F reverses that sign solely to display attenuation on a positive axis. Figure 1E retains the original contrast direction. FDR values equal to zero at finite permutation resolution should be reported as below resolution, not as exact zero.

## Bulk cohorts and clinical models

TCGA uses the public Xena HiSeqV2 log-expression matrix, GDC primary-cancer clinical data, and published ABSOLUTE purity. Global PC1 is calculated from standardized expression of the 5,000 most variable genes. The primary comparator models adjust for age, sex, purity, and PC1. The ClearCode34 implementation uses the authors' public training matrix, median factors, sample-wise R-style scaling, and threshold-zero nearest-centroid prediction; predicted labels are exported. A strict parser extracts all 380 reference calls and posterior probabilities from the checksum-locked Brooks publisher Table S3 DOC file; agreement and score-by-class summaries are recomputed.

CPTAC uses LinkedOmics log2-FPKM and proteome matrices and the Clark 2019 publisher Table S1. The histology-defined clear-cell cohort has 103 cases. Published age bins are represented by their lower bound; this approximation is explicit in model column names. Primary clinical models adjust for age and sex.

E-MTAB-1980 uses the published quantile-normalized array matrix, fixed probes, and Sato clinical supplement. The fixed mappings are BAK1 A_23_P145357; BAX A_33_P3341144; BID mean of A_24_P187948 and A_33_P3227041; SLC7A11 A_33_P3242623. Models adjust for age and sex; T category is distinguished from overall stage.

TRACERx uses the deposited `tumour_tpm.rds` and author primary-region annotation. The canonical reconstruction standardizes the deposited regional expression values, averages regional gene-z scores, then averages regions within patients and standardizes patient scores. Alternative architectures and log-transform sensitivities are labelled separately. Grade/stage models adjust for age and purity. Evolution models adjust for age, purity, grade, stage, WGI and genomic ITH. In the BAK1-versus-Rest3 joint model, complete-case gene components are standardized but their Rest3 mean is not restandardized; its coefficient is per mean-component unit. Regional dispersion analyses adjust for region number and mean purity.

Eight death/stress programs are tested separately, jointly, through comparator PC1, and after removing the four score genes. The same complete-case patients are used for each comparison family. Standalone program coefficients and VIFs are also exported. Gene definitions and mapped coverage are recorded.

TCGA Hallmark score tests use 10,000 Freedman–Lane residual-row permutations, common patient permutations across pathways, and +1 correction. Competitive tests use HC3 gene statistics, exclusion of the four score genes, 1,000 grade-stratified label permutations, same-sign null normalization, and BH correction. They answer different statistical questions. The 45-gene program is evaluated as both a mean-gene-z score and rank-normalized ssGSEA, with and without four-gene overlap.

Reconstructed earlier-revision sensitivities in Table S3 retain age/sex/purity adjustment for TCGA, without global PC1. Primary comparator models in Table S1 include PC1. These specifications must not be interchanged. Common-patient component models and immune-taxonomy adjustments are recomputed from the same public matrices and published subtype assignments.

## Immune composition and cellular analyses

TCGA CIBERSORT relative fractions and leukocyte fractions come from the GDC PanImmune publication files. Aliquots are restricted to primary tumors and averaged per patient. Tissue proxies multiply relative fractions by leukocyte fraction. Missing constituent fractions remain missing, including for grouped proxies.

EPIC TRef is applied to linear TPM. TCGA STAR TPM is restored from the explicitly documented Xena log2(TPM+1) values and mapped through the public GENCODE v36 probe map. CPTAC Broad transcript TPM is summed by gene with GENCODE v34 and matched through public sample descriptions. Only EPIC convergence-code-zero samples enter models. Fractions, fit diagnostics, cohort coverage, and R session information are retained. The original clinical covariates are used for these Table S3 sensitivities.

GSE159115 public matrices are matched exactly to author cell identifiers. Raw counts are summed by patient and author label (minimum 20 cells), normalized to CPM, and log2(CPM+1) transformed. Paired macrophage/tumor tests enumerate signed-rank sign permutations. The broad-compartment sensitivity uses the available tumor annotations, excludes author filter/doublet flags, merges related labels, requires 30 cells, uses population SD, and bootstraps paired patient differences. It is a tumor-only reconstruction, distinct from a normal-tissue reference analysis.

The UMAP is calculated from public counts using the parameters exported in `results/GSE159115/UMAP_methods.json`. LIANA rank aggregate, library-size-adjusted macrophage splitting, and the CellChat method are run per patient. Resource versions affect interaction membership; full results and the versioned Python environment are retained. CellChat here denotes the method implementation available through LIANA, not an implicit claim to have run the separate R package.

## Spatial definition sensitivity

GSE210041 counts are normalized to 10,000 per spot, log1p transformed, and gene standardized per slide. Equal4 and neighboring MYD88/TNF/TLR2/TLR4 scores are compared using six nearest myeloid-enriched spots and compartment-preserving permutations. All marker panels, thresholds, seeds and attempted definitions are exported.

Figure S6F compares the recoverable archived epithelial-marker definition with the explicitly specified hypoxia-target definition. It does not reconstruct the undocumented submitted marker definition. Grid intervals describe variation across 1,512 definitions per tumor, not confidence intervals. For random marker families, both attempted and evaluable counts are reported; undefined/undersized compartments are not assigned artificial correlations.

## Protein and phosphosite analyses

CPTAC RNA/protein layer correlations use the 103-case clear-cell reference cohort. The 11 mTOR-context targets are fixed in `proteomics.py`. The public University of Michigan phosphoproteome is converted to reference-channel-relative log abundance, aliquots are mapped to patients, and non-tumor/QC rows are excluded. Peptide modifications must map uniquely to a canonical UniProt sequence; export numbering is not assumed to be canonical. No missing phosphosite abundance is imputed.

`config/phosphosites.json` defines the reconstructed panel. Six recoverable historical sites are retained. The historical MTOR T2474 label is absent from this export; MTOR S2448 is analysed as a distinct, explicitly selected sensitivity, not as a renaming of T2474. Raw and total-protein-adjusted rank associations each have a seven-test BH family. Sample sizes differ by site. Claims about the historical panel cannot be transferred to this reconstructed panel.

Kinase substrate analysis uses checksum-locked OmniPath phosphorylation edges and UniProt sequences, unique modified-peptide mapping, total-protein residualization, at least 80% RNA-cohort coverage per site, at least five sites, and at least three substrate proteins. Sites are averaged within proteins before protein-level averaging. Complete-case activity models and eligibility counts are exported; ineligible kinases are not fitted.

## Clinical trials

JAVELIN uses the publisher's clinical and log2-TPM workbook. PFS events are `1 − PFS_P_CNSR`, according to the released censoring convention. The submitted subgroup and formal interaction specifications are retained; supplementary within-arm models use age, sex, PD-L1 and pooled score scaling/cutoff. Cluster-adjusted models exclude missing cluster labels and retain the observed `Other` category explicitly. The nine-score interaction family is exported in full.

CheckMate uses exact canonicalized identifiers and the released CM-025 RNA/clinical intersection. Its released `OS_CNSR`/`PFS_CNSR` fields are interpreted according to the study-specific event convention implemented and documented in `checkmate.py`, rather than assuming the JAVELIN convention. Continuous and within-arm-median models, pooled-median sensitivities, all 24 score-by-treatment tests, and declared FDR families are exported.
