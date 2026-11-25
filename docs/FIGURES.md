# Figure and table map

All plots read calculated numerical outputs. Historical numerical-directory names are retained only for the submitted notebooks; the output filenames use current manuscript numbering.

| Output | Stage | Generated file |
|---|---|---|
| Figure 1 | `figure1` | `figures/Figure_1.pdf` |
| Figure 2 | `figure2` | `figures/Figure_2.pdf` |
| Figure 3 | `figure3` | `figures/Figure_3.pdf` |
| Figure 4 | `figure4` | `figures/Figure_4.pdf` |
| Figure 5 | `figure5` | `figures/Figure_5.pdf` |
| Figure 6 | `figure6` | `figures/Figure_6.pdf` |
| Figure 7 | `figure7` | `figures/Figure_7.pdf` |
| Figure S1 | `supplement1` | `figures/Figure_S1.pdf` |
| Figure S2 | `supplement2` | `figures/Figure_S2.pdf` |
| Figure S3 | `supplement3` | `figures/Figure_S3.pdf` |
| Figure S4 | `supplement4` | `figures/Figure_S4.pdf` |
| Figure S5 | `cellular_plots` | `figures/Figure_S5.pdf` |
| Figure S6 | `supplement6` | `figures/Figure_S6.pdf` |
| Figure S7 | `umap_plots` | `figures/Figure_S7.pdf` |
| Figure S8 | `supplement8` | `figures/Figure_S8.pdf` |
| Table S1 | `tables` | `tables/Supplementary_Table_S1.xlsx` |
| Table S2 | `tables` | `tables/Supplementary_Table_S2.xlsx` |
| Table S3 | `tables` | `tables/Supplementary_Table_S3.xlsx` |

Run a selected output with `python scripts/reproduce.py --download --stages STAGE`. Tables S1–S3 are divided into numbered sections in each workbook’s Table index. `config/table_sources.json` provides the full sheet-to-result mapping. Full source DE, spatial-definition and interaction results are also exported as CSV files.

Figure S6F uses the archived-versus-hypoxia comparison, and Figure 1F displays the sign-reversed Torin contrast as attenuation. Figure S1B separates linear coefficients and odds ratios onto different axes. Consult the Methods before adapting captions.

## Reference styling restored (23 November 2026)

Each figure has an executable plotting notebook in `notebooks/figures/`.
The notebooks call shared plotting functions, with consistent fonts, palette,
light grids, alternating forest-plot bands, and result boxes defined in
`src/mitoxy/plot_style.py`. All labels and estimates use the current numerical
outputs. Figure 2A contains four cohort cards with data-derived counts and no
score-definition block, gene badges, or connecting lines.

After computing the inputs, render every figure without rerunning analysis:

```bash
python scripts/render_figures.py
# Or render selected figures:
python scripts/render_figures.py --figures 2 6 S6
```

For the delivered cached-output bundle, copy its `data/` and `results/` folders
into the repository before rendering. S6 also requires the checksum-verified
GSE210041 spatial archives under `data/raw/GSE210041/`; these are included in the
updated output bundle and retrieved by the normal public-data workflow.

Restored visual structures include the Figure 5 full-width overview/regional
panels and stacked variance bar, Figure 6 patient component dots and cluster
violins, S3 architecture trajectories, S5 horizontal significance bars, S6
histology overlays and summary table, and the seven-panel S8 layout. Figure 1E
again shows the reference's nine pathway rows using current GSEA values.
S8E shows the reference ligand panel using current main/UMI-adjusted results.
The current top-ten-per-direction interaction selection in S8D is retained.

Analytical corrections are retained: S1B separates coefficient and odds-ratio
axes; Figure 4D retains the correct predictor-unit label; S6F retains the
reconstructable archived-versus-hypoxia definitions; Figure 1F retains the
sign-reversed attenuation contrast. Figure 3 retains all current comparator
rows, Figure 5E/F retain the available current metrics/models, and Figure 7D
retains the ORR, PFS, and OS interaction estimates. Different current results
can therefore change the plotted positions, intervals, counts, and row content
relative to the historical reference. S4D displays the existing log-odds
estimates and intervals exponentiated as odds ratios, without refitting.

The manuscript and calculated analysis tables are unchanged by this styling
revision. PNG and PDF are exported at 300 dpi (vector elements remain vector
in PDF). Obsolete TIFF output is not distributed.

Validation: all 15 plotting notebooks passed notebook-schema validation and direct code-cell execution. All 15 figures rendered successfully. The 218 pre-existing processed/result input files were verified unchanged byte-for-byte. Main and supplementary layouts were visually inspected, and PDF text confirmed removal of the Figure 2A score-definition block. Interactive Jupyter kernel startup was unavailable in the validation environment; code-cell execution did not require a kernel server.
