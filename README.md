# Mitoxy transcriptional state in clear-cell renal cell carcinoma

Reproducible analysis code for the Mitoxy ccRCC manuscript, based on the [submitted repository](https://github.com/Andy-Ring/mitoxy-ccrcc), supplied revision code, and documented reconstructed methods. The workflow downloads public study data and computes the statistics, seven main figures, eight supplementary figures, and Tables S1–S3.

This version documents the analysis reconstructed from public inputs. Its estimates follow the specifications below and may differ from earlier versions. See [Methods](docs/METHODS.md) for score units and the spatial and phosphosite specifications.

## Install and run

Use Linux, Python 3.12, and R 4.4 (base R is sufficient). A 32 GB compute node is recommended. Public inputs are cached locally; leave at least 20 GB free for inputs and intermediate files.

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m pip install -e .

# Full public-data workflow, including figures and tables.
python scripts/reproduce.py --download --all

# Or run one figure and its scientific dependencies.
python scripts/reproduce.py --download --stages figure4

python -m unittest discover -s tests -v
python scripts/reproduce.py --check-coverage
```

Alternatively, create the combined Python/R environment from the repository root with `micromamba create -f environment.yml` and activate `mitoxy-ccrcc`.

If `Rscript` is not on PATH, set `MITOXY_RSCRIPT` to its full path. The EPIC source release is downloaded with a fixed checksum and loaded by the R script; no R package installation is performed during analysis. Set `MITOXY_THREADS` to the allocated CPU count. A Slurm template is provided in `scripts/run_slurm.sh`.

The runner lists available stages with `--list`, resolves dependencies, and reports failed branches. `--all` also builds the labelled supplementary workbooks and checks declared output coverage. A selected-stage run certifies only that selection. No saved manuscript estimates or artwork are used as analysis inputs.

## Repository and outputs

| Location | Contents |
|---|---|
| `notebooks/` | Ordered entry points and portable analytical notebooks |
| `src/mitoxy/` | Cohort construction, statistics, figures, and workflow |
| `scripts/cell_communication/` | Patient-level LIANA and CellChat-method analyses |
| `config/` | Public URLs/checksums, fixed definitions, figure/table mapping, code provenance |
| `docs/` | Methods, data sources, and validation |
| `tests/` | Statistical and data-integrity checks |
| `figures/` | Generated PDF/PNG Figures 1–7 and S1–S8 |
| `tables/` | Generated Tables S1–S3, section index, and source checksums |
| `results/`, `data/derived/` | Numerical results and analysis-ready data |

Generated outputs and downloaded data are excluded from Git. Full differential-expression tables are under `results/source_de/`; complete spatial definitions and ligand–receptor results remain available as CSV files alongside the supplementary summaries. [Figure mapping](docs/FIGURES.md) connects outputs to commands and source tables.

The input studies are GSE235046, TCGA-KIRC, CPTAC ccRCC, E-MTAB-1980, TRACERx Renal, GSE159115, GSE210041, JAVELIN Renal 101, and CheckMate 025. Consult [Data sources](docs/DATA_SOURCES.md) and cite the originating studies and software. Code provenance is recorded in `config/code_sources.json`; the submitted base commit is `fce9ade29ef6de731407a7fcf14098f30a272cd3`.

## Plotting only

Use `notebooks/figures/` or `python scripts/render_figures.py` to render current
computed results in the restored manuscript style. See `docs/FIGURES.md` for
cached-input setup and the styling-only changes.
