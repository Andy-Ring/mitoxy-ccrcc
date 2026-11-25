# Reproducibility

Run `python scripts/reproduce.py --download --all` from an installed checkout. The dependency graph builds analysis-ready data, model tables, figures, and labelled supplementary workbooks. Failed dependencies block downstream branches and cause a nonzero exit. `--check-coverage` verifies that all declared reconstructed outputs exist; it does not assert equality to an earlier manuscript.

The source manifest locks fixed public files. Each stage records its status, runtime, Python version, and code hash in `results/run_records.json`. Original notebook executions are retained under `results/notebooks/`. Supplementary table indexes link each sheet to the numerical source file and its checksum.

Reproducible here means recomputation from the documented deposited matrices and annotations. It does not imply raw-read realignment, reconstruction of unavailable experimental metadata, or identical rendering across software versions. Seeds and package versions are fixed where supported. Floating-point/optimizer differences and finite-permutation resolution can affect final digits.

The reconstructed spatial comparison, current phosphosite panel, exact trial missing-category handling, and tumor-only broad compartments are specified in [Methods](METHODS.md). Their results should be interpreted with those definitions. In particular, the old undocumented submitted spatial definition and the absent MTOR T2474 measurement are not silently recreated from retained numerical estimates.
