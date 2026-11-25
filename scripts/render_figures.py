#!/usr/bin/env python3
"""Render figures from existing computed tables without rerunning analyses."""
import argparse
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from mitoxy import figures as main, supplementary_plots as supp, reproduction_plots as cell
JOBS = {
    '1':main.source, '2':main.clinical, '3':main.specificity,
    '4':main.architecture, '5':main.evolution, '6':supp.immune_context,
    '7':supp.treatment, 'S1':supp.clinical_survival, 'S2':supp.comparators,
    'S3':supp.architecture, 'S4':supp.evolution, 'S5':cell.cellular,
    'S6':supp.spatial, 'S7':cell.umap, 'S8':supp.communication,
}
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--figures', nargs='+', choices=list(JOBS), default=list(JOBS))
a = p.parse_args()
for number in a.figures:
    print('Rendering Figure ' + number, flush=True)
    JOBS[number]()
