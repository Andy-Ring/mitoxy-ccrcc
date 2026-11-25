"""Statistical equivalence and public-input reconstruction checks."""
import unittest
import numpy as np
import pandas as pd
from scipy.sparse import load_npz
from statsmodels.duration.hazard_regression import PHReg
from statsmodels.stats.multitest import multipletests
from mitoxy.data import ROOT
from mitoxy.survival import cox,bh

class SurvivalTests(unittest.TestCase):
    def test_breslow_matches_independent_implementation(self):
        rng=np.random.default_rng(31)
        x=rng.normal(size=(80,2));t=np.ceil(rng.exponential(np.exp(-x[:,0]) * 5)+.1)
        e=rng.binomial(1,.7,size=len(t))
        actual=cox(t,e,x);reference=PHReg(t,x,status=e,ties='breslow').fit()
        np.testing.assert_allclose(actual['beta'],reference.params,atol=1e-7)
        np.testing.assert_allclose(actual['se'],reference.bse,atol=1e-7)

    def test_fractional_events_are_rejected(self):
        with self.assertRaises(ValueError):cox([1,2,3],[1,.5,0],[1,2,3])

    def test_bh_matches_statsmodels(self):
        p=np.array([.02,.001,.5,.02,.7])
        np.testing.assert_allclose(bh(p),multipletests(p,method='fdr_bh')[1])

class PublicCellTests(unittest.TestCase):
    def test_annotation_library_sizes_match_counts(self):
        d=ROOT/'data/derived/GSE159115'
        if not (d/'GSE159115_ccRCC_cell_gene_counts.npz').exists():self.skipTest('Run cellular stage first')
        x=load_npz(d/'GSE159115_ccRCC_cell_gene_counts.npz')
        meta=pd.read_csv(d/'GSE159115_ccRCC_cell_metadata.csv.gz')
        np.testing.assert_allclose(np.asarray(x.sum(axis=1)).ravel(),meta.total_UMI)
        self.assertTrue(meta.cell.is_unique)

    def test_checkmate_q_values_use_declared_families(self):
        p=ROOT/'results/CheckMate025/within_arm_models.csv'
        if not p.exists():self.skipTest('Run CheckMate stage first')
        d=pd.read_csv(p)
        for _,g in d.groupby('family'):
            np.testing.assert_allclose(g.q,multipletests(g.p,method='fdr_bh')[1])
        interaction=pd.read_csv(p.parent/'interaction_models.csv')
        np.testing.assert_allclose(interaction.q_BH_24,multipletests(interaction.p,method='fdr_bh')[1])

if __name__=='__main__':unittest.main()

class ClinicalDiagnosisTests(unittest.TestCase):
    def test_index_cancer_is_independent_of_gdc_list_order(self):
        from mitoxy.clinical import primary_diagnoses
        prior={'classification_of_tumor':'Prior primary','tumor_grade':'G1'}
        kidney={'classification_of_tumor':'Primary','tumor_grade':'G4'}
        for records in [[prior,kidney],[kidney,prior]]:
            self.assertEqual(primary_diagnoses({'submitter_id':'case','diagnoses':records}),[kidney])

    def test_missing_primary_cannot_silently_use_prior_cancer(self):
        from mitoxy.clinical import primary_diagnoses
        with self.assertRaises(ValueError):primary_diagnoses({'diagnoses':[{'classification_of_tumor':'Prior primary'}]})
