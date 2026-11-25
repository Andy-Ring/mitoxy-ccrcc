"""Structural/statistical checks; no tests force a desired scientific estimate."""
import json
import unittest
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import norm
from mitoxy.analysis import ROOT,GENES,model,zscore
from mitoxy.data import check_source

class AnalysisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not (ROOT/'data/derived/TRACERx_all_patients.csv').exists():raise unittest.SkipTest('Run tracerx stage first')
        cls.p=pd.read_csv(ROOT/'data/derived/TRACERx_all_patients.csv',index_col=0)
        cls.d=pd.read_csv(ROOT/'data/derived/TRACERx_analysis_patients.csv',index_col=0)
        cls.models=pd.read_csv(ROOT/'results/TRACERx_models.csv').set_index('model_id')

    def test_patient_scaling(self):
        np.testing.assert_allclose(self.p[GENES+['Equal4']].std(ddof=0),1,atol=1e-12)
        np.testing.assert_allclose(self.p[GENES+['Equal4']].mean(),0,atol=1e-12)

    def test_regional_score_not_patient_component_reweighted(self):
        r=pd.read_csv(ROOT/'data/derived/TRACERx_regions.csv')
        reconstructed=zscore(r.groupby('Patient')[GENES].mean().mean(axis=1))
        np.testing.assert_allclose(reconstructed.loc[self.p.index],self.p.Equal4,atol=1e-12)

    def test_rest3_units_and_joint_component_scaling(self):
        expected=zscore(self.d[GENES[1:]]).mean(axis=1)
        np.testing.assert_allclose(expected,self.d.Rest3,atol=1e-12)
        self.assertFalse(np.isclose(self.d.Rest3.std(ddof=0),1))
        np.testing.assert_allclose(self.d.BAK1_joint.std(ddof=0),1)

    def test_slc7a11_hc3_formula_and_scale(self):
        d=self.d; x=sm.add_constant(d[['SLC7A11','purity','Age']]).to_numpy();y=d.grade.to_numpy()
        bread=np.linalg.inv(x.T@x);coef=bread@x.T@y;resid=y-x@coef
        h=np.sum((x@bread)*x,axis=1)
        covariance=bread@x.T@np.diag((resid/(1-h))**2)@x@bread
        actual=self.models.loc['TRACERx.grade.SLC7A11']
        self.assertAlmostEqual(actual.se_HC3,np.sqrt(covariance[1,1]),places=10)
        self.assertAlmostEqual(actual.ci_high-actual.beta,norm.ppf(.975)*actual.se_HC3,places=10)
        scaled=d.copy();scaled['SLC7A11']=d.SLC7A11*7
        other=model(scaled,'grade',['SLC7A11','purity','Age'],'SLC7A11','scale_invariance')[0]
        self.assertAlmostEqual(other['beta']*7,actual.beta,places=10)
        self.assertAlmostEqual(other['se_HC3']*7,actual.se_HC3,places=10)

    def test_proliferation_comparisons_use_same_patients(self):
        if not (ROOT/'results/CPTAC_proliferation_models.csv').exists():self.skipTest('Run CPTAC stage first')
        d=pd.read_csv(ROOT/'results/CPTAC_proliferation_models.csv')
        self.assertTrue(d.groupby('outcome').patient_sha256.nunique().eq(1).all())
        self.assertTrue(d.groupby('outcome').n.nunique().eq(1).all())
        self.assertTrue((d.ci_low<=d.beta).all() and (d.beta<=d.ci_high).all())
        self.assertTrue(d.p.between(0,1).all())

    def test_present_public_input_integrity(self):
        checked=0
        for source in json.loads((ROOT/'config/public_sources.json').read_text()):
            path=ROOT/'data/raw'/source['path']
            if path.exists():check_source(path,source);checked+=1
        if not checked:self.skipTest('No downloaded public inputs')

if __name__=='__main__':unittest.main()
