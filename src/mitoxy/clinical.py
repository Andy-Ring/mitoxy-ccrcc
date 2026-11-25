"""Clinical cohort construction from deposited study metadata."""
import re
import pandas as pd
from .data import ROOT

def cptac_clinical(expression_ids):
    source=ROOT/'data/raw/CPTAC/Clark2019_TableS1.xlsx'
    d=pd.read_excel(source,sheet_name='ccrcc_clinical_characteristics',dtype=str)
    d.columns=d.columns.str.strip()
    def column(aliases):
        normalize=lambda s:re.sub('[^a-z0-9]','',str(s).lower())
        names={normalize(c):c for c in d}
        for alias in aliases:
            if normalize(alias) in names:return names[normalize(alias)]
        raise ValueError(f'Clinical column missing or ambiguous: {aliases}; columns={list(d)}')
    h=column(['histologic_type','histologic type'])
    d=d.loc[d[h].str.strip().str.lower().eq('clear cell renal cell carcinoma')].copy()
    ids=set(map(str,expression_ids));overlap={c:len(set(d[c].dropna().str.strip())&ids) for c in d}
    best=max(overlap.values());candidates=[c for c,n in overlap.items() if n==best]
    if best<100 or len(candidates)!=1:raise ValueError('CPTAC case identifier cannot be uniquely matched')
    d.index=d[candidates[0]].str.strip();d.index.name='Case_ID'
    if not d.index.is_unique or len(d)!=103:raise ValueError('Unexpected CPTAC clear-cell cohort')
    out=pd.DataFrame(index=d.index)
    out['Age']=d[column(['Age','Age at diagnosis'])]
    out['Gender']=d[column(['Gender','Sex'])]
    out['Grade']=d[column(['Grade','Tumor grade','Fuhrman grade'])]
    out['Tumor_Stage_Pathological']=d[column(['Tumor_Stage_Pathological','Pathologic stage','Stage'])]
    def grade(x):
        m=re.search(r'([1-4])',str(x));return float(m.group(1)) if m else float('nan')
    def stage(x):
        s=str(x).upper().replace('STAGE','').strip()
        for roman,n in [('IV',4),('III',3),('II',2),('I',1)]:
            if s.startswith(roman):return float(n)
        return grade(s)
    out['grade_num']=out.Grade.map(grade);out['stage_num']=out.Tumor_Stage_Pathological.map(stage)
    out.to_csv(ROOT/'data/derived/CPTAC_clinical_from_public.csv')
    return out

def primary_diagnoses(case):
    """Use the index cancer, excluding prior/synchronous cancers and recurrence."""
    records=case.get('diagnoses') or []
    primary=[d for d in records if str(d.get('classification_of_tumor','')).strip().lower()=='primary']
    if not primary:
        raise ValueError('No explicitly primary GDC diagnosis for '+str(case.get('submitter_id')))
    return primary
