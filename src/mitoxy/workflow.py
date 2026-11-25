"""Dependency-aware manuscript computation with provenance and explicit failures."""
import hashlib,json,os,platform,subprocess,sys,time
from datetime import datetime,timezone
from .data import ROOT,sha256,download_public
from .notebooks import execute_notebook

# Dependencies describe scientific inputs, not cached figure files.
STAGES={
 'source_scores':dict(groups=['GSE235046'],notebook='source'),
 'source_de':dict(groups=['GSE235046']),
 'source_pathways':dict(depends=['source_de']),
 'tcga':dict(groups=['TCGA-KIRC'],notebook='tcga'),
 'tracerx':dict(groups=['TRACERx']),
 'cptac':dict(groups=['CPTAC','CPTAC clinical','genesets']),
 'clinical':dict(groups=['TCGA-KIRC','CPTAC','CPTAC ccRCC','CPTAC clinical','E-MTAB-1980','TRACERx'],notebook='clinical'),
 'survival':dict(depends=['clinical'],notebook='survival'),
 'regional':dict(depends=['tracerx'],notebook='evolution'),
 'evolution_sensitivity':dict(depends=['regional'],notebook='evolution_sensitivity'),
 'tcga_comparators':dict(depends=['tcga'],groups=['genesets']),
 'cptac_comparators':dict(depends=['cptac'],groups=['genesets']),
 'tracerx_comparators':dict(depends=['tracerx'],groups=['genesets']),
 'evolution':dict(depends=['regional','tracerx_comparators']),
 'cellular':dict(groups=['GSE159115']),
 'umap':dict(depends=['cellular']),
 'spatial':dict(groups=['GSE210041']),
 'checkmate':dict(groups=['CheckMate025']),
 'communication':dict(depends=['cellular']),
 'trials':dict(groups=['JAVELIN'],notebook='trials'),
 'immune':dict(depends=['tcga_comparators'],groups=['TCGA-immune']),
 'signatures':dict(depends=['tcga_comparators','cptac_comparators','tracerx_comparators']),
 'membership':dict(depends=['signatures']),
 'immune_subtypes':dict(depends=['tcga_comparators','cptac_comparators'],groups=['CPTAC supplementary','CPTAC-package']),
 'component_models':dict(depends=['tcga_comparators','cptac_comparators','tracerx_comparators']),
 'compartments':dict(depends=['cellular']),
 'kinase_activity':dict(depends=['proteomics'],groups=['kinases']),
 'tcga_sensitivities':dict(depends=['tcga_comparators'],groups=['ClearCode34']),
 'javelin_context':dict(depends=['trials']),
 'architecture':dict(depends=['clinical','tracerx','javelin_context','tcga_sensitivities']),
 'proteomics':dict(depends=['tcga_comparators','cptac_comparators'],groups=['CPTAC-package','UniProt','CPTAC ccRCC']),
 'epic':dict(depends=['tcga_comparators','cptac_comparators'],groups=['TCGA-TPM','CPTAC-package','EPIC']),
 'pathway_scores':dict(depends=['tcga_comparators']),
 'competitive_gsea':dict(depends=['tcga_comparators']),
 'figure6':dict(depends=['cellular','immune','javelin_context']),
 'figure7':dict(depends=['trials','checkmate']),
 'supplement1':dict(depends=['survival']),
 'supplement3':dict(depends=['architecture']),
 'supplement6':dict(depends=['spatial']),
 'supplement2':dict(depends=['tcga_comparators','cptac_comparators']),
 'supplement4':dict(depends=['evolution_sensitivity','evolution']),
 'supplement8':dict(depends=['communication']),
 'figure2':dict(depends=['clinical']),
 'figure1':dict(depends=['source_scores','source_pathways']),
 'figure3':dict(depends=['tcga_comparators','cptac_comparators']),
 'figure4':dict(depends=['tracerx']),
 'figure5':dict(depends=['evolution']),
 'cellular_plots':dict(depends=['cellular']),
 'umap_plots':dict(depends=['umap']),
}

STAGES['tables']=dict(depends=list(STAGES))

def order(requested):
    result=[];visiting=set()
    def visit(stage):
        if stage not in STAGES:raise ValueError('Unknown stage: '+stage)
        if stage in result:return
        if stage in visiting:raise ValueError('Cyclic dependency: '+stage)
        visiting.add(stage)
        for parent in STAGES[stage].get('depends',[]):visit(parent)
        visiting.remove(stage);result.append(stage)
    for stage in requested:visit(stage)
    return result

def run(stage):
    spec=STAGES[stage]
    if spec.get('notebook'):
        from .recovered_inputs import prepare
        prepare(include_gdc=stage in ['tcga','clinical','molecular','survival'])
        return execute_notebook('cohorts/'+spec['notebook']+'/analysis.ipynb')
    if stage=='immune':
        from .immune import analyze
        return analyze()
    if stage in ['signatures','membership','tcga_sensitivities','javelin_context','architecture','immune_subtypes','component_models']:
        from . import extensions
        return getattr(extensions, 'javelin' if stage=='javelin_context' else stage)()
    if stage=='compartments':
        from .cellular import compartments
        return compartments()
    if stage=='kinase_activity':
        from .proteomics import kinase_activity
        return kinase_activity()
    if stage=='proteomics':
        from .proteomics import analyze
        return analyze()
    if stage=='epic':
        from .deconvolution import analyze
        return analyze()
    if stage in ['pathway_scores','competitive_gsea']:
        from . import pathways
        return (pathways.patient_scores if stage=='pathway_scores' else pathways.competitive)()
    if stage in ['figure6','supplement1','supplement3','supplement6']:
        from . import supplementary_plots as p
        return {'figure6':p.immune_context,'supplement1':p.clinical_survival,'supplement3':p.architecture,'supplement6':p.spatial}[stage]()
    if stage=='tracerx':
        from .analysis import tracerx
        return tracerx()
    if stage=='cptac':
        from .analysis import cptac
        return cptac()  # strict public clinical source; no silent snapshot fallback
    if stage=='source_de':
        from .source import differential_expression
        return differential_expression()
    if stage=='source_pathways':
        from .source import pathways
        return pathways()
    if stage.endswith('_comparators'):
        from .bulk import analyze
        return analyze({'tcga':'TCGA','cptac':'CPTAC','tracerx':'TRACERx'}[stage.split('_')[0]])
    if stage=='evolution':
        from .evolution import analyze
        return analyze()
    if stage=='cellular':
        from .cellular import build_counts,pseudobulk
        return pseudobulk()
    if stage=='umap':
        from .cellular import embedding
        return embedding()
    if stage=='spatial':
        from .spatial import analyze
        return analyze()
    if stage=='checkmate':
        from .checkmate import analyze
        return analyze()
    if stage=='communication':
        scripts=ROOT/'scripts/cell_communication';base=['--input-dir',str(ROOT/'data/derived/GSE159115')]
        for name,extra in [('main',[]),('umi_adjusted',['--umi-adjusted'])]:
            subprocess.run([sys.executable,str(scripts/'liana_gse159115.py'),*base,'--output-dir',str(ROOT/'results/LIANA'/name),*extra],check=True,cwd=ROOT)
        return subprocess.run([sys.executable,str(scripts/'cellchat_confirm.py'),*base,'--output-dir',str(ROOT/'results/LIANA/cellchat')],check=True,cwd=ROOT)
    if stage in ['figure7','supplement2','supplement4','supplement8']:
        from . import supplementary_plots
        return {'figure7':supplementary_plots.treatment,'supplement2':supplementary_plots.comparators,'supplement4':supplementary_plots.evolution,'supplement8':supplementary_plots.communication}[stage]()
    if stage=='tables':
        return subprocess.run([sys.executable,str(ROOT/'scripts/export_tables.py')],check=True,cwd=ROOT)
    if stage=='figure2':return execute_notebook('figures/figure2.ipynb')
    if stage.startswith('figure'):
        from . import figures
        return {'figure1':figures.source,'figure3':figures.specificity,'figure4':figures.architecture,'figure5':figures.evolution}[stage]()
    if stage.endswith('_plots'):
        from . import reproduction_plots
        return getattr(reproduction_plots,stage.removesuffix('_plots'))()
    raise ValueError('Unimplemented stage: '+stage)

def execute(stage,download=False):
    for name in ['results','data/derived','data/provenance','figures']:(ROOT/name).mkdir(parents=True,exist_ok=True)
    start=time.time();status='failed';error=None
    try:
        groups=STAGES[stage].get('groups',[])
        if download and groups:download_public(groups=groups)
        run(stage);status='computed'
    except Exception as exc:error=f'{type(exc).__name__}: {exc}';raise
    finally:
        output=ROOT/'results/run_records.json'
        records=json.loads(output.read_text()) if output.exists() else []
        code={str(p.relative_to(ROOT)):sha256(p) for folder in ['src','scripts','config','notebooks']
              for p in (ROOT/folder).rglob('*') if p.is_file() and '__pycache__' not in p.parts}
        records.append(dict(stage=stage,status=status,error=error,seconds=time.time()-start,
            utc=datetime.now(timezone.utc).isoformat(),python=platform.python_version(),
            code_sha256=hashlib.sha256(json.dumps(code,sort_keys=True).encode()).hexdigest()))
        output.write_text(json.dumps(records,indent=2)+'\n')

def coverage(require_complete=False):
    components=json.loads((ROOT/'config/manuscript_outputs.json').read_text())
    for row in components:
        row['missing']=[path for path in row['files'] if not (ROOT/path).is_file() or (ROOT/path).stat().st_size==0]
        row['verification']='computed' if not row['missing'] else 'not_computed'
    unresolved=[r for r in components if r['missing']]
    report={'complete':not unresolved,'meaning':'Declared reconstructed outputs exist; this is not equality to historical estimates.','components':components,'unresolved':unresolved}
    if require_complete and unresolved:raise RuntimeError('Missing reconstructed outputs: '+', '.join(r['output'] for r in unresolved))
    return report
