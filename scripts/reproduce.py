#!/usr/bin/env python3
"""Run manuscript analyses from public inputs; failed branches cannot pass silently."""
import argparse,sys,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from mitoxy.workflow import STAGES,order,execute,coverage
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--stages',nargs='+',choices=list(STAGES))
p.add_argument('--all',action='store_true',help='Run every integrated stage; unresolved manuscript coverage gives a nonzero exit')
p.add_argument('--download',action='store_true',help='Fetch and verify each stage\'s available public inputs')
p.add_argument('--list',action='store_true',help='Show stages and dependencies')
p.add_argument('--check-coverage',action='store_true')
a=p.parse_args()
if a.list:
 print(json.dumps(STAGES,indent=2));sys.exit(0)
if a.check_coverage:
 report=coverage();print(json.dumps(report,indent=2));sys.exit(0 if report['complete'] else 2)
if not a.all and not a.stages:p.error('Choose --stages, --all, --list, or --check-coverage')
failures={}
for stage in order(list(STAGES) if a.all else a.stages):
 blocked=[d for d in STAGES[stage].get('depends',[]) if d in failures]
 if blocked:
  failures[stage]='Blocked by '+', '.join(blocked);print(stage+': '+failures[stage],flush=True);continue
 try:print('Computing '+stage,flush=True);execute(stage,download=a.download)
 except Exception as e:failures[stage]=str(e);print(stage+' FAILED: '+str(e),flush=True)
if failures:
 print(json.dumps({'failed_or_blocked':failures},indent=2));sys.exit(1)
if a.all:coverage(require_complete=True)
print('Selected stages completed. This does not certify unselected manuscript outputs.')
