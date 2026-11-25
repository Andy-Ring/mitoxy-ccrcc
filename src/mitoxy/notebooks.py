"""Execute each notebook in an isolated Python process, including on HPC nodes.

No Jupyter server or listening socket is required. IPython cell semantics and
stdout/stderr are retained; failed cells stop execution and preserve a traceback.
"""
import os,sys,subprocess,traceback
import nbformat
from IPython.core.interactiveshell import InteractiveShell
from IPython.utils.capture import capture_output
from .data import ROOT

def execute_notebook(path):
    env=os.environ.copy();env['PYTHONPATH']=str(ROOT/'src')+os.pathsep+env.get('PYTHONPATH','')
    subprocess.run([sys.executable,'-m','mitoxy.notebooks',path],cwd=ROOT,env=env,check=True)

def worker(relative):
    path=ROOT/'notebooks'/relative;nb=nbformat.read(path,as_version=4)
    out=ROOT/'results/notebooks'/relative;out.parent.mkdir(parents=True,exist_ok=True)
    shell=InteractiveShell();os.chdir(ROOT);execution=0
    try:
        for cell in nb.cells:
            if cell.cell_type!='code':continue
            execution+=1;cell.execution_count=execution
            with capture_output() as captured:
                result=shell.run_cell(cell.source)
            cell.outputs=[nbformat.v4.new_output('stream',name=k,text=v) for k,v in
                          [('stdout',captured.stdout),('stderr',captured.stderr)] if v]
            if result.error_before_exec or result.error_in_exec:
                error=result.error_before_exec or result.error_in_exec
                cell.outputs.append(nbformat.v4.new_output('error',ename=type(error).__name__,
                    evalue=str(error),traceback=traceback.format_exception(error)))
                raise error
    finally:nbformat.write(nb,out)

if __name__=='__main__':worker(sys.argv[1])
