"""Shared reference-figure styling; numerical values always come from inputs."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.ticker import MaxNLocator, ScalarFormatter

BLUE = '#0072B2'
ORANGE = '#D55E00'
TEAL = '#009E73'
PURPLE = '#7A5195'
INK = '#27313D'
GRID = '#DCE3E8'
BAND = '#F0F5F7'
C = [BLUE, ORANGE, TEAL, PURPLE, INK]

def style():
    plt.rcParams.update({
        'font.family':'DejaVu Sans', 'font.size':9, 'axes.labelsize':9,
        'axes.titlesize':9, 'xtick.labelsize':8, 'ytick.labelsize':8,
        'axes.spines.top':False, 'axes.spines.right':False,
        'axes.labelcolor':INK, 'text.color':INK, 'xtick.color':INK,
        'ytick.color':INK, 'axes.edgecolor':'#B7C1CB', 'axes.linewidth':.65,
        'axes.grid':True, 'grid.color':GRID, 'grid.alpha':.45, 'grid.linewidth':.6,
        'axes.axisbelow':True, 'legend.frameon':False, 'legend.fontsize':8,
        'figure.facecolor':'white', 'axes.facecolor':'white',
        'savefig.facecolor':'white', 'pdf.fonttype':42, 'svg.fonttype':'none',
        'figure.constrained_layout.w_pad':.16, 'figure.constrained_layout.h_pad':.16,
        'figure.constrained_layout.wspace':.08, 'figure.constrained_layout.hspace':.12,
    })

def ptext(value):
    if not np.isfinite(value): return 'NA'
    if value == 0: return '< resolution'
    return f'{value:.1e}' if value < .001 else f'{value:.3f}'

def box(ax, x, y, text, color=BLUE, transform=None, **kwargs):
    return ax.text(x,y,text,transform=transform or ax.transData,
        va='center',ha='left',fontsize=7,color=color,linespacing=1.3,
        bbox=dict(boxstyle='round,pad=.3,rounding_size=.12',
                  fc='white',ec=to_rgba(color,.55),lw=.6),**kwargs)

def forest(ax,d,labels,color=BLUE,xlabel='Coefficient (95% CI)',value='beta',
           null=0,symbol='β',log=False,annotate=True,colors=None):
    d=d.reset_index(drop=True)
    palette=colors if colors is not None else [color]*len(d)
    lower=np.nanmin(d.ci_low); upper=np.nanmax(d.ci_high)
    lo=min(lower,null);hi=max(upper,null)
    if log:
        if lower<=0: raise ValueError('Ratio interval must be positive')
        lo=np.log(min(lower,null));hi=np.log(max(upper,null))
    span=max(hi-lo,.1)
    left=lo-.06*span;right=hi+(.65 if annotate else .08)*span
    if log:
        ax.set_xscale('log');ax.set_xlim(np.exp(left),np.exp(right))
        ax.xaxis.set_major_formatter(ScalarFormatter())
    else:ax.set_xlim(left,right);ax.xaxis.set_major_locator(MaxNLocator(5))
    for i,(_,r) in enumerate(d.iterrows()):
        c=palette[i]
        if i%2==0: ax.axhspan(i-.44,i+.44,color=BAND,zorder=0)
        ax.errorbar(r[value],i,xerr=[[r[value]-r.ci_low],[r.ci_high-r[value]]],
                    fmt='o',capsize=3,markersize=4,mec='white',mew=.5,
                    elinewidth=1.2,color=c,zorder=3)
        if annotate:
            lines=[f'{symbol} = {r[value]:.2f}']
            p=next((r[k] for k in ['p_value','p','p_exact'] if k in r),None)
            if p is not None:lines.append('P '+('' if p==0 else '= ')+ptext(p))
            q=next((r[k] for k in ['fdr','q_BH_24','q'] if k in r),None)
            if q is not None: lines.append('q '+('' if q==0 else '= ')+ptext(q))
            box(ax,.76,i,'\n'.join(lines),c,ax.get_yaxis_transform())
    ax.set(yticks=range(len(d)),yticklabels=list(labels),xlabel=xlabel,
           ylim=(len(d)-.45,-.55))
    ax.axvline(null,color='#AAB5C1',ls=(0,(3,3)),lw=.8)
    ax.grid(axis='y',visible=False);ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y',length=0,pad=8)

def letter(ax,text):
    ax.text(-.12,1.07,text,transform=ax.transAxes,weight='bold',fontsize=16,
            ha='left',va='bottom',clip_on=False)
