import csv, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

M=300; L=300
BA=np.full((M+1,L+1),np.nan)
FA=np.full((M+1,L+1),np.nan)
src=np.full((M+1,L+1),np.nan)  # 0 tie, 1 localization gain, 2 detection gain
violB=0; violF=0; tieAB=0; tieAF=0
with open('/data/build/grid.csv') as f:
    r=csv.DictReader(f)
    for row in r:
        mu=int(row['mu']); lam=int(row['lam'])
        An=int(row['An']); Bn=int(row['Bn']); Fn=int(row['Fn'])
        o=int(row['o']); Adet=int(row['Adet']); Bdet=int(row['Bdet'])
        BA[mu,lam]=Bn-An
        FA[mu,lam]=Fn-An
        if Bn-An<0: violB+=1
        if Fn-An<0: violF+=1
        if An==Bn: tieAB+=1
        if An==Fn: tieAF+=1
        if An==Bn:
            src[mu,lam]=0
        elif o>0:
            src[mu,lam]=1
        elif Adet<Bdet:
            src[mu,lam]=2
        else:
            src[mu,lam]=1  # fallback (shouldn't happen)
print('violations B-A<0:',violB,' F-A<0:',violF,' tieAB:',tieAB,' tieAF:',tieAF)

def heat(mat,title,fname,vmaxq=99):
    fig,ax=plt.subplots(figsize=(4.2,3.6))
    data=mat[:, 1:].T  # rows=lam, cols=mu
    vmax=np.nanpercentile(data,vmaxq)
    im=ax.imshow(data,origin='lower',aspect='auto',cmap='viridis',
                 extent=[0,M,1,L],vmin=0,vmax=vmax)
    ax.set_xlabel(r'$\mu$'); ax.set_ylabel(r'$\lambda$'); ax.set_title(title,fontsize=10)
    cb=fig.colorbar(im,ax=ax); cb.ax.tick_params(labelsize=7)
    fig.tight_layout(); fig.savefig(fname); plt.close(fig)
    print('wrote',fname,'min',np.nanmin(data),'max',np.nanmax(data))

heat(BA, r'$B_{\rm next}-A_{\rm next}$ (Brent $-$ proposed)', '/data/build/figs/brent_minus_nextval_next.pdf')
heat(FA, r'$F_{\rm next}-A_{\rm next}$ (Floyd $-$ proposed)', '/data/build/figs/floyd_minus_nextval_next.pdf')

# source-of-advantage categorical map
fig,ax=plt.subplots(figsize=(5.0,4.0))
cmap=ListedColormap(['#d9d9d9','#1b9e77','#d95f02'])
norm=BoundaryNorm([-0.5,0.5,1.5,2.5],cmap.N)
data=src[:,1:].T
im=ax.imshow(data,origin='lower',aspect='auto',cmap=cmap,norm=norm,extent=[0,M,1,L])
ax.set_xlabel(r'$\mu$'); ax.set_ylabel(r'$\lambda$')
ax.set_title('Source of strict advantage over Brent',fontsize=10)
cb=fig.colorbar(im,ax=ax,ticks=[0,1,2])
cb.ax.set_yticklabels(['tie','localization gain ($o>0$)','detection gain'],fontsize=8)
fig.tight_layout(); fig.savefig('/data/build/figs/strict_dominance_source.pdf'); plt.close(fig)
print('wrote source fig')
