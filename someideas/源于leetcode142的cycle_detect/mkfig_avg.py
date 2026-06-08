import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt, math
phis,cA=np.load('cAphi.npy')
mean=cA.mean()
fig,ax=plt.subplots(figsize=(6.2,3.6))
ax.plot(phis,cA,color='#1f4e79',lw=1.8,label=r'$c_A(\varphi)$ exact')
ax.axhline(mean,color='#c0504d',ls='--',lw=1.0,label=r'$\bar c_A=%.6f$'%mean)
# 1-harmonic fit
C=np.fft.rfft(cA-mean)/len(cA); a1=2*abs(C[1]); p1=math.atan2(C[1].imag,C[1].real)
ax.plot(phis,mean+a1*np.cos(2*math.pi*phis+p1),color='#9bbb59',ls=':',lw=1.4,label='1st harmonic')
ax.set_xlabel(r'$\varphi=\{\frac{1}{2}\log_2 n\}$'); ax.set_ylabel(r'$c_A$')
ax.set_title(r'Average successor-eval constant: log-periodic ripple (amp $%.2f\times10^{-4}$)'%(a1*1e4))
ax.legend(fontsize=8,loc='upper right'); ax.grid(alpha=.3)
fig.tight_layout(); fig.savefig('figs/fig_avg_ripple.pdf'); print('wrote fig_avg_ripple.pdf  mean=%.6f amp=%.4e phase=%.3f pp=%.3e'%(mean,a1,p1,cA.max()-cA.min()))
