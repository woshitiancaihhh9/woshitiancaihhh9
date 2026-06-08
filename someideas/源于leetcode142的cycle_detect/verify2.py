import math, numpy as np
D = np.loadtxt('grid.csv', delimiter=',', skiprows=1, dtype=np.int64)
mu, lam, An, Bn, Fn, o, Adet, Bdet = (D[:,i] for i in range(8))
def ceillog2(x):
    x=np.asarray(x,dtype=np.int64); return np.ceil(np.log2(np.maximum(x,1))).astype(np.int64)
r0=ceillog2(mu+1); rlam=ceillog2(-(-lam//3)); rmin=np.maximum(r0,rlam)
pow_rmin=(1<<rmin).astype(np.int64); k=mu+lam
A_pred=np.where(lam<=pow_rmin,2*k+1,2*k+pow_rmin-1)
bad=np.where(A_pred!=An)[0]
print('mismatch cells (mu,lam,An,A_pred):', [(int(mu[i]),int(lam[i]),int(An[i]),int(A_pred[i])) for i in bad])
half=(1<<np.maximum(rmin-1,0)).astype(np.int64)
ok=((o==0)|(o==half-1))
print('o in {0, 2^{rmin-1}-1} violations:', int(np.sum(~ok)))
sel=np.where(lam>pow_rmin)[0][:6]
print('sample lam-dom (mu,lam,o,2^{rmin-1}-1):',[(int(mu[i]),int(lam[i]),int(o[i]),int(half[i]-1)) for i in sel])

SQRT2PI=math.sqrt(2*math.pi)
def Phi(t): return 0.5*(1+math.erf(t/math.sqrt(2)))
def npdf(t):return math.exp(-0.5*t*t)/SQRT2PI
def Aanti(x,e):t=x+e; return t*Phi(t)+npdf(t)
def H(a,b,c,d):
    if b<=a or d<=c:return 0.0
    return SQRT2PI*((Aanti(b,d)-Aanti(a,d))-(Aanti(b,c)-Aanti(a,c)))
def F_of_phi(phi):
    tot=0.0
    for l in range(-60,9):
        Q=2.0**(l-phi); cy=3.0*2.0**(l-1-phi); dy=3.0*2.0**(l-phi)
        tot+=Q*H(0.0,2.0**(l-phi),cy,dy)
        M2=2.0**(l+1-phi); tot+=M2*H(2.0**(l-phi),M2,max(cy,M2),dy)
    return tot
Ek=math.sqrt(math.pi/2)
print('\nexact c_A validation table:')
for n in [10**6,10**8,10**10,10**12,10**14]:
    phi=(0.5*math.log2(n))%1.0; print(f'  n=1e{int(round(math.log10(n))):2d}  phi={phi:.3f}  c_A={2.0+F_of_phi(phi)/Ek:.6f}')

# chunked Monte-Carlo
rng=np.random.default_rng(7)
tot=dict(A=0.0,B=0.0,F=0.0,k=0.0,leB=0,leF=0,n=0)
CH=2_000_000; NCH=10
def cl2(z): return np.ceil(np.log2(np.maximum(z,1.0)))
for _ in range(NCH):
    s=np.sqrt(-2*np.log(rng.random(CH))); V=rng.random(CH); x=s*V; y=s*(1-V)
    m=rng.integers(20,40,CH).astype(np.float64); sqn=2.0**m
    mu_=x*sqn; lam_=np.maximum(y*sqn,1.0); kk=mu_+lam_
    rmin_=np.maximum(cl2(mu_+1),cl2(np.ceil(lam_/3))); prm=2.0**rmin_
    A_=2*kk+np.where(lam_<=prm,1.0,prm-1.0)
    B_=2.0**cl2(np.maximum(lam_,mu_+1))-1+2*lam_+2*mu_
    F_=3*(lam_*np.ceil(mu_/lam_))+2*mu_
    tot['A']+=A_.sum(); tot['B']+=B_.sum(); tot['F']+=F_.sum(); tot['k']+=kk.sum()
    tot['leB']+=int(np.sum(A_<=B_)); tot['leF']+=int(np.sum(A_<=F_)); tot['n']+=CH
print('\nMonte-Carlo average constants (N=%d):'%tot['n'])
print('  c_A=',round(tot['A']/tot['k'],5),' c_B=',round(tot['B']/tot['k'],5),' c_F=',round(tot['F']/tot['k'],5))
print('  A<=B frac=',round(tot['leB']/tot['n'],5),' A<=F frac=',round(tot['leF']/tot['n'],5))

NP=1024; phis=np.arange(NP)/NP
cAphi=np.array([2.0+F_of_phi(p)/Ek for p in phis])
np.save('cAphi.npy',np.vstack([phis,cAphi]))
print('\nsaved cAphi.npy mean=%.6f pp=%.3e'%(cAphi.mean(),cAphi.max()-cAphi.min()))
