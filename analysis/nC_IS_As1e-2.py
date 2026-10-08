#!/usr/bin/env python3
"""n(C_m) from IS at As=1e-2 (NL=256,npeak=16 == kp of 512/32), local compaction-peak rule (direct code) in window.
Event file gbc/ev_As1e-2_bc<bc>.csv: E,seed,lnW,x,y,z,r,val ; S,seed,lnW,nev.
n(bin)=sum_{events in bin, |x|_inf<=m}/(N |R|)*Wbar/dC, Wbar=exp(mean lnW+var/2). Bootstrap over samples.
usage: nC_IS_As1e-2.py out.csv bc:m [bc:m ...]   (R cube |x|_inf<=m; |R|=(2m+1)^3)"""
import sys, numpy as np
DC=0.025; EDGES=0.05+DC*np.arange(0,27)
def load(bc):
    S={};E=[]
    for l in open('../gbc/ev_As1e-2_bc%s.csv'%bc):
        p=l.strip().split(',')
        if p[0]=='S': S[int(p[1])]=float(p[2])
        elif p[0]=='E': E.append([float(x) for x in p[1:]])
    return S,np.array(E)
def est(seeds,lnWs,E,m,nmin=10):
    # E: seed,lnW,x,y,z,r,val ; seeds possibly resampled (list with repeats)
    cnt=dict(zip(*np.unique(seeds,return_counts=True)))
    w=np.array([cnt.get(int(s),0) for s in E[:,0]])
    inR=(np.abs(E[:,2:5]).max(1)<=m)&(w>0)
    N=len(seeds); nR=(2*m+1)**3; out=[]
    for a,b in zip(EDGES[:-1],EDGES[1:]):
        s=inR&(E[:,6]>=a)&(E[:,6]<b); k=w[s].sum()
        if s.sum()<nmin: out.append((0.5*(a+b),np.nan,int(s.sum())));continue
        lw=np.repeat(E[s,1],w[s]); Wb=np.exp(lw.mean()+lw.var(ddof=1)/2)
        out.append((0.5*(a+b),k/(N*nR)*Wb/DC,int(s.sum())))
    return out
if __name__=='__main__':
    outfn=sys.argv[1]; rows=[]; rng=np.random.default_rng(3)
    for spec in sys.argv[2:]:
        bc,m=spec.split(':'); m=int(m); S,E=load(bc); seeds=np.array(sorted(S)); lnWs=np.array([S[s] for s in seeds])
        base=est(seeds,lnWs,E,m)
        bs=np.array([[x[1] for x in est(rng.choice(seeds,len(seeds)),None,E,m)] for _ in range(100)])
        lo,hi=np.nanpercentile(bs,[16,84],axis=0)
        for (c,n,k),l,h in zip(base,lo,hi):
            if np.isfinite(n): rows.append((float(bc),m,c,n,l,h,k))
        print(bc,m,len(seeds),'samples')
    with open(outfn,'w') as f:
        f.write('bc,m,C,n,n_lo,n_hi,nev\n')
        for r in rows: f.write('%g,%d,%.4f,%.4e,%.4e,%.4e,%d\n'%r)
