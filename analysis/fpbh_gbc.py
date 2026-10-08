#!/usr/bin/env python3
"""fPBH from IS samples (Gaussian_LN_GBC): Laplacian-peak vs compaction-peak definition.
cols: 0 seed,1 mu2,2 k3,3 k*rm(L),4 zeta(L),5 C(L),6 lnW,7-9 L pos,10 Cpk,11 k*rm(pk),12 zeta(pk),13-15 C-peak pos
Estimators (same per-bin lognormal W-bar as fPBH_IS.nb):
  n_bin = (sum_i w_i over samples in bin)/(N dlnM) * Wbar_bin,  Wbar=exp(mean lnW + var/2)
  w_i = k3^3            (existing IS estimator)             -> 'k3'
  w_i = 1[peak at origin site]   (exact single-site identity) -> 'site0'
fPBH(M)=M n / 3.95e-20 (M in 1e20 g), M=Mks*xm^2 exp(2 zeta)(C-Cth)^0.36
"""
import sys, numpy as np
Cth, Gam, Mks, fac = 0.533, 0.36, 0.024079612549999338, 1/3.95e-20

def load(fn):
    return np.loadtxt(fn, delimiter=',')

def mass(xm, z, C):
    return Mks * xm**2 * np.exp(2*z) * (C-Cth)**Gam

def fpbh(d, kind, dlnM=0.5, lo=None, hi=None, idx=None):
    """returns (lnM centres, f(M)=dfPBH/dlnM, total fPBH) for peak kind in {'L','C'} and weight in kind string"""
    pk, wt = kind.split('-')
    N = len(d)
    if pk == 'L':
        C, xm, z, pos = d[:,5], d[:,3], d[:,4], d[:,7:10]
    else:
        C, xm, z, pos = d[:,10], d[:,11], d[:,12], d[:,13:16]
    if wt == 'k3':
        w = d[:,2]**3
    elif wt == 'site0':
        w = (np.abs(pos).sum(axis=1) == 0).astype(float)
    sel = (C > Cth)
    lnM = np.log(mass(xm, z, C))
    lnW = d[:,6]
    return C, lnM, w, lnW, sel, N

def estimate(d, kind, edges):
    C, lnM, w, lnW, sel, N = fpbh(d, kind)
    ctr, f = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        m = sel & (lnM >= a) & (lnM < b)
        if m.sum() < 2: 
            ctr.append(0.5*(a+b)); f.append(0.0); continue
        Wbar = np.exp(lnW[m].mean() + lnW[m].var(ddof=1)/2)
        nb = w[m].sum()/(N*(b-a)) * Wbar
        M = np.exp(0.5*(a+b))
        ctr.append(0.5*(a+b)); f.append(M*nb*fac)
    ctr, f = np.array(ctr), np.array(f)
    return ctr, f, (f*np.diff(edges)).sum()

if __name__ == '__main__':
    d = load(sys.argv[1] if len(sys.argv)>1 else '../gbc/muk_bc12.5.csv')
    print('N =', len(d))
    allM = []
    for k in ['L-k3','L-site0','C-site0']:
        C, lnM, w, lnW, sel, N = fpbh(d, k)
        allM.append(lnM[sel])
    allM = np.concatenate(allM)
    edges = np.arange(np.floor(allM.min()), np.ceil(allM.max())+0.5, 0.5)
    for k in ['L-k3','L-site0','C-site0']:
        c, f, tot = estimate(d, k, edges)
        print(k, 'fPBH_tot = %.4e' % tot)

def run_all(fn, nboot=300, dlnM=0.5, seed=1):
    d = load(fn); N = len(d)
    kinds = ['L-k3','L-site0','C-site0']
    allM = np.concatenate([fpbh(d,k)[1][fpbh(d,k)[4]] for k in kinds])
    edges = np.arange(np.floor(allM.min()), np.ceil(allM.max())+dlnM, dlnM)
    res = {k: estimate(d, k, edges) for k in kinds}
    rng = np.random.default_rng(seed)
    tot = {k: [] for k in kinds}
    for b in range(nboot):
        db = d[rng.integers(0, N, N)]
        for k in kinds:
            tot[k].append(estimate(db, k, edges)[2])
    tot = {k: np.array(v) for k, v in tot.items()}
    return d, edges, res, tot

if __name__ == '__main__' and len(sys.argv) > 2 and sys.argv[2] == 'full':
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    d, edges, res, tot = run_all(sys.argv[1])
    print('N =', len(d))
    for k in res:
        t = tot[k]; print('%-8s fPBH_tot = %.3e  (bootstrap 16-50-84%%: %.2e %.2e %.2e)' % (k, res[k][2], *np.percentile(t,[16,50,84])))
    r = tot['C-site0']/tot['L-site0']
    print('ratio C-site0 / L-site0 : %.2f  (16-84%%: %.2f - %.2f)' % (res['C-site0'][2]/res['L-site0'][2], *np.percentile(r,[16,84])))
    r2 = tot['C-site0']/tot['L-k3']
    print('ratio C-site0 / L-k3    : %.2f  (16-84%%: %.2f - %.2f)' % (res['C-site0'][2]/res['L-k3'][2], *np.percentile(r2,[16,84])))
    r3 = tot['L-site0']/tot['L-k3']
    print('ratio L-site0 / L-k3    : %.2f  (16-84%%: %.2f - %.2f)' % (res['L-site0'][2]/res['L-k3'][2], *np.percentile(r3,[16,84])))
    fig, ax = plt.subplots(figsize=(6,4.2))
    sty = {'L-k3':('k','--','Laplacian peak (k3$^3$ estimator; existing IS)'),
           'L-site0':('b','-','Laplacian peak (site-0 estimator)'),
           'C-site0':('r','-','compaction peak (site-0 estimator)')}
    for k,(c,ls,lab) in sty.items():
        ctr, f, t = res[k]; m = f>0
        ax.plot(np.exp(ctr[m]), f[m], color=c, ls=ls, marker='o', ms=3, label=lab)
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel(r'$M_{\rm PBH}$ [$10^{20}$ g]'); ax.set_ylabel(r'$df_{\rm PBH}/d\ln M$')
    ax.legend(fontsize=7); fig.tight_layout(); fig.savefig('fpbh_gbc.png', dpi=150)
    np.savetxt('fpbh_gbc_curves.csv', np.column_stack([np.exp(res['L-k3'][0]), res['L-k3'][1], res['L-site0'][1], res['C-site0'][1]]),
               delimiter=',', header='M[1e20g],f_L_k3,f_L_site0,f_C_site0', comments='')
