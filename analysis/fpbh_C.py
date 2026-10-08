#!/usr/bin/env python3
"""Total fPBH and n(C), C-binned lognormal-W estimator, Laplacian-peak vs compaction-peak (site0 and k3 weights)."""
import sys, numpy as np
Cth, Gam, Mks, fac = 0.533, 0.36, 0.024079612549999338, 1/3.95e-20
KINDS = ['L-k3','L-site0','C-site0']

def prep(d, kind):
    pk, wt = kind.split('-')
    if pk == 'L': C, xm, z, pos = d[:,5], d[:,3], d[:,4], d[:,7:10]
    else:         C, xm, z, pos = d[:,10], d[:,11], d[:,12], d[:,13:16]
    w = d[:,2]**3 if wt == 'k3' else (np.abs(pos).sum(1) == 0).astype(float)
    sel = C > Cth
    M = np.where(sel, Mks*xm**2*np.exp(2*z)*np.clip(C-Cth,0,None)**Gam, 0.0)
    return C, M, w, sel

def est(d, kind, cb, nmin=10):
    C, M, w, sel = prep(d, kind); lnW = d[:,6]; N = len(d)
    ctot, nC = 0.0, []
    for a, b in zip(cb[:-1], cb[1:]):
        m = (C >= a) & (C < b)
        if m.sum() < nmin:
            nC.append(np.nan); continue
        Wb = np.exp(lnW[m].mean() + lnW[m].var(ddof=1)/2)
        nC.append(w[m].sum()/N*Wb/(b-a))
        if a >= Cth: ctot += (w[m]*M[m]).sum()/N*Wb*fac
    return ctot, np.array(nC)

if __name__ == '__main__':
    fn = sys.argv[1]; db = float(sys.argv[2]) if len(sys.argv) > 2 else 0.01
    d = np.loadtxt(fn, delimiter=','); d = d[np.isfinite(d).all(1)]
    cb = np.arange(Cth, 0.70, db)
    print('N =', len(d), ' bin width', db)
    rng = np.random.default_rng(1); res = {}
    for k in KINDS:
        t, nC = est(d, k, cb)
        bs = np.array([est(d[rng.integers(0,len(d),len(d))], k, cb)[0] for _ in range(200)])
        res[k] = (t, bs)
        print('%-8s fPBH_tot=%.3e  bootstrap 16/50/84%%: %.2e %.2e %.2e' % (k, t, *np.percentile(bs,[16,50,84])))
    for a, b in [('C-site0','L-site0'),('C-site0','L-k3'),('L-site0','L-k3')]:
        r = res[a][1]/res[b][1]
        print('ratio %s/%s = %.2f   (16-84%%: %.2f - %.2f)' % (a, b, res[a][0]/res[b][0], *np.percentile(r,[16,84])))
    # n(C) table
    print('  C     n_L(k3)    n_L(site0)  n_C(site0)')
    tabs = {k: est(d,k,cb)[1] for k in KINDS}
    for i, a in enumerate(cb[:-1]):
        print(' %.3f  %.2e   %.2e   %.2e' % (a+db/2, tabs['L-k3'][i], tabs['L-site0'][i], tabs['C-site0'][i]))
