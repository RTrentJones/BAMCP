# Round-4 power, revised: primary = mean over 4 pos/neg pairs of balanced accuracy (anchors excluded).
# Stratum-specific baselines and effects; nested arms (rho=.5, q=.3); 3 models offsets -0.9/0/+0.9;
# item sd 1.2, item x cond sd 1.0. Test: sign-flip on item diffs weighted by 1/(2*n_class*4).
# Also cluster design effect: items in clusters of size 2 sharing a cluster effect (sd 0.6).
import numpy as np, math, sys
rng=np.random.default_rng(44)
expit=lambda x:1/(1+np.exp(-x)); logit=lambda p:math.log(p/(1-p))
beta=np.array([-0.9,0.0,0.9]); su,sv=1.2,1.0
# marginal->conditional intercept per stratum (solve numerically)
z=rng.standard_normal(100000)*math.sqrt(su**2+sv**2+0.36)
def mu_for(p):
    lo,hi=-8,8
    for _ in range(40):
        m=(lo+hi)/2; q=np.mean([expit(m+b+z).mean() for b in beta]); lo,hi=(m,hi) if q<p else (lo,m)
    return (lo+hi)/2
B={1:.62,2:.68,3:.72,4:.66,5:.68,6:.62,7:.70,8:.82}   # revised expected B (section 3)
pairs=[(1,2),(3,4),(5,6),(7,8)]
SCEN={'null':{},
 'uniform+0.10':{k:.10 for k in B},
 'C1-predicted':{1:0,2:.15,3:0,4:.15,5:-.10,6:.05,7:.10,8:-.05},
 'C2-predicted':{1:.05,2:.05,3:0,4:0,5:0,6:.15,7:.15,8:-.05},
 'concentrated+0.15(2 pairs)':{1:.15,2:.15,3:.15,4:.15}}
MU0={k:mu_for(B[k]) for k in B}
def run(nper,scen,nmod=3,nsim=2000,nperm=1000,alpha=0.025,margin=0.07):
    eff=SCEN[scen]; MU1={k:mu_for(min(.97,max(.03,B[k]+eff.get(k,0)))) for k in B}
    strata=np.repeat(list(B),nper); n=len(strata)
    w=np.ones(n)/(2*nper*4)  # balanced acc weight: per class mean, per pair mean, mean over pairs
    bb=beta if nmod==3 else beta[1:]
    m0=np.array([MU0[s] for s in strata]); m1=np.array([MU1[s] for s in strata])
    true=np.mean([ (expit(0)*0) for _ in [0]])
    sup=0; eq=0; est=[]
    for _ in range(nsim):
        cl=rng.standard_normal(n//2).repeat(2)*0.6; rng.shuffle(cl)
        u=rng.standard_normal(n)*su+cl
        e0=rng.standard_normal(n); e1=.5*e0+math.sqrt(.75)*rng.standard_normal(n)
        P0=expit(m0[:,None]+u[:,None]+sv*e0[:,None]+bb); P1=expit(m1[:,None]+u[:,None]+sv*e1[:,None]+bb)
        U0=rng.random(P0.shape); U1=np.where(rng.random(P0.shape)<.3,U0,rng.random(P0.shape))
        d=((U1<P1).mean(1)-(U0<P0).mean(1))*w
        obs=d.sum(); est.append(obs)
        s=rng.choice([-1,1],size=(nperm,n)); null=(s*d).sum(1)
        p=(np.sum(np.abs(null)>=abs(obs)-1e-12)+1)/(nperm+1)
        sup+= (p<alpha) and obs>0
        se=null.std(); eq+= (obs-1.96*se>-margin) and (obs+1.96*se<margin)
    est=np.array(est)
    return sup/nsim, eq/nsim, est.mean(), est.std()
for nper,nmod in ((15,3),(17,3),(20,3),(22,2)):
    for sc in SCEN:
        p,e,m,s=run(nper,sc,nmod)
        print(f"items/stratum={nper} pair_items={nper*8} models={nmod} runs={nper*8*3*nmod}+anchors scen={sc:28s} meanDelta={m:+.3f} SE={s:.3f} power={p:.2f} TOST(+-0.07)={e:.2f}",flush=True)
