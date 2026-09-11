#!/usr/bin/env python3
"""Independent NumPy/SciPy screen of the full exceptional Albert action.

This is a numerical screening experiment, not a proof and not a Pfaffian test.
The repository's sign_free_albert executable supplies the separate signed-Pf test.
All 26 traceless Herm_3(O) directions participate in random time-dependent slices.
Example: OPENBLAS_NUM_THREADS=1 python albert_independent_screen.py --samples 500
"""
import argparse
import json
import numpy as np
from scipy.linalg import eigh, svdvals

def conjugate(a):
    b=a.copy();b[1:]*=-1;return b
def cd(a,b):
    if len(a)==1:return a*b
    h=len(a)//2;a0,a1=a[:h],a[h:];b0,b1=b[:h],b[h:]
    return np.r_[cd(a0,b0)-cd(conjugate(b1),a1),cd(b1,a0)+cd(a1,conjugate(b0))]
oct_basis=np.eye(8)
table=np.array([[cd(a,b) for b in oct_basis] for a in oct_basis])
def product(x,y):return np.einsum('ika,kjb,abc->ijc',x,y,table,optimize=True)
def jordan(x,y):return (product(x,y)+product(y,x))/2
def inner(x,y):return np.sum(x*y)
def determinant(x):
    a,b,c=x[0,0,0],x[1,1,0],x[2,2,0]
    triple=cd(cd(x[0,1],x[1,2]),x[2,0])[0]
    return a*b*c+2*triple-a*np.dot(x[1,2],x[1,2])-b*np.dot(x[0,2],x[0,2])-c*np.dot(x[0,1],x[0,1])
ts=[]
for i in range(3):
    t=np.zeros((3,3,8));t[i,i,0]=1;ts.append(t)
for i,j in [(0,1),(0,2),(1,2)]:
    for a in range(8):
        t=np.zeros((3,3,8));t[i,j,a]=1/np.sqrt(2);t[j,i]=conjugate(t[i,j]);ts.append(t)
ts=np.array(ts)
ks=np.array([(ts[0]-ts[1])/np.sqrt(2),(ts[0]+ts[1]-2*ts[2])/np.sqrt(6),*ts[3:]])
def action(k):return np.array([[inner(a,jordan(k,b)) for b in ts] for a in ts])
ls=np.array([action(k) for k in ks]);I=np.eye(27)
def exp_symmetric(h):
    w,v=eigh(h);return (v*np.exp(w))@v.T
def decode(v):return np.einsum('a,aijk->ijk',v,ts)
def encode(x):return np.einsum('aijk,ijk->a',ts,x)
def exp_jordan(k):
    power=np.zeros((3,3,8));power[range(3),range(3),0]=1
    s=max(0,int(np.ceil(np.log2(max(1.,np.linalg.norm(k)*2)))))
    k=k/2**s;out=power.copy()
    for n in range(1,81):
        power=jordan(k,power)/n;out+=power
        if np.linalg.norm(power)<2e-17:break
    else:raise AssertionError('Jordan exponential did not converge')
    for _ in range(s):out=jordan(out,out)
    return out
def exp_action_oracle(k):
    y=exp_jordan(k/2);yy=jordan(y,y)
    return np.array([encode(2*jordan(y,jordan(y,t))-jordan(yy,t)) for t in ts]).T

p=argparse.ArgumentParser();p.add_argument('--samples',type=int,default=500)
p.add_argument('--pulses',type=int,default=12);p.add_argument('--invariance',type=int,default=1000)
p.add_argument('--seed',type=int,default=9120401);args=p.parse_args()
rng=np.random.default_rng(args.seed)
assert np.max(np.abs(ls-ls.transpose(0,2,1)))<1e-14
unit=np.zeros((3,3,8));unit[range(3),range(3),0]=1
identity_error=float(np.max(np.abs(action(unit)-I)))
gram=ts.reshape(27,-1)@ts.reshape(27,-1).T
assert np.max(np.abs(gram-I))<1e-14 and identity_error<1e-14
comm=np.array([(ls[i]@ls[j]-ls[j]@ls[i]).ravel() for i in range(26) for j in range(i)])
sing=svdvals(comm);rank=int(np.sum(sing>1e-10))
assert rank==52
jordan_error=cubic_error=oracle_error=0.
for n in range(args.invariance):
    x=decode(rng.normal(size=27));y=decode(rng.normal(size=27));xx=jordan(x,x)
    e=np.linalg.norm(jordan(xx,jordan(x,y))-jordan(x,jordan(xx,y)))
    jordan_error=max(jordan_error,float(e))
    coeff=rng.uniform(-.15,.15,26);k=np.einsum('a,aijk->ijk',coeff,ks)
    e=exp_symmetric(np.einsum('a,aij->ij',coeff,ls));out=decode(e@encode(x))
    error=abs(determinant(out)-determinant(x))/max(1.,abs(determinant(x)))
    cubic_error=max(cubic_error,float(error))
    if n<24:oracle_error=max(oracle_error,float(np.linalg.norm(e-exp_action_oracle(k))))
min_logdet=float('inf');max_logdet=-float('inf');negative=zero=0
for n in range(args.samples):
    b=I.copy();logR=0.
    for step in range(args.pulses):
        coeff=rng.uniform(-.25,.25,26)
        b=exp_symmetric(np.einsum('a,aij->ij',coeff,ls))@b
        logR+=rng.uniform(-.3,.3)
    sign,logdet=np.linalg.slogdet(I+np.exp(logR)*b)
    negative+=int(sign<0);zero+=int(sign==0)
    min_logdet=min(min_logdet,float(logdet));max_logdet=max(max_logdet,float(logdet))
assert negative==0 and zero==0 and cubic_error<1e-9 and oracle_error<1e-10
print(json.dumps({'status':'PASS numerical screen (not proof)','seed':args.seed,'samples':args.samples,
 'pulses':args.pulses,'directions':26,'commutator_span_numerical_rank':rank,
 'negative_samples':negative,'zero_samples':zero,'min_log_det':min_logdet,'max_log_det':max_logdet,
 'identity_error':identity_error,'max_jordan_identity_absolute_error':jordan_error,
 'max_cubic_relative_error':cubic_error,'max_independent_exponential_error':oracle_error},indent=2))
