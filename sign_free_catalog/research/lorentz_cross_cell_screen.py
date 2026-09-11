#!/usr/bin/env python3
"""Reject a fixed-sign intercell hopping added to independently driven Lorentz cells.
This is a determinant counterexample search, not an accepted sign-free model.
"""
import json
import numpy as np
from scipy.linalg import eigh
from decimal import Decimal as D, getcontext
def dz(n):return [[D(0) for _ in range(n)] for _ in range(n)]
def di(n):
    a=dz(n)
    for i in range(n):a[i][i]=D(1)
    return a
def dm(a,b):return [[sum((a[i][k]*b[k][j] for k in range(len(a))),D(0)) for j in range(len(a))] for i in range(len(a))]
def dd(a):
    a=[r[:] for r in a];out=D(1)
    for i in range(len(a)):
        p=max(range(i,len(a)),key=lambda j:abs(a[j][i]))
        if p!=i:a[i],a[p]=a[p],a[i];out=-out
        out*=a[i][i]
        for j in range(i+1,len(a)):
            t=a[j][i]/a[i][i]
            for k in range(i+1,len(a)):a[j][k]-=t*a[i][k]
    return out
rng=np.random.default_rng(9121101)
n=8
T=np.zeros((n,n));T[0,4]=T[4,0]=1
def exp(h):
    w,v=eigh(h);return (v*np.exp(w))@v.T
fixed=exp(.6*T)
for sample in range(10000):
    b=np.eye(n);fields=[]
    for step in range(12):
        h=np.zeros((n,n));co=rng.normal(0,.8,size=(2,4));fields.append(co.tolist())
        for r in range(2):
            for a in range(3):h[4*r,4*r+a+1]=h[4*r+a+1,4*r]=co[r,a]
            h[4*r:4*r+4,4*r:4*r+4]+=co[r,3]*np.eye(4)
        b=fixed@exp(h)@b
    sign,ld=np.linalg.slogdet(np.eye(n)+b)
    if sign<0:
        getcontext().prec=70
        e=D('.6').exp();ei=1/e;fm=di(n);fm[0][0]=fm[4][4]=(e+ei)/2
        fm[0][4]=fm[4][0]=(e-ei)/2
        bm=di(n);sum_mu=D(0)
        for co in fields:
            em=di(n)
            for r in range(2):
                v=[D(str(x)) for x in co[r][:3]];mu=D(str(co[r][3]));sum_mu+=mu
                norm=sum(x*x for x in v).sqrt();s=dz(4)
                for a in range(3):s[0][a+1]=s[a+1][0]=v[a]
                s2=dm(s,s);en=norm.exp();eni=1/en;f1=(en-eni)/(2*norm);f2=((en+eni)/2-1)/(norm*norm)
                block=[[mu.exp()*(D(i==j)+f1*s[i][j]+f2*s2[i][j]) for j in range(4)] for i in range(4)]
                for i in range(4):
                    for j in range(4):em[4*r+i][4*r+j]=block[i][j]
            bm=dm(dm(fm,em),bm)
        det70=dd([[bm[i][j]+D(i==j) for j in range(n)] for i in range(n)]);centered=(-2*sum_mu).exp()*det70
        assert det70<0
        print(json.dumps({'status':'negative determinant counterexample','seed':9121101,
          'sample':sample,'fixed_center_center_hopping':.6,'slices':12,
          'log_abs_det':float(ld),'determinant_70_digit':str(det70),
          'centered_fock_trace_70_digit':str(centered),
          'fields_boost1_boost2_boost3_density':fields},indent=2))
        break
else:print(json.dumps({'status':'no negative found (not proof)','samples':10000}))
