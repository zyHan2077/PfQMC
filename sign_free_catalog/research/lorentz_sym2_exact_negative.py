#!/usr/bin/env python3
"""Exact rational negative control for the Lorentz traceless Sym^2 14-orbital model.

No numerical matrix exponentials and no SymPy are used.  A rational (not
orthonormal) basis is sufficient because the final determinant is basis invariant.
Allowed boost A=L_(E01+E10); relaxed B=Pi L_E00 Pi.  With a=log(4),
b=5 log(5/4), exp(aA) exp(bB) exp(-aA) exp(-bB) has a negative Fock trace.
The orthogonal projection Pi Y=Y-Tr(JY)J/5 makes B a self-adjoint physical
one-particle operator, but B does not belong to the Lorentz action algebra.
"""
from fractions import Fraction as F
import json

N = 5
def zero(n, m=None):
    return [[F(0) for _ in range(n if m is None else m)] for _ in range(n)]
def eye(n):
    x = zero(n)
    for i in range(n): x[i][i] = F(1)
    return x
def add(a,b): return [[x+y for x,y in zip(r,s)] for r,s in zip(a,b)]
def scale(a,c): return [[x*c for x in r] for r in a]
def mul(a,b):
    return [[sum((a[i][k]*b[k][j] for k in range(len(b))),F(0))
             for j in range(len(b[0]))] for i in range(len(a))]
def transpose(a): return list(map(list,zip(*a)))
def determinant(a):
    a=[r[:] for r in a]; d=F(1)
    for i in range(len(a)):
        pivot=next((j for j in range(i,len(a)) if a[j][i]),None)
        if pivot is None: return F(0)
        if pivot != i: a[i],a[pivot]=a[pivot],a[i]; d=-d
        p=a[i][i];d*=p
        for j in range(i+1,len(a)):
            t=a[j][i]/p
            for k in range(i+1,len(a)): a[j][k]-=t*a[i][k]
            a[j][i]=F(0)
    return d

J=eye(N)
for i in range(1,N): J[i][i]=-1
def project(y):
    trace=sum(J[i][i]*y[i][i] for i in range(N))
    return add(y,scale(J,-trace/N))
basis=[]
for i in range(1,N):
    t=zero(N);t[0][0]=t[i][i]=F(1);basis.append(t)
pairs=[(i,j) for i in range(N) for j in range(i+1,N)]
for i,j in pairs:
    t=zero(N);t[i][j]=t[j][i]=F(1);basis.append(t)
def coordinates(y): return [y[i][i] for i in range(1,N)]+[y[i][j] for i,j in pairs]
def action(k):
    return transpose([coordinates(project(add(mul(k,t),mul(t,transpose(k))))) for t in basis])
K=zero(N);K[0][1]=K[1][0]=F(1)
E=zero(N);E[0][0]=F(1)
A,B=action(K),action(E)
M=len(A);I=eye(M)
def spectral_function(a, nodes, values):
    out=zero(M)
    for node,value in zip(nodes,values):
        p=I
        for other in nodes:
            if other != node: p=mul(p,scale(add(a,scale(I,-other)),1/(node-other)))
        out=add(out,scale(p,value))
    return out
an=list(map(F,[-2,-1,0,1,2]));bn=[F(0),F(1),F(8,5)]
av=[F(4)**int(x) for x in an];bv=[F(1),F(5,4)**5,F(5,4)**8]
EA=spectral_function(A,an,av);EB=spectral_function(B,bn,bv)
EAi=spectral_function(A,an,[1/x for x in av]);EBi=spectral_function(B,bn,[1/x for x in bv])
assert mul(EA,EAi)==I and mul(EB,EBi)==I
g=mul(mul(mul(EA,EB),EAi),EBi)
z=determinant(add(I,g))
expected=F(-95232669932822733339285069632243,4611686018427387904000000000000)
assert z == expected
print(json.dumps({'orbitals':M,'signed_trace_exact':str(z),'signed_trace_float':float(z),
                  'a':'log(4)','b':'5 log(5/4)','inverse_checks':'exact',
                  'status':'PASS negative relaxed-J control'},indent=2))
