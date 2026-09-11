import numpy as np
from scipy.linalg import expm
from itertools import product

# Independent dense references, no repository imports.
def pfaff(a):
    a=np.array(a,complex,copy=True); z=1+0j
    for k in range(0,len(a),2):
        p=k+1+np.argmax(abs(a[k,k+1:])); v=a[k,p]
        if abs(v)<1e-13: return 0j
        if p!=k+1:
            a[[p,k+1],:]=a[[k+1,p],:]; a[:,[p,k+1]]=a[:,[k+1,p]]; z=-z
        v=a[k,k+1]; z*=v
        x=a[k,k+2:].copy(); y=a[k+1,k+2:].copy()
        a[k+2:,k+2:]+=(np.outer(y,x)-np.outer(x,y))/v
    return z

def number_a(h):
    n=len(h); a=np.zeros((2*n,2*n),complex)
    for i in range(n):
        for j in range(n):
            a[2*i,2*j]=a[2*i+1,2*j+1]=.5*(h[i,j]-h[j,i])
            a[2*i,2*j+1]=.5j*(h[i,j]+h[j,i]); a[2*i+1,2*j]=-.5j*(h[i,j]+h[j,i])
    return a

def gauss(a):
    d=len(a); n=d//2; eye=np.eye(d); b=expm(a)
    g=2*np.linalg.inv(eye+b)-eye; g=(g-g.T)/2
    s=(expm(-a/4)-expm(a/4))/np.sqrt(2)
    tr=2**n*(-1)**n*pfaff(np.block([[s,-eye],[eye,s]]))
    return tr,g

def trace_gauss(fs):
    n=len(fs[0][1])//2; d=2*n
    if len(fs)%2: fs=fs+[(2**n,np.zeros((d,d))) ]
    L=len(fs); a=np.zeros((d*L,d*L),complex); z=2**n*(-1)**(n*(L//2))
    for i,(tr,g) in enumerate(fs[::-1]):
        a[i*d:(i+1)*d,i*d:(i+1)*d]=g; z*=tr/2**n
        for j in range(i+1,L):
            a[i*d:(i+1)*d,j*d:(j+1)*d]=(-1)**(j-i)*np.eye(d)
            a[j*d:(j+1)*d,i*d:(i+1)*d]=-(-1)**(j-i)*np.eye(d)
    return z*pfaff(a)

def annih(n,i):
    a=np.zeros((2**n,2**n),complex)
    for v in range(2**n):
        if v>>i&1: a[v^(1<<i),v]=(-1)**((v&((1<<i)-1)).bit_count())
    return a

def number_fock(h):
    n=len(h); cs=[annih(n,i) for i in range(n)]
    out=-np.trace(h)/2*np.eye(2**n,dtype=complex)
    for i in range(n):
        for j in range(n): out+=h[i,j]*cs[i].conj().T@cs[j]
    return out

def fermion_weight(N,edges,events,beta=1.2,t=.8,V=.5,hvals=None,stagger=True):
    hvals=hvals if hvals is not None else [.6]*N
    hop=np.zeros((N,N),complex)
    for i,j in edges: hop[i,j]=hop[j,i]=t
    q=number_fock(hop); exact=np.eye(2**N,dtype=complex); fs=[]; scalar=1+0j; last=0.
    for tm,kind,loc in events+[(beta,'end',0)]:
        dt=tm-last; last=tm; exact=expm(dt*q)@exact; fs.append(gauss(number_a(dt*hop)))
        if kind=='end': break
        scalar*=-1
        if kind=='f': sites=edges[loc]; coeff=V
        else: sites=[loc]; coeff=hvals[loc]*((-1)**loc if stagger else 1)
        scalar*=coeff; exact=coeff*(-1)*exact
        for site in sites:
            nh=np.zeros((N,N),complex); nh[site,site]=1
            exact=number_fock(nh)@exact
            fs.extend([gauss(number_a(.5j*np.pi*nh))]*2); scalar*=(-.5j)
    return scalar*trace_gauss(fs), np.trace(exact)

def spin_z(v,i): return (((v>>(2*i))&1)+((v>>(2*i+1))&1)-1)/2

def spin_path(N,edges,events,v,beta=1.2,J=.7):
    old=v; path=[v]; amp=1.; last=0.
    for tm,kind,loc in events+[(beta,'end',0)]:
        energy=-J*sum(spin_z(v,i)*spin_z(v,j) for i,j in edges)
        amp*=np.exp(-(tm-last)*energy); last=tm
        if kind=='s':
            occ=(v>>(2*loc))&3
            if occ not in (0,3): return 0.,path
            v^=3<<(2*loc); amp*=.5; path.append(v)
    return (amp if v==old else 0.),path

def projector(n,v):
    g=np.zeros((2*n,2*n),complex)
    for i in range(n): g[2*i,2*i+1]=1j*(-1 if v>>i&1 else 1); g[2*i+1,2*i]=-g[2*i,2*i+1]
    return 1.,g

def rotation(n,site,create):
    i,j=2*site,2*site+1; u=np.zeros(2*n,complex); v=np.zeros(2*n,complex)
    u[2*i]=.5; u[2*i+1]=-.5j; v[2*j]=.5; v[2*j+1]=-.5j
    a=2*(np.outer(u,v)-np.outer(v,u)); a=np.pi/2*(a-a.conj().T)*(1 if create else -1)
    return gauss(a)

def spin_pf(N,path,wrong=False):
    fs=[]; n=2*N
    for ix,(a,b) in enumerate(zip(path[:-1],path[1:])):
        site=((a^b).bit_length()-1)//2
        create=((a>>(2*site))&3)==0
        if wrong and ix==0: create=not create
        fs += [projector(n,a),rotation(n,site,create)]
    return trace_gauss(fs)

cases=[('two_site_density_spin',2,[(0,1)],[(.13,'f',0),(.27,'s',0),(.49,'s',1),(.78,'s',0),(1.03,'s',1)],[.6,-.4]),
       ('four_site_square',4,[(0,1),(1,2),(2,3),(3,0)],[(.08,'f',0),(.23,'s',0),(.41,'s',1),(.63,'f',2),(.82,'s',0),(1.01,'s',1)],[.6,-.4,.3,.7])]
for name,N,ed,ev,hv in cases:
    wp,wf=fermion_weight(N,ed,ev,hvals=hv)
    zs=0.; count=0; maxerr=0.; chosen=None
    for v in range(4**N):
        amp,path=spin_path(N,ed,ev,v)
        if amp:
            zs+=amp; count+=1
            if count<=8:
                sp=spin_pf(N,path); maxerr=max(maxerr,abs(sp-1))
                chosen=path
    print(name,'fermion_pf=',wp,'fock=',wf,'spin_sum=',zs,'joint_weight=',wp*zs,'nonzero_spin_paths=',count,'spin_pf_max_error=',maxerr)
    if chosen: print('signed_rotation_negative_control=',spin_pf(N,chosen,True))
ev=[(.3,'s',0),(.7,'s',1)]
wp,wf=fermion_weight(2,[(0,1)],ev,hvals=[.6,.4],stagger=False)
correct,_=fermion_weight(2,[(0,1)],ev,hvals=[.6,.4],stagger=True)
print('raw_density_negative_control=',wp,'fock=',wf,'staggered_fermion=',correct,'spin_sum=',sum(spin_path(2,[(0,1)],ev,v)[0] for v in range(16)))
# Independent original 6-mode Fock construction and rotated CT-INT word.
N=2; n=6; ed=[(0,1)]; t=.8; V=.5; J=.7; hv=[.6,-.4]; beta=1.2
cs=[annih(n,i) for i in range(n)]; eye=np.eye(2**n); ns=[c.conj().T@c for c in cs]
dn=[ns[i]-.5*eye for i in range(N)]; sz=[]; sx=[]
for i in range(N):
    a,b=cs[N+2*i],cs[N+2*i+1]
    sz.append((ns[N+2*i]+ns[N+2*i+1]-eye)/2)
    sx.append((a.conj().T@b.conj().T+b@a)/2)
H0=-t*(cs[0].conj().T@cs[1]+cs[1].conj().T@cs[0])-J*sz[0]@sz[1]
H=H0+V*dn[0]@dn[1]+sum(hv[i]*dn[i]@sx[i] for i in range(N))
HU=H0+V*dn[0]@dn[1]+sum(hv[i]*(-1)**i*dn[i]@sx[i] for i in range(N))
U=expm(1j*np.pi*sz[1]); print('full_fock_H_hermiticity=',np.linalg.norm(H-H.conj().T),'rotation_error=',np.linalg.norm(U.conj().T@H@U-HU))
ev=cases[0][3]; out=eye.astype(complex); last=0.
for tm,kind,loc in ev+[(beta,'end',0)]:
    out=expm(-(tm-last)*H0)@out; last=tm
    if kind=='f': out=(-V*dn[0]@dn[1])@out
    elif kind=='s': out=(-hv[loc]*(-1)**loc*dn[loc]@sx[loc])@out
print('full_fock_word_trace=',np.trace(out))
def embed_factor(f,offset,total):
    tr,g=f; n=len(g)//2; G=np.zeros((2*total,2*total),complex)
    G[2*offset:2*(offset+n),2*offset:2*(offset+n)]=g
    return tr*2**(total-n),G

def global_path_weight(N,ed,ev,v,hv,beta=1.2,t=.8,V=.5,J=.7):
    total=3*N; hop=np.zeros((N,N),complex)
    for i,j in ed: hop[i,j]=hop[j,i]=t
    fs=[embed_factor(projector(2*N,v),N,total)]; scalar=1+0j; last=0.; start=v
    for tm,kind,loc in ev+[(beta,'end',0)]:
        dt=tm-last; last=tm
        scalar*=np.exp(dt*J*sum(spin_z(v,i)*spin_z(v,j) for i,j in ed))
        fs.append(embed_factor(gauss(number_a(dt*hop)),0,total))
        if kind=='end': break
        scalar*=-1
        if kind=='f': sites=ed[loc]; scalar*=V
        else: sites=[loc]; scalar*=hv[loc]*(-1)**loc*.5
        for site in sites:
            nh=np.zeros((N,N),complex); nh[site,site]=1
            fac=embed_factor(gauss(number_a(.5j*np.pi*nh)),0,total)
            fs.extend([fac,fac]); scalar*=(-.5j)
        if kind=='s':
            occ=(v>>(2*loc))&3
            if occ not in (0,3): return 0j
            fs.append(embed_factor(projector(2*N,v),N,total))
            fs.append(embed_factor(rotation(2*N,loc,occ==0),N,total)); v^=3<<(2*loc)
    if v!=start:return 0j
    return scalar*trace_gauss(fs)
name,N,ed,ev,hv=cases[0]
wp,_=fermion_weight(N,ed,ev,hvals=hv); amp,path=spin_path(N,ed,ev,0)
print('single_combined_full_6mode_Pf=',global_path_weight(N,ed,ev,0,hv),'factorized_path=',wp*amp)
rng=np.random.default_rng(160507420); positive=0; zero=0; largest=0.; minpos=np.inf
for N,ed in [(2,[(0,1)]),(4,[(0,1),(1,2),(2,3),(3,0)])]:
    for sample in range(16):
        hv=rng.uniform(.2,.9,N)*rng.choice([-1,1],N)
        vv=[('f',int(rng.integers(len(ed)))) for _ in range(int(rng.integers(0,4)))]
        for _ in range(int(rng.integers(1,4))):
            site=int(rng.integers(N)); vv.extend([('s',site),('s',site)])
        rng.shuffle(vv); times=np.sort(rng.uniform(0,1.2,len(vv)))
        ev=[(float(tm),kind,loc) for tm,(kind,loc) in zip(times,vv)]
        pf,exact=fermion_weight(N,ed,ev,hvals=hv)
        zs=sum(spin_path(N,ed,ev,v)[0] for v in range(4**N))
        wp,we=pf*zs,exact*zs; largest=max(largest,abs(wp-we)/max(1.,abs(we)))
        assert abs(wp.imag)<1e-10 and wp.real>=-1e-10
        if we.real>1e-12: positive+=1; minpos=min(minpos,we.real); assert abs(wp/abs(wp)-1)<1e-8
        else: zero+=1
print('random_audit_seed=160507420 count=32 positive=',positive,'small_or_zero=',zero,'min_positive=',minpos,'max_absolute_scaled_Pf_Fock_error=',largest)
