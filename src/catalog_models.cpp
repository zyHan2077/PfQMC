#include "catalog/catalog.h"
#include <map>
#include <numeric>

namespace catalog {
static MatType identity(int n) { return MatType::Identity(n,n); }
static MatType zeroA(int n) { return MatType::Zero(2*n,2*n); }
static MatType density(int n,int i) { MatType h=MatType::Zero(n,n); h(i,i)=1.; return numberMajorana(h); }
static MatType hop(int n,int i,int j,DataType value=1.) {
    MatType h=MatType::Zero(n,n); h(i,j)=value; h(j,i)=std::conj(value); return numberMajorana(h);
}
static MatType embed(const MatType& a,int modes,int offset) {
    MatType out=zeroA(modes); out.block(2*offset,2*offset,a.rows(),a.cols())=a; return out;
}
static std::vector<std::pair<int,int>> grid(int lx,int ly) {
    require(lx>0 && ly>0,"Invalid lattice dimensions");
    std::vector<std::pair<int,int>> edges;
    for(int x=0;x<lx;++x) for(int y=0;y<ly;++y) {
        int r=x*ly+y;
        if(x+1<lx) edges.emplace_back(r,(x+1)*ly+y);
        if(y+1<ly) edges.emplace_back(r,r+1);
    }
    return edges;
}
Hamiltonian openHopping(int n) {
    require(n>=2,"Open hopping requires at least two modes");
    Hamiltonian h; h.id="open_hopping"; h.modes=n; h.q0=zeroA(n);
    for(int i=0;i<n;++i) {
        h.q0+=(.13-.05*i)*density(n,i);
        h.channels.push_back({"density_"+std::to_string(i),density(n,i),.15});
    }
    for(int i=0;i<n-1;++i) {
        h.q0+=.4*hop(n,i,i+1,std::exp(.19i*(i+1.)));
        h.channels.push_back({"density_bond_"+std::to_string(i),density(n,i)-density(n,i+1),.2});
    }
    h.description="Open hopping chain with fixed edge phases and complete negative density squares";
    h.validate(); return h;
}
Hamiltonian oddPairRing(int n,int pairs) {
    require(n>=3 && pairs>=1 && pairs<=n && pairs%2==1,"Ring needs an odd number of pairing edges");
    Hamiltonian h; h.id="odd_pair_ring"; h.modes=n; h.q0=zeroA(n); h.conserving=false;
    for(int e=0;e<n;++e) {
        int i=e,j=(e+1)%n; if(i>j) std::swap(i,j);
        DataType phase=std::exp(.17i*(e+1.));
        bool paired=(e==n-1 || (e<pairs-1));
        h.q0+=.3*(paired ? pairMajorana(n,i,j,phase):hop(n,i,j,phase));
        h.q0+=.04*(e-1)*density(n,e);
        h.channels.push_back({"density_bond_"+std::to_string(e),density(n,i)+((e%2) ? -1.:1.)*density(n,j),.18});
    }
    h.description="One fixed channel per edge; odd pairing parity; shared real density HS fields";
    h.validate(); return h;
}
static std::vector<MatType> fourBasis() {
    return {density(4,0)+density(4,1),density(4,1)+density(4,2),
            hop(4,1,2)-hop(4,3,0)+pairMajorana(4,3,0,1.),
            -hop(4,0,1)+hop(4,3,2,std::exp(-1.i*std::acos(-1.)/3.))
                -pairMajorana(4,3,2,std::exp(-1.i*std::acos(-1.)/3.))};
}
Hamiltonian fourMode() {
    auto b=fourBasis(); Hamiltonian h; h.id="four_mode"; h.modes=4; h.conserving=false;
    h.q0=.2*b[0]-.15*b[1]+.4*b[2]+.3*b[3];
    h.channels={{"O1",b[0],.7},{"O2",b[1],1.1}}; h.scalarEnergy=(.7+1.1)/4.;
    h.description="Explicit four-mode Hermitian pairing model of classification 3.2";
    h.validate(); return h;
}
static std::vector<MatType> hermitianTracelessBasis(int d) {
    std::vector<MatType> b;
    for(int k=1;k<d;++k) {
        MatType t=MatType::Zero(d,d); double c=1./std::sqrt(double(k*(k+1)));
        for(int i=0;i<k;++i) t(i,i)=c; t(k,k)=-k*c; b.push_back(t);
    }
    for(int i=0;i<d;++i) for(int j=i+1;j<d;++j) {
        MatType t=MatType::Zero(d,d); t(i,j)=t(j,i)=1./std::sqrt(2.); b.push_back(t);
        t(i,j)=-1.i/std::sqrt(2.); t(j,i)=1.i/std::sqrt(2.); b.push_back(t);
    }
    return b;
}
Hamiltonian adjointOdd(int d) {
    require(d>=3 && d%2==1,"Adjoint pairing needs odd d >= 3");
    auto basis=hermitianTracelessBasis(d);
    auto ad=[&](const MatType& x) {
        int size=d*d-1; MatType a(size,size);
        for(int j=0;j<size;++j) { MatType comm=x*basis[j]-basis[j]*x;
            for(int i=0;i<size;++i) a(i,j)=(basis[i]*comm).trace(); }
        return a;
    };
    Hamiltonian h; h.id="adjoint_d"+std::to_string(d); h.modes=(d*d-1)/2; h.conserving=false;
    MatType x=MatType::Zero(d,d);
    for(int i=0;i<d;++i) x(i,i)=.08*(i-.5*(d-1));
    for(int i=0;i<d-1;++i) x(i,i+1)=x(i+1,i)=.3;
    h.q0=ad(x);
    for(int i=0;i<d-1;++i) h.channels.push_back({"Cartan_"+std::to_string(i),ad(basis[i]),.2});
    h.description="Odd adjoint representation; only diagonal real HS fields, fixed positive adjacent auxiliary edges";
    h.validate(); return h;
}
MatType spinX(int m) {
    require(m>=0,"Negative spin"); MatType x=MatType::Zero(m+1,m+1);
    for(int i=0;i<m;++i) x(i,i+1)=x(i+1,i)=.5*std::sqrt(double((m-i)*(i+1)));
    return x;
}
MatType spinZ(int m) {
    MatType z=MatType::Zero(m+1,m+1); for(int i=0;i<=m;++i) z(i,i)=.5*m-i; return z;
}
Hamiltonian spinCell(int m,bool mu) {
    require(m>=0 && ((!mu && m%4==3) || m%2==0),"This spin/fugacity combination is not sign-free");
    Hamiltonian h; h.id="spin_cell_m"+std::to_string(m); h.modes=m+1;
    auto x=numberMajorana(spinX(m)),z=numberMajorana(spinZ(m)),d=numberMajorana(identity(m+1));
    h.q0=.23*x-.14*z+(mu ? .11:0.)*d;
    h.channels={{"Sx",x,.2},{"Sz",z,.3}};
    if(mu) h.channels.push_back({"D",d,.25});
    h.description="Physical orbital spin representation with unprojected full Fock space";
    h.validate(); return h;
}
static MatType kron(const MatType& a,const MatType& b) {
    MatType out(a.rows()*b.rows(),a.cols()*b.cols());
    for(int i=0;i<a.rows();++i) for(int j=0;j<a.cols();++j) out.block(i*b.rows(),j*b.cols(),b.rows(),b.cols())=a(i,j)*b;
    return out;
}
Hamiltonian spinTensor(int jx,int jy) {
    require(jx>=1 && jy>=1,"Integer positive tensor spins required");
    int nx=2*jx+1,ny=2*jy+1,n=nx*ny;
    std::vector<MatType> a={numberMajorana(identity(n)),numberMajorana(kron(spinX(2*jx),identity(ny))),
        numberMajorana(kron(spinZ(2*jx),identity(ny))),numberMajorana(kron(identity(nx),spinX(2*jy))),
        numberMajorana(kron(identity(nx),spinZ(2*jy)))};
    Hamiltonian h; h.id="spin_tensor_"+std::to_string(nx)+"x"+std::to_string(ny); h.modes=n;
    h.q0=.11*a[0]+.2*a[1]+.17*a[2]-.12*a[3]+.13*a[4];
    for(int i=0;i<5;++i) h.channels.push_back({"global_"+std::to_string(i),a[i],.08});
    h.description="2D integer-spin tensor hopping with five correlated global fields; interactions need not be spatially local";
    h.validate(); return h;
}
Hamiltonian spinLattice(int j,int lx,int ly,bool mixed) {
    require(j>=1,"Integer spin j>=1 required"); int cells=lx*ly,q=2*j+1,n=cells*q;
    auto edges=grid(lx,ly);
    std::vector<MatType> local={numberMajorana(identity(q)),numberMajorana(spinX(2*j)),numberMajorana(spinZ(2*j))};
    Hamiltonian h; h.id=mixed ? "spin_lattice_mixed":"spin_lattice_j"+std::to_string(j); h.modes=n; h.q0=zeroA(n);
    for(int r=0;r<cells;++r) {
        for(int a=0;a<3;++a) {
            MatType o=embed(local[a],n,q*r);
            h.q0+=(a==0 ? .08*(r%2 ? -1.:1.):a==1 ? .17:-.12)*o;
            h.channels.push_back({"onsite_"+std::to_string(r)+"_"+std::to_string(a),o,.12+.02*a});
        }
    }
    int e=0;
    for(auto edge:edges) {
        for(int a=1;a<=2;++a) {
            int b=mixed ? 3-a:a; double eta=e%2 ? -1.:1.;
            h.channels.push_back({"edge_"+std::to_string(e)+"_"+std::to_string(a),
                embed(local[a],n,q*edge.first)+eta*embed(local[b],n,q*edge.second),.1});
        }
        ++e;
    }
    h.description="Finite-range negative squares, including onsite compensation and shared bond fields";
    h.validate(); return h;
}
Hamiltonian fiveOrbital(int lx,int ly) {
    int cells=lx*ly,n=5*cells; auto base=spinLattice(1,lx,ly); Hamiltonian h;
    h.id="five_orbital"; h.modes=n; h.q0=embed(base.q0,n,0);
    // Order all 3-orbital cells first, then the two identical itinerant blocks.
    for(const auto& c:base.channels) h.channels.push_back({c.name,embed(c.a,n,0),c.coupling});
    for(int r=0;r<cells;++r) {
        MatType c=density(n,3*cells+r)+density(n,4*cells+r);
        h.q0+=.09*c;
        MatType z=embed(numberMajorana(spinZ(2)),n,3*r);
        h.channels.push_back({"shared_itinerant_"+std::to_string(r),c+z,.25});
    }
    for(auto e:grid(lx,ly)) for(int spin=0;spin<2;++spin) h.q0+=.35*hop(n,(3+spin)*cells+e.first,(3+spin)*cells+e.second);
    h.description="Dynamical five-orbital local model: two itinerant copies share exactly the same HS field with Sz";
    h.validate(); return h;
}
Hamiltonian pairingLattice(int lx,int ly) {
    int cells=lx*ly,n=4*cells; auto base=fourMode(); auto b=fourBasis(); Hamiltonian h;
    h.id="pairing_lattice"; h.modes=n; h.q0=zeroA(n); h.conserving=false; h.scalarEnergy=cells*base.scalarEnergy;
    for(int r=0;r<cells;++r) {
        h.q0+=embed(base.q0,n,4*r);
        for(const auto& c:base.channels) h.channels.push_back({c.name+"_"+std::to_string(r),embed(c.a,n,4*r),c.coupling});
    }
    int i=0;
    for(auto e:grid(lx,ly)) {
        h.channels.push_back({"pairing_bond_"+std::to_string(i),embed(b[0],n,4*e.first)+(i%2 ? -1.:1.)*embed(b[1],n,4*e.second),.2}); ++i;
    }
    h.description="Array of four-mode pairing cells coupled through O1/O2 squares";
    h.validate(); return h;
}
static std::vector<unsigned> subsets(int d,int k) {
    require(d>0 && d<=16 && k>=0 && k<=d,"Exterior basis limits: d<=16, 0<=k<=d");
    std::vector<unsigned> out;
    for(unsigned b=0;b<(1u<<d);++b) if(__builtin_popcount(b)==k) out.push_back(b);
    return out;
}
MatType wedgeAction(const MatType& x,int k) {
    int d=x.rows(); auto basis=subsets(d,k); std::map<unsigned,int> index;
    for(int i=0;i<int(basis.size());++i) index[basis[i]]=i;
    MatType out=MatType::Zero(basis.size(),basis.size());
    for(int col=0;col<int(basis.size());++col) for(int j=0;j<d;++j) if(basis[col]&(1u<<j)) {
        unsigned removed=basis[col]^(1u<<j); int s1=__builtin_popcount(basis[col]&((1u<<j)-1));
        for(int i=0;i<d;++i) if(!(removed&(1u<<i))) {
            int s2=__builtin_popcount(removed&((1u<<i)-1));
            out(index.at(removed|(1u<<i)),col)+=((s1+s2)%2 ? -1.:1.)*x(i,j);
        }
    }
    return out;
}
MatType primitiveBasis(int n,int k) {
    require(k>=1 && k<=n,"Primitive exterior k must be in [1,n]");
    auto big=subsets(2*n,k); if(k<2) return identity(big.size());
    auto small=subsets(2*n,k-2); std::map<unsigned,int> rows;
    for(int i=0;i<int(small.size());++i) rows[small[i]]=i;
    Eigen::MatrixXd c=Eigen::MatrixXd::Zero(small.size(),big.size());
    for(int col=0;col<int(big.size());++col) for(int i=0;i<n;++i) {
        unsigned b=big[col]; if(!(b&(1u<<i)) || !(b&(1u<<(n+i)))) continue;
        int s1=__builtin_popcount(b&((1u<<i)-1)); b^=1u<<i;
        int s2=__builtin_popcount(b&((1u<<(n+i))-1)); b^=1u<<(n+i);
        c(rows.at(b),col)+=(s1+s2)%2 ? -1.:1.;
    }
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(c,Eigen::ComputeFullV);
    require(svd.rank()==int(small.size()),"Symplectic contraction rank mismatch");
    MatType p=svd.matrixV().rightCols(big.size()-small.size()).cast<DataType>();
    require((c.cast<DataType>()*p).norm()<1e-10,"Primitive contraction kernel failed");
    return p;
}
static MatType slSymmetric(int d,int which) {
    MatType a=MatType::Zero(d,d);
    for(int i=0;i<d;++i) for(int j=i;j<d;++j) a(i,j)=a(j,i)=std::sin((i+1.)*(j+2.)+which*.71)/std::sqrt(double(d));
    a-=a.trace()/double(d)*identity(d); return a;
}
Hamiltonian exteriorSL(int d,int k) {
    require(d==2*k && k>=2,"Only SL middle exterior powers k>=2 are globally sign-free");
    MatType a=wedgeAction(slSymmetric(d,0),k),b=wedgeAction(slSymmetric(d,1),k);
    Hamiltonian h; h.id="sl_exterior_"+std::to_string(d)+"_"+std::to_string(k); h.modes=a.rows();
    h.q0=numberMajorana(.14*a-.11*b);
    h.channels={{"rho_X",numberMajorana(a),.08},{"rho_Y",numberMajorana(b),.08}};
    h.description="One physical orbital per k-subset; traceless correlated fields, no additional fugacity";
    h.validate(); return h;
}
static long choose(int n,int k) { if(k<0 || k>n) return 0; long v=1; for(int i=1;i<=k;++i) v=v*(n-i+1)/i; return v; }
static MatType spSymmetric(int n,int seed) {
    MatType a=MatType::Zero(n,n),b=a;
    for(int i=0;i<n;++i) for(int j=i;j<n;++j) {
        a(i,j)=a(j,i)=std::sin((i+1.)*(j+2.)+seed)/std::sqrt(double(n));
        b(i,j)=b(j,i)=std::cos((i+2.)*(j+1.)+seed)/std::sqrt(double(n));
    }
    MatType x(2*n,2*n); x.topLeftCorner(n,n)=a; x.bottomRightCorner(n,n)=-a;
    x.topRightCorner(n,n)=x.bottomLeftCorner(n,n)=b; return x;
}
Hamiltonian exteriorSp(int n,int k) {
    require(k>=1 && k<=n && (k%2==0 || choose(n,(k-1)/2)%2==0),"Primitive symplectic power fails sign-free parity criterion");
    MatType p=primitiveBasis(n,k),a=p.adjoint()*wedgeAction(spSymmetric(n,0),k)*p,b=p.adjoint()*wedgeAction(spSymmetric(n,1),k)*p;
    Hamiltonian h; h.id="sp_exterior_"+std::to_string(n)+"_"+std::to_string(k); h.modes=p.cols();
    h.q0=numberMajorana(.12*a-.09*b);
    h.channels={{"rho_X",numberMajorana(a),.045},{"rho_Y",numberMajorana(b),.045}};
    h.description="Physical orbitals span the symplectic contraction kernel; no arbitrary chemical potential";
    h.validate(); return h;
}
Hamiltonian sl2Orthogonal(const std::vector<int>& mult) {
    int n=0,parity=0;
    for(int m=0;m<int(mult.size());++m) { require(mult[m]>=0 && mult[m]%2==0,"SL2 multiplicities must be nonnegative even integers"); n+=mult[m]/2*(m+1); if(m%4==1) parity+=mult[m]/2; }
    require(n>0 && parity%2==0,"SL2 Spin-lift parity criterion failed");
    MatType x=MatType::Zero(n,n),z=x; int pos=0;
    for(int m=0;m<int(mult.size());++m) for(int copy=0;copy<mult[m]/2;++copy) {
        x.block(pos,pos,m+1,m+1)=spinX(m); z.block(pos,pos,m+1,m+1)=spinZ(m); pos+=m+1;
    }
    Hamiltonian h; h.id="sl2_orthogonal"; h.modes=n; h.q0=numberMajorana(.2*x+.13*z);
    h.channels={{"common_Sx",numberMajorana(x),.12},{"common_Sz",numberMajorana(z),.12}};
    h.description="Explicit number-conserving representative of the orthogonal SL2 family, sharing fields across multiplicities";
    h.validate(); return h;
}
Hamiltonian transformedHopping(int n) {
    auto h=openHopping(n); h.id="local_unitary_hopping";
    for(int i=0;i+1<n;i+=2) h.controlledPhases.emplace_back(i,i+1,.37);
    h.description="Correlated hopping H'=U H U^dagger with fixed depth-one controlled-density phases; sample in the inverse-transformed basis";
    return h;
}
std::vector<std::string> modelNames() {
    return {"open_hopping","odd_pair_ring","odd_pair_ring_three","four_mode","adjoint_d3","adjoint_d5","adjoint_d7",
      "spin_cell_j1","spin_cell_j2","spin_tensor_3x3","spin_tensor_5x5","spin_lattice_j1","spin_lattice_j2",
      "spin_lattice_mixed","five_orbital","pairing_lattice","pairing_lattice_2d","sl_exterior_4_2","sl_exterior_6_3",
      "sl_exterior_8_4","sp_exterior_2_2","sp_exterior_4_3","sp_exterior_5_5","sl2_u2","sl2_u3","sl2_mixed","local_unitary_hopping"};
}
Hamiltonian makeModel(const std::string& name) {
    Hamiltonian h;
    if(name=="open_hopping") h=openHopping();
    else if(name=="odd_pair_ring") h=oddPairRing();
    else if(name=="odd_pair_ring_three") h=oddPairRing(5,3);
    else if(name=="four_mode") h=fourMode();
    else if(name=="adjoint_d3") h=adjointOdd(3);
    else if(name=="adjoint_d5") h=adjointOdd(5);
    else if(name=="adjoint_d7") h=adjointOdd(7);
    else if(name=="spin_cell_j1") h=spinCell(2);
    else if(name=="spin_cell_j2") h=spinCell(4);
    else if(name=="spin_tensor_3x3") h=spinTensor();
    else if(name=="spin_tensor_5x5") h=spinTensor(2,2);
    else if(name=="spin_lattice_j1") h=spinLattice();
    else if(name=="spin_lattice_j2") h=spinLattice(2);
    else if(name=="spin_lattice_mixed") h=spinLattice(1,2,2,true);
    else if(name=="five_orbital") h=fiveOrbital();
    else if(name=="pairing_lattice") h=pairingLattice();
    else if(name=="pairing_lattice_2d") h=pairingLattice(2,2);
    else if(name=="sl_exterior_4_2") h=exteriorSL();
    else if(name=="sl_exterior_6_3") h=exteriorSL(6,3);
    else if(name=="sl_exterior_8_4") h=exteriorSL(8,4);
    else if(name=="sp_exterior_2_2") h=exteriorSp(2,2);
    else if(name=="sp_exterior_4_3") h=exteriorSp(4,3);
    else if(name=="sp_exterior_5_5") h=exteriorSp(5,5);
    else if(name=="sl2_u2") h=sl2Orthogonal({0,0,2});
    else if(name=="sl2_u3") h=sl2Orthogonal({0,0,0,2});
    else if(name=="sl2_mixed") h=sl2Orthogonal({0,0,2,2});
    else if(name=="local_unitary_hopping") h=transformedHopping();
    else throw std::runtime_error("Unknown model: "+name);
    h.id=name; return h;
}
MatType StaticModel::generator(const std::vector<int>& s) const {
    int n=lx*ly; require(int(s.size())==n,"Static background size"); MatType h=MatType::Zero(n,n);
    for(int i=0;i<n;++i) { require(s[i]==1 || s[i]==-1,"Background must be +/-1"); h(i,i)=.5*U*s[i]; }
    for(auto e:grid(lx,ly)) {
        DataType value=-hopping*std::exp(1.i*flux*double(e.first%ly)); h(e.first,e.second)=value; h(e.second,e.first)=std::conj(value);
    }
    MatType a=-numberMajorana(h);
    if(pairing!=0.) for(auto e:grid(lx,ly)) a-=pairMajorana(n,e.first,e.second,pairing*std::exp(.2i));
    return a;
}
double StaticModel::logScalar(const std::vector<int>& s,double beta) const { return .5*beta*muF*std::accumulate(s.begin(),s.end(),0); }
std::vector<std::string> generatorNames() {
    return {"creation_only_chain","creation_only_spin","directed_two_mode","directed_dag","directed_scc"};
}
Hamiltonian makeGenerator(const std::string& name) {
    Hamiltonian h;
    if(name=="creation_only_chain" || name=="creation_only_spin") {
        h=name=="creation_only_chain" ? openHopping(4):spinCell(2);
        h.conserving=false;h.hermitian=false;
        h.q0+=pairMajorana(h.modes,0,h.modes-1,DataType(.3,.2),false);
    } else {
        require(name=="directed_two_mode" || name=="directed_dag" || name=="directed_scc","Unknown non-Hermitian generator");
        h.modes=name=="directed_two_mode" ? 2:4;h.hermitian=false;
        MatType q=MatType::Zero(h.modes,h.modes);
        for(int i=0;i<h.modes;++i)q(i,i)=.13-.09*i;
        if(name=="directed_two_mode")q(0,1)=.3i;
        if(name=="directed_dag") {q(0,1)=.3i;q(0,2)=.4;q(1,3)=-.2+.3i;q(2,3)=-.3i;}
        if(name=="directed_scc") {q(0,1)=q(1,0)=.3;q(2,3)=q(3,2)=.2;q(1,2)=.4i;q(0,3)=-.2+.1i;}
        h.q0=numberMajorana(q);
        for(int i=0;i<h.modes;++i)h.channels.push_back({"real_density_"+std::to_string(i),density(h.modes,i),.2});
    }
    h.id=name;h.description="Non-Hermitian triangular extension: evaluate signed traces; no Hermitian HS walker claim";
    h.validate();return h;
}
}
