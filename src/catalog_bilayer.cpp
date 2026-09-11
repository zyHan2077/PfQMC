#include "catalog/bilayer_models.h"
#include <algorithm>
#include <array>
#include <cmath>

// Alet, Damle, Pujari, arXiv:1511.01586v2, Eqs.(1),(3), closed-path parity proof.
// The three parameter faces form ONE construction. Finite-graph balance in
// the validator is evidence, not the all-size proof.
// Gaussian dimer basis construction follows catalog_frustrated_impurity.cpp.
namespace catalog {
namespace {
using Spin=std::array<MatType,3>;
MatType tensor(const MatType& a,const MatType& b) {
    MatType out(a.rows()*b.rows(),a.cols()*b.cols());
    for(int i=0;i<a.rows();++i) for(int j=0;j<a.cols();++j)
        out.block(i*b.rows(),j*b.cols(),b.rows(),b.cols())=a(i,j)*b;
    return out;
}
MatType annihilation(int n,int i) {
    MatType out=MatType::Zero(1<<n,1<<n);
    for(unsigned v=0;v<(1u<<n);++v) if(v&(1u<<i))
        out(v^(1u<<i),v)=__builtin_popcount(v&((1u<<i)-1))%2 ? -1.:1.;
    return out;
}
std::vector<MatType> majoranas(int n) {
    std::vector<MatType> out;
    for(int i=0;i<n;++i) {auto a=annihilation(n,i);out.push_back(a+a.adjoint());out.push_back(1.i*(a.adjoint()-a));}
    return out;
}
std::vector<Spin> majoranaSpins(const std::vector<MatType>& g) {
    std::vector<Spin> out;
    for(size_t i=0;i<g.size()/3;++i) out.push_back({-.5i*g[3*i+1]*g[3*i+2],-.5i*g[3*i+2]*g[3*i],-.5i*g[3*i]*g[3*i+1]});
    return out;
}
Spin pauli() {
    Spin s{MatType::Zero(2,2),MatType::Zero(2,2),MatType::Zero(2,2)};
    s[0](0,1)=s[0](1,0)=.5;s[1](0,1)=-.5i;s[1](1,0)=.5i;s[2](0,0)=.5;s[2](1,1)=-.5;return s;
}
MatType dimerRotation() {
    MatType u=MatType::Zero(4,4);double a=1./std::sqrt(2.);
    u(1,0)=a;u(2,0)=-a;u(1,1)=u(2,1)=a;u(0,2)=1.;u(3,3)=-1.;return u;
}
// Remove representation roundoff BEFORE multiplying physical couplings. This
// avoids interpreting forbidden vertices as tiny nonzero ones, without cutting
// off a genuinely small exchange coupling in validConfiguration().
MatType exactZeros(MatType a) {
    for(int i=0;i<a.rows();++i) for(int j=0;j<a.cols();++j)
        if(std::abs(a(i,j))<1e-14) a(i,j)=0.;
    return a;
}
std::array<Spin,2> cellSpins(bool dimer) {
    auto p=pauli();MatType id=MatType::Identity(2,2),u=dimer ? dimerRotation():MatType::Identity(4,4);
    std::array<Spin,2> s;
    for(int a=0;a<3;++a) {s[0][a]=exactZeros(u.adjoint()*tensor(p[a],id)*u);s[1][a]=exactZeros(u.adjoint()*tensor(id,p[a])*u);}
    return s;
}
MatType dot(const Spin& a,const Spin& b) {
    MatType out=MatType::Zero(a[0].rows(),a[0].cols());for(int k=0;k<3;++k) out+=a[k]*b[k];return out;
}
struct Transition {int p=-1,q=-1;DataType amplitude=0.;Gaussian pulse;};
struct Basis {
    MatType vectors;
    std::array<MatType,4> covariance;
    std::array<std::array<Transition,4>,4> transitions;
};
Basis makeBasis(bool dimer,int parity) {
    const auto g=majoranas(3);const auto s=majoranaSpins(g);
    MatType id=MatType::Identity(8,8),pz=MatType::Zero(8,8);
    for(int v=0;v<8;++v) pz(v,v)=__builtin_popcount(unsigned(v))%2 ? -1.:1.;
    MatType pp=.5*(id+double(parity)*pz);std::array<MatType,4> projectors;
    if(dimer) {
        projectors[0]=.25*id-dot(s[0],s[1]);MatType z=s[0][2]+s[1][2];
        projectors[1]=(id-projectors[0])*(id-z*z);
        projectors[2]=(.5*id+s[0][2])*(.5*id+s[1][2]);
        projectors[3]=(.5*id-s[0][2])*(.5*id-s[1][2]);
    } else for(int v=0;v<4;++v) projectors[v]=(.5*id+(v&2 ? -1.:1.)*s[0][2])*(.5*id+(v&1 ? -1.:1.)*s[1][2]);
    Basis b;b.vectors=MatType::Zero(8,4);
    for(int v=0;v<4;++v) {
        MatType p=pp*projectors[v];Eigen::SelfAdjointEigenSolver<MatType> es(.5*(p+p.adjoint()));
        require(es.info()==Eigen::Success && std::abs(es.eigenvalues()(7)-1.)<1e-12,"Bilayer basis projector failed");
        b.vectors.col(v)=es.eigenvectors().col(7);
    }
    auto wanted=cellSpins(dimer);std::array<DataType,4> phases{1.,0.,0.,0.};std::vector<int> todo{0};
    for(size_t k=0;k<todo.size();++k) {
        int i=todo[k];
        for(int spin=0;spin<2;++spin) for(int a=0;a<3;++a) {
            MatType actual=b.vectors.adjoint()*s[spin][a]*b.vectors;
            for(int j=0;j<4;++j) if(std::abs(phases[j])==0. && std::abs(wanted[spin][a](j,i))>1e-10) {
                phases[j]=phases[i]*actual(j,i)/wanted[spin][a](j,i);phases[j]/=std::abs(phases[j]);todo.push_back(j);
            }
        }
    }
    for(int v=0;v<4;++v) {require(std::abs(phases[v])>.9,"Disconnected bilayer basis phases");b.vectors.col(v)*=phases[v];}
    for(int spin=0;spin<2;++spin) for(int a=0;a<3;++a)
        require((b.vectors.adjoint()*s[spin][a]*b.vectors-wanted[spin][a]).norm()<1e-11,"Bilayer cell spin matrix mismatch");
    for(int v=0;v<4;++v) {
        b.covariance[v]=MatType::Zero(6,6);
        for(int i=0;i<6;++i) for(int j=i+1;j<6;++j) {
            b.covariance[v](i,j)=(b.vectors.col(v).adjoint()*g[i]*g[j]*b.vectors.col(v))(0,0);
            b.covariance[v](j,i)=-b.covariance[v](i,j);
        }
        require((b.covariance[v]*b.covariance[v]-MatType::Identity(6,6)).norm()<1e-11,"Bilayer local state is not pure Gaussian");
    }
    for(int before=0;before<4;++before) for(int after=0;after<4;++after) if(before!=after) {
        auto& tr=b.transitions[before][after];
        for(int p=0;p<6 && tr.p<0;++p) for(int q=p+1;q<6 && tr.p<0;++q) {
            DataType u=(b.vectors.col(after).adjoint()*g[p]*g[q]*b.vectors.col(before))(0,0);
            if(std::abs(u)>.1) {
                tr.p=p;tr.q=q;tr.amplitude=u;MatType a=MatType::Zero(6,6);
                // Two pulses give exp(pi/2 gamma_p gamma_q)=gamma_p gamma_q.
                a(p,q)=.5*std::acos(-1.);a(q,p)=-a(p,q);tr.pulse=gaussian(a);
            }
        }
        require(tr.p>=0,"No Gaussian bilayer transition found");
    }
    return b;
}
const Basis& basis(bool dimer,int parity) {
    static const std::array<Basis,4> all{makeBasis(false,1),makeBasis(false,-1),makeBasis(true,1),makeBasis(true,-1)};
    return all[(dimer ? 2:0)+(parity<0 ? 1:0)];
}
bool validState(const std::vector<int>& state,int cells) {
    if(int(state.size())!=cells) return false;for(int s:state) if(s<0 || s>=4) return false;return true;
}
int localIndex(const BilayerTerm& t,const std::vector<int>& state) {
    int out=0;for(int i:t.cells) out=4*out+state[i];return out;
}
double energy(const std::vector<BilayerTerm>& ts,const std::vector<int>& s) {
    double out=0.;for(const auto& t:ts) {int i=localIndex(t,s);out+=t.matrix(i,i).real();}return out;
}

struct SpinProduct {int first,second,axis;double coupling;};
std::vector<std::vector<SpinProduct>> physicalTerms(const AnisotropicBilayer& m) {
    std::vector<std::vector<SpinProduct>> out;
    for(auto edge:m.edges) {
        std::vector<SpinProduct> term;
        for(int a=0;a<3;++a) for(int l=0;l<2;++l) for(int k=0;k<2;++k)
            term.push_back({2*edge.first+l,2*edge.second+k,a,
                (a==2?m.Jz:m.Jperp)+(l==k?1.:-1.)*(a==2?m.Kz:m.Kperp)});
        out.push_back(term);
    }
    for(int i=0;i<m.sites;++i) {
        std::vector<SpinProduct> term;
        for(int a=0;a<3;++a) term.push_back({2*i,2*i+1,a,a==2?m.Dz:m.Dperp});out.push_back(term);
    }
    return out;
}
cVecType physicalAction(const std::vector<SpinProduct>& terms,const cVecType& v) {
    cVecType out=cVecType::Zero(v.size());
    // Direct full occupation action of FOUR Majoranas, independent of spin
    // tensor matrices, dimer matrices, and Gaussian covariance contractions.
    for(auto t:terms) for(unsigned start=0;start<unsigned(v.size());++start) if(v(start)!=DataType(0.)) {
        unsigned state=start;DataType amp=-.25*t.coupling*v(start);
        int ia=3*t.first+(t.axis+1)%3,ib=3*t.first+(t.axis+2)%3;
        int ja=3*t.second+(t.axis+1)%3,jb=3*t.second+(t.axis+2)%3;
        for(int g:{jb,ja,ib,ia}) {
            int mode=g/2,occ=(state>>mode)&1;
            amp *= __builtin_popcount(state&((1u<<mode)-1))%2 ? -1.:1.;
            if(g%2)amp*=1.i*double(1-2*occ);state^=1u<<mode;
        }
        out(state)+=amp;
    }
    return out;
}
cVecType fullState(const AnisotropicBilayer& m,const std::vector<int>& state,const std::vector<int>& parities,bool dimer) {
    MatType v=MatType::Ones(1,1);
    for(int i=m.cells()-1;i>=0;--i)v=tensor(v,basis(dimer,parities[i]).vectors.col(state[i]));
    return v.col(0);
}
}
void AnisotropicBilayer::validate() const {
    require(sites>=1&&sites<=100,"Bilayer needs 1..100 sites");
    for(double c:{Jz,Jperp,Kz,Kperp,Dz,Dperp})require(std::isfinite(c),"Nonfinite bilayer coupling");
    std::vector<std::pair<int,int>> sorted;
    for(auto e:edges) {
        require(e.first>=0&&e.second>=0&&e.first<sites&&e.second<sites&&e.first!=e.second,"Invalid bilayer edge");
        if(e.first>e.second)std::swap(e.first,e.second);
        require(std::find(sorted.begin(),sorted.end(),e)==sorted.end(),"Repeated bilayer edge");sorted.push_back(e);
    }
}
bool AnisotropicBilayer::bipartite() const {
    validate();std::vector<int> color(sites,-1);
    for(int r=0;r<sites;++r)if(color[r]<0) {
        color[r]=0;std::vector<int> todo{r};
        for(size_t k=0;k<todo.size();++k)for(auto e:edges) {
            int y=e.first==todo[k]?e.second:e.second==todo[k]?e.first:-1;if(y<0)continue;
            if(color[y]<0){color[y]=1-color[todo[k]];todo.push_back(y);}
            else if(color[y]==color[todo[k]])return false;
        }
    }
    return true;
}
bool AnisotropicBilayer::signFreeParameters() const {
    return bipartite()&&(Kz==0.||Kperp==0.||Jperp==0.);
}
std::vector<BilayerTerm> AnisotropicBilayer::terms(bool dimer) const {
    validate();auto s=cellSpins(dimer);Spin l,n;
    for(int a=0;a<3;++a){l[a]=exactZeros(s[0][a]+s[1][a]);n[a]=exactZeros(s[0][a]-s[1][a]);}
    MatType h=Jz*exactZeros(tensor(l[2],l[2]))+Kz*exactZeros(tensor(n[2],n[2]));
    h+=Jperp*exactZeros(tensor(l[0],l[0])+tensor(l[1],l[1]));
    h+=Kperp*exactZeros(tensor(n[0],n[0])+tensor(n[1],n[1]));
    std::vector<BilayerTerm> ts;for(auto e:edges)ts.push_back({{e.first,e.second},h});
    MatType onsite=Dz*exactZeros(s[0][2]*s[1][2])+Dperp*exactZeros(s[0][0]*s[1][0]+s[0][1]*s[1][1]);
    for(int i=0;i<sites;++i)ts.push_back({{i},onsite});return ts;
}
bool AnisotropicBilayer::validConfiguration(const BilayerConfiguration& c,bool dimer) const {
    validate();if(!std::isfinite(c.beta) || c.beta<=0. || int(c.parities.size())!=cells() || c.path.size()!=c.events.size()+1) return false;
    for(int p:c.parities) if(p!=1 && p!=-1) return false;
    for(const auto& s:c.path) if(!validState(s,cells())) return false;
    if(c.path.front()!=c.path.back()) return false;
    auto ts=terms(dimer);double last=0.;
    for(size_t k=0;k<c.events.size();++k) {
        auto e=c.events[k];if(!std::isfinite(e.time) || e.time<last || e.time>c.beta || e.term<0 || e.term>=int(ts.size())) return false;
        last=e.time;const auto& t=ts[e.term];const auto& before=c.path[k];const auto& after=c.path[k+1];
        if(before==after) return false;
        for(int i=0;i<cells();++i) if(before[i]!=after[i] && std::find(t.cells.begin(),t.cells.end(),i)==t.cells.end()) return false;
        if(t.matrix(localIndex(t,after),localIndex(t,before))==DataType(0.)) return false;
    }
    return true;
}
double AnisotropicBilayer::diagonalEnergy(const std::vector<int>& state,bool dimer) const {
    validate();require(validState(state,cells()),"Invalid bilayer cell state");return energy(terms(dimer),state);
}
Gaussian AnisotropicBilayer::stateProjector(const std::vector<int>& state,const std::vector<int>& parities,bool dimer) const {
    validate();require(validState(state,cells()) && int(parities.size())==cells(),"Invalid bilayer projector state");
    Gaussian out;out.trace=1.;out.g=MatType::Zero(2*modes(),2*modes());
    for(int i=0;i<cells();++i) {
        require(parities[i]==1 || parities[i]==-1,"Invalid bilayer cell parity");
        out.g.block(6*i,6*i,6,6)=basis(dimer,parities[i]).covariance[state[i]];
    }
    return out;
}
DataType AnisotropicBilayer::configurationPfaffian(const BilayerConfiguration& c,bool dimer) const {
    if(!validConfiguration(c,dimer)) return 0.;auto ts=terms(dimer);
    std::vector<Gaussian> fs{stateProjector(c.path.front(),c.parities,dimer)};DataType scalar=1.;double last=0.;
    for(size_t k=0;k<c.events.size();++k) {
        auto e=c.events[k];scalar*=std::exp(-(e.time-last)*energy(ts,c.path[k]));last=e.time;
        const auto& t=ts[e.term];scalar*= -t.matrix(localIndex(t,c.path[k+1]),localIndex(t,c.path[k]));
        for(int i:t.cells) if(c.path[k][i]!=c.path[k+1][i]) {
            const auto& tr=basis(dimer,c.parities[i]).transitions[c.path[k][i]][c.path[k+1][i]];
            scalar/=tr.amplitude;Gaussian f;f.g=MatType::Zero(2*modes(),2*modes());
            f.g.block(6*i,6*i,6,6)=tr.pulse.g;f.trace=std::ldexp(1.,modes()-3)*tr.pulse.trace;
            fs.push_back(f);fs.push_back(f);
        }
        fs.push_back(stateProjector(c.path[k+1],c.parities,dimer));
    }
    scalar*=std::exp(-(c.beta-last)*energy(ts,c.path.back()));
    require(std::isfinite(std::abs(scalar)) && std::abs(scalar)>0.,"Bilayer CT scalar outside floating-point range");
    DataType result=scalar*traceGaussianFactors(fs);require(std::isfinite(std::abs(result)),"Nonfinite bilayer Pfaffian");return result;
}
DataType AnisotropicBilayer::configurationLocalReference(const BilayerConfiguration& c,bool dimer) const {
    if(!validConfiguration(c,dimer)) return 0.;auto ts=terms(dimer);DataType scalar=1.;double last=0.;
    for(size_t k=0;k<c.events.size();++k) {
        auto e=c.events[k];scalar*=std::exp(-(e.time-last)*energy(ts,c.path[k]));last=e.time;const auto& t=ts[e.term];
        scalar*= -t.matrix(localIndex(t,c.path[k+1]),localIndex(t,c.path[k]));
    }
    return scalar*std::exp(-(c.beta-last)*energy(ts,c.path.back()));
}

DataType AnisotropicBilayer::configurationFockReference(const BilayerConfiguration& c,bool dimer) const {
    if(!validConfiguration(c,dimer))return 0.;require(sites<=4,"Bilayer full occupation oracle limited to 12 modes");
    auto ts=physicalTerms(*this);std::vector<SpinProduct> all;
    for(const auto& t:ts)all.insert(all.end(),t.begin(),t.end());
    DataType out=1.;double last=0.;
    for(size_t k=0;k<c.events.size();++k) {
        auto before=fullState(*this,c.path[k],c.parities,dimer),after=fullState(*this,c.path[k+1],c.parities,dimer);auto e=c.events[k];
        double energy=(before.adjoint()*physicalAction(all,before))(0,0).real();
        out*=std::exp(-(e.time-last)*energy);last=e.time;
        out*=-(after.adjoint()*physicalAction(ts[e.term],before))(0,0);
    }
    auto v=fullState(*this,c.path.back(),c.parities,dimer);
    return out*std::exp(-(c.beta-last)*(v.adjoint()*physicalAction(all,v))(0,0).real());
}
MatType AnisotropicBilayer::exactHamiltonian() const {
    validate();require(sites<=2,"Dense physical bilayer H limited to six modes");int dim=1<<modes();
    auto ts=physicalTerms(*this);std::vector<SpinProduct> all;
    for(const auto& t:ts)all.insert(all.end(),t.begin(),t.end());
    MatType h(dim,dim);for(int i=0;i<dim;++i){cVecType v=cVecType::Zero(dim);v(i)=1.;h.col(i)=physicalAction(all,v);}return h;
}
MatType AnisotropicBilayer::exactSpinHamiltonian(bool dimer) const {
    validate();require(sites<=4,"Dense logical bilayer H limited to eight spins");int dim=1<<(2*sites);MatType h=MatType::Zero(dim,dim);
    for(const auto& t:terms(dimer)) for(int start=0;start<dim;++start) {
        std::vector<int> s(sites);int v=start;for(int i=sites-1;i>=0;--i){s[i]=v%4;v/=4;}
        int before=localIndex(t,s);
        for(int after=0;after<t.matrix.rows();++after)if(t.matrix(after,before)!=DataType(0.)) {
            auto next=s;int v=after;
            for(int k=int(t.cells.size())-1;k>=0;--k){next[t.cells[k]]=v%4;v/=4;}
            int target=0;for(int q:next)target=4*target+q;h(target,start)+=t.matrix(after,before);
        }
    }
    return h;
}
AnisotropicBilayer squareBilayer(int face) {
    require(face>=0&&face<=3,"Bilayer face must be 0..3");AnisotropicBilayer m;
    m.sites=4;m.edges={{0,1},{1,2},{2,3},{3,0}};m.Kz=.4;
    if(face==0)m.Kz=0.;else if(face==1)m.Kperp=0.;else if(face==2)m.Jperp=0.;
    return m;
}
}
