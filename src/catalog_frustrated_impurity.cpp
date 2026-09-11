#include "catalog/frustrated_impurity.h"
#include <algorithm>
#include <array>
#include <cmath>

// Hann, Huffman, Chandrasekharan, Ann.Phys.376 (2017) 63-75:
// https://arxiv.org/html/1608.05144v2, Eqs (1), (23)-(25), Table 4.
// The all-length proof is the paired-loop crossing/turning parity argument
// in Sections IV-V, not a numerical balance test of a small configuration graph.
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
        require(es.info()==Eigen::Success && std::abs(es.eigenvalues()(7)-1.)<1e-12,"Impurity basis projector failed");
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
    for(int v=0;v<4;++v) {require(std::abs(phases[v])>.9,"Disconnected impurity basis phases");b.vectors.col(v)*=phases[v];}
    for(int spin=0;spin<2;++spin) for(int a=0;a<3;++a)
        require((b.vectors.adjoint()*s[spin][a]*b.vectors-wanted[spin][a]).norm()<1e-11,"Impurity cell spin matrix mismatch");
    for(int v=0;v<4;++v) {
        b.covariance[v]=MatType::Zero(6,6);
        for(int i=0;i<6;++i) for(int j=i+1;j<6;++j) {
            b.covariance[v](i,j)=(b.vectors.col(v).adjoint()*g[i]*g[j]*b.vectors.col(v))(0,0);
            b.covariance[v](j,i)=-b.covariance[v](i,j);
        }
        require((b.covariance[v]*b.covariance[v]-MatType::Identity(6,6)).norm()<1e-11,"Impurity local state is not pure Gaussian");
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
        require(tr.p>=0,"No Gaussian impurity transition found");
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
int localIndex(const ImpurityTerm& t,const std::vector<int>& state) {
    int out=0;for(int i:t.cells) out=4*out+state[i];return out;
}
double energy(const std::vector<ImpurityTerm>& ts,const std::vector<int>& s) {
    double out=0.;for(const auto& t:ts) {int i=localIndex(t,s);out+=t.matrix(i,i).real();}return out;
}
std::vector<MatType> fullFockTerms(const FrustratedImpurity& m) {
    require(m.length==2,"Impurity full Fock oracle supports length=2");
    auto s=majoranaSpins(majoranas(m.modes()));
    return {(m.chain+m.armImbalance)*dot(s[0],s[2])+m.chain*dot(s[1],s[3]),m.chain*dot(s[4],s[5]),
            m.impurity*(dot(s[0],s[1])+dot(s[0],s[4])+dot(s[1],s[4]))};
}
cVecType fullState(const FrustratedImpurity& m,const std::vector<int>& state,const std::vector<int>& parities,bool dimer) {
    MatType v=MatType::Ones(1,1);
    // Fock mode cell 0 is the least significant tensor factor.
    for(int cell=m.cells()-1;cell>=0;--cell) v=tensor(v,basis(dimer && cell<m.length,parities[cell]).vectors.col(state[cell]));
    return v.col(0);
}
}
void FrustratedImpurity::validate() const {
    require(length>=2 && length%2==0 && length<=100,"Impurity requires positive even chain length 2..100");
    require(std::isfinite(chain) && chain>=0. && std::isfinite(impurity) && impurity>=0.
        && std::isfinite(armImbalance) && chain+armImbalance>=0.,"Impurity requires finite nonnegative exchange couplings");
}
bool FrustratedImpurity::signFreeParameters() const {validate();return armImbalance==0.;}
std::vector<ImpurityTerm> FrustratedImpurity::terms(bool dimer) const {
    validate();auto l=cellSpins(dimer),r=cellSpins(false);MatType id=MatType::Identity(4,4);std::vector<ImpurityTerm> ts;
    for(int i=0;i<length-1;++i) {
        MatType first=MatType::Zero(16,16),second=first;
        for(int a=0;a<3;++a) {first+=tensor(l[0][a],l[0][a]);second+=tensor(l[1][a],l[1][a]);}
        MatType h=chain*exactZeros(first+second)+armImbalance*exactZeros(first);
        ts.push_back({{i,i+1},h});
    }
    for(int i=0;i<length-1;++i) {
        if(i%2==0) ts.push_back({{length+i/2},chain*exactZeros(dot(r[0],r[1]))});
        else {
            MatType h=MatType::Zero(16,16);for(int a=0;a<3;++a) h+=tensor(r[1][a],r[0][a]);
            ts.push_back({{length+i/2,length+i/2+1},chain*exactZeros(h)});
        }
    }
    MatType h=tensor(dot(l[0],l[1]),id);
    for(int a=0;a<3;++a) h+=tensor(l[0][a]+l[1][a],r[0][a]);
    ts.push_back({{0,length},impurity*exactZeros(h)});return ts;
}
bool FrustratedImpurity::validConfiguration(const ImpurityConfiguration& c,bool dimer) const {
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
double FrustratedImpurity::diagonalEnergy(const std::vector<int>& state,bool dimer) const {
    validate();require(validState(state,cells()),"Invalid impurity cell state");return energy(terms(dimer),state);
}
Gaussian FrustratedImpurity::stateProjector(const std::vector<int>& state,const std::vector<int>& parities,bool dimer) const {
    validate();require(validState(state,cells()) && int(parities.size())==cells(),"Invalid impurity projector state");
    Gaussian out;out.trace=1.;out.g=MatType::Zero(2*modes(),2*modes());
    for(int i=0;i<cells();++i) {
        require(parities[i]==1 || parities[i]==-1,"Invalid impurity cell parity");
        out.g.block(6*i,6*i,6,6)=basis(dimer && i<length,parities[i]).covariance[state[i]];
    }
    return out;
}
DataType FrustratedImpurity::configurationPfaffian(const ImpurityConfiguration& c,bool dimer) const {
    if(!validConfiguration(c,dimer)) return 0.;auto ts=terms(dimer);
    std::vector<Gaussian> fs{stateProjector(c.path.front(),c.parities,dimer)};DataType scalar=1.;double last=0.;
    for(size_t k=0;k<c.events.size();++k) {
        auto e=c.events[k];scalar*=std::exp(-(e.time-last)*energy(ts,c.path[k]));last=e.time;
        const auto& t=ts[e.term];scalar*= -t.matrix(localIndex(t,c.path[k+1]),localIndex(t,c.path[k]));
        for(int i:t.cells) if(c.path[k][i]!=c.path[k+1][i]) {
            const auto& tr=basis(dimer && i<length,c.parities[i]).transitions[c.path[k][i]][c.path[k+1][i]];
            scalar/=tr.amplitude;Gaussian f;f.g=MatType::Zero(2*modes(),2*modes());
            f.g.block(6*i,6*i,6,6)=tr.pulse.g;f.trace=std::ldexp(1.,modes()-3)*tr.pulse.trace;
            fs.push_back(f);fs.push_back(f);
        }
        fs.push_back(stateProjector(c.path[k+1],c.parities,dimer));
    }
    scalar*=std::exp(-(c.beta-last)*energy(ts,c.path.back()));
    require(std::isfinite(std::abs(scalar)) && std::abs(scalar)>0.,"Impurity CT scalar outside floating-point range");
    DataType result=scalar*traceGaussianFactors(fs);require(std::isfinite(std::abs(result)),"Nonfinite impurity Pfaffian");return result;
}
DataType FrustratedImpurity::configurationLocalReference(const ImpurityConfiguration& c,bool dimer) const {
    if(!validConfiguration(c,dimer)) return 0.;auto ts=terms(dimer);DataType scalar=1.;double last=0.;
    for(size_t k=0;k<c.events.size();++k) {
        auto e=c.events[k];scalar*=std::exp(-(e.time-last)*energy(ts,c.path[k]));last=e.time;const auto& t=ts[e.term];
        scalar*= -t.matrix(localIndex(t,c.path[k+1]),localIndex(t,c.path[k]));
    }
    return scalar*std::exp(-(c.beta-last)*energy(ts,c.path.back()));
}
MatType FrustratedImpurity::exactHamiltonian() const {
    validate();auto ts=fullFockTerms(*this);MatType h=MatType::Zero(ts[0].rows(),ts[0].cols());for(const auto& t:ts) h+=t;return h;
}
MatType FrustratedImpurity::exactSpinHamiltonian() const {
    validate();require(length==2,"Impurity spin spectrum oracle supports length=2");auto p=pauli();std::vector<Spin> s;
    for(int i=0;i<6;++i) {
        Spin si;
        for(int a=0;a<3;++a) {MatType v=MatType::Ones(1,1);for(int j=0;j<6;++j) v=tensor(v,i==j ? p[a]:MatType::Identity(2,2));si[a]=v;}
        s.push_back(si);
    }
    return (chain+armImbalance)*dot(s[0],s[2])+chain*(dot(s[1],s[3])+dot(s[4],s[5]))
        +impurity*(dot(s[0],s[1])+dot(s[0],s[4])+dot(s[1],s[4]));
}
DataType FrustratedImpurity::configurationFockReference(const ImpurityConfiguration& c,bool dimer) const {
    if(!validConfiguration(c,dimer)) return 0.;auto ts=fullFockTerms(*this);MatType h=ts[0];for(size_t i=1;i<ts.size();++i) h+=ts[i];
    DataType result=1.;double last=0.;
    for(size_t k=0;k<c.events.size();++k) {
        auto v=fullState(*this,c.path[k],c.parities,dimer),next=fullState(*this,c.path[k+1],c.parities,dimer);auto e=c.events[k];
        double diagonal=(v.adjoint()*h*v)(0,0).real();result*=std::exp(-(e.time-last)*diagonal);last=e.time;
        result*= -(next.adjoint()*ts[e.term]*v)(0,0);
    }
    auto v=fullState(*this,c.path.back(),c.parities,dimer);return result*std::exp(-(c.beta-last)*(v.adjoint()*h*v)(0,0).real());
}
}
