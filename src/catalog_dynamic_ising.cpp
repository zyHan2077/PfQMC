#include "catalog/dynamic_ising.h"
#include <algorithm>
#include <cmath>
#include <queue>
#include <set>

// The CT-INT + spin-worldline proof is from Huffman and Chandrasekharan,
// https://arxiv.org/html/1605.07420v2, Section II. The pair-fermion embedding
// retains the complete Fock space: odd pair parity is a static spin-zero sector.
namespace catalog {
namespace {
constexpr double pi=3.1415926535897932384626433832795;
double spinZ(std::uint64_t v,int i) {
    return (double((v>>(2*i))&1u)+double((v>>(2*i+1))&1u)-1.)/2.;
}
Gaussian embed(const Gaussian& f,int offset,int total) {
    const int n=f.g.rows()/2;
    Gaussian out; out.g=MatType::Zero(2*total,2*total);
    out.g.block(2*offset,2*offset,2*n,2*n)=f.g;
    out.trace=std::ldexp(1.,total-n)*f.trace;
    return out; // Only g and trace are used by the spacetime Pfaffian.
}
MatType annihilation(int n,int i) {
    MatType out=MatType::Zero(1<<n,1<<n);
    for(unsigned v=0;v<(1u<<n);++v) if(v&(1u<<i))
        out(v^(1u<<i),v)=__builtin_popcount(v&((1u<<i)-1))%2 ? -1.:1.;
    return out;
}
struct FockOperators {
    MatType h0;
    std::vector<MatType> dn,sx;
};
// Independently construct the original Fock operators from creation/annihilation
// matrices, rather than reusing the Majorana conversion in the Pfaffian path.
FockOperators fullOperators(const DynamicIsing& model) {
    require(model.sites<=3,"DynamicIsing full Fock oracle supports sites<=3");
    int n=3*model.sites,dim=1<<n;
    MatType id=MatType::Identity(dim,dim);
    std::vector<MatType> c,number,sz;
    for(int i=0;i<n;++i) {c.push_back(annihilation(n,i));number.push_back(c.back().adjoint()*c.back());}
    FockOperators f; f.h0=MatType::Zero(dim,dim);
    for(int i=0;i<model.sites;++i) {
        int a=model.sites+2*i,b=a+1;
        f.dn.push_back(number[i]-.5*id);
        sz.push_back(.5*(number[a]+number[b]-id));
        f.sx.push_back(.5*(c[a].adjoint()*c[b].adjoint()+c[b]*c[a]));
    }
    for(auto [i,j]:model.edges) {
        f.h0-=model.hopping*(c[i].adjoint()*c[j]+c[j].adjoint()*c[i]);
        f.h0-=model.ising*sz[i]*sz[j];
    }
    return f;
}
}
DynamicIsing DynamicIsing::square() {
    DynamicIsing m; m.sites=4; m.edges={{0,1},{1,2},{2,3},{3,0}};
    m.transverse={.6,-.4,.3,.7}; return m;
}
std::vector<int> DynamicIsing::sublattices() const {
    require(sites>0 && sites<=30,"DynamicIsing requires 1<=sites<=30 (occupation mask limit)");
    std::vector<std::vector<int>> neighbors(sites); std::set<std::pair<int,int>> seen;
    for(auto [i,j]:edges) {
        require(i>=0 && j>=0 && i<sites && j<sites && i!=j,"Invalid DynamicIsing edge");
        require(seen.insert(std::minmax(i,j)).second,"Duplicate DynamicIsing edge");
        neighbors[i].push_back(j); neighbors[j].push_back(i);
    }
    std::vector<int> color(sites,0); std::queue<int> todo;
    for(int start=0;start<sites;++start) if(!color[start]) {
        color[start]=1;todo.push(start);
        while(!todo.empty()) {
            int i=todo.front();todo.pop();
            for(int j:neighbors[i]) {
                if(!color[j]) {color[j]=-color[i];todo.push(j);}
                require(color[j]==-color[i],"DynamicIsing graph must be bipartite");
            }
        }
    }
    return color;
}
void DynamicIsing::validate() const {
    (void)sublattices();
    require(std::isfinite(hopping) && std::isfinite(repulsion) && repulsion>=0.
        && std::isfinite(ising),"DynamicIsing requires real finite t,J and V>=0");
    require(int(transverse.size())==sites,"DynamicIsing transverse field size mismatch");
    for(double h:transverse) require(std::isfinite(h),"Nonfinite DynamicIsing h_i");
}
bool DynamicIsing::validEvents(double beta,const std::vector<IsingEvent>& events) const {
    if(!std::isfinite(beta) || beta<=0.) return false;
    double last=0.;
    for(const auto& e:events) {
        if(!std::isfinite(e.time) || e.time<last || e.time>beta || e.index<0) return false;
        if(e.kind==IsingVertex::DensityBond) {if(e.index>=int(edges.size())) return false;}
        else if(e.kind==IsingVertex::SpinFlip) {if(e.index>=sites) return false;}
        else return false;
        last=e.time;
    }
    return true;
}
bool DynamicIsing::validConfiguration(const IsingConfiguration& c) const {
    validate();
    if(!validEvents(c.beta,c.events) || c.spinPath.size()!=c.events.size()+1) return false;
    const auto bound=std::uint64_t(1)<<(2*sites);
    for(auto v:c.spinPath) if(v>=bound) return false;
    if(c.spinPath.front()!=c.spinPath.back()) return false;
    for(size_t k=0;k<c.events.size();++k) {
        auto before=c.spinPath[k],after=c.spinPath[k+1]; const auto& e=c.events[k];
        if(e.kind==IsingVertex::DensityBond) {if(before!=after) return false;}
        else {
            auto pair=(before>>(2*e.index))&3u;
            if(pair!=0 && pair!=3) return false;
            if(after!=(before^(std::uint64_t(3)<<(2*e.index)))) return false;
        }
    }
    return true;
}
MatType DynamicIsing::hoppingMatrix() const {
    MatType h=MatType::Zero(sites,sites);
    for(auto [i,j]:edges) h(i,j)=h(j,i)=hopping;
    return h; // Q=-H_kinetic; trace(h)=0.
}
double DynamicIsing::spinEnergy(std::uint64_t v) const {
    double value=0.;for(auto [i,j]:edges) value-=ising*spinZ(v,i)*spinZ(v,j);return value;
}
double DynamicIsing::spinPathAmplitude(const IsingConfiguration& c) const {
    double out=1.,last=0.;
    for(size_t k=0;k<c.events.size();++k) {
        out*=std::exp(-(c.events[k].time-last)*spinEnergy(c.spinPath[k]));
        if(c.events[k].kind==IsingVertex::SpinFlip) out*=.5;
        last=c.events[k].time;
    }
    return out*std::exp(-(c.beta-last)*spinEnergy(c.spinPath.back()));
}
Gaussian DynamicIsing::spinProjector(std::uint64_t v) const {
    validate();require(v<(std::uint64_t(1)<<(2*sites)),"DynamicIsing projector outside Fock space");
    Gaussian f; f.trace=1.;f.g=MatType::Zero(4*sites,4*sites);
    for(int i=0;i<2*sites;++i) {
        f.g(2*i,2*i+1)=(v&(std::uint64_t(1)<<i)) ? -1.i:1.i;
        f.g(2*i+1,2*i)=-f.g(2*i,2*i+1);
    }
    return f;
}
Gaussian DynamicIsing::spinFlipRotation(int site,bool create) const {
    validate();require(site>=0 && site<sites,"DynamicIsing rotation site outside graph");
    MatType a=pairMajorana(2*sites,2*site,2*site+1,1.,false);
    return gaussian((create ? .5*pi:-.5*pi)*(a-a.adjoint()).eval());
}
DataType DynamicIsing::spinPathPfaffian(const IsingConfiguration& c) const {
    if(!validConfiguration(c)) return 0.;
    std::vector<Gaussian> fs{spinProjector(c.spinPath.front())};
    for(size_t k=0;k<c.events.size();++k) if(c.events[k].kind==IsingVertex::SpinFlip) {
        int i=c.events[k].index;
        fs.push_back(spinProjector(c.spinPath[k]));
        fs.push_back(spinFlipRotation(i,((c.spinPath[k]>>(2*i))&3u)==0));
    }
    return traceGaussianFactors(fs);
}
DataType DynamicIsing::configurationPfaffian(const IsingConfiguration& c) const {
    if(!validConfiguration(c)) return 0.;
    const int total=3*sites;auto color=sublattices();MatType hop=hoppingMatrix();
    std::vector<Gaussian> fs{embed(spinProjector(c.spinPath.front()),sites,total)};
    DataType scalar=1.;double last=0.;
    for(size_t k=0;k<=c.events.size();++k) {
        double time=k<c.events.size() ? c.events[k].time:c.beta,dt=time-last;last=time;
        scalar*=std::exp(-dt*spinEnergy(c.spinPath[k]));
        if(dt!=0.) fs.push_back(embed(gaussian(numberMajorana(dt*hop)),0,total));
        if(k==c.events.size()) break;
        const auto& e=c.events[k];scalar=-scalar;
        std::vector<int> insertions;
        if(e.kind==IsingVertex::DensityBond) {
            if(repulsion==0.) return 0.;
            auto [i,j]=edges[e.index];insertions={i,j};scalar*=repulsion;
        } else {
            if(transverse[e.index]==0.) return 0.;
            insertions={e.index};scalar*=.5*color[e.index]*transverse[e.index];
        }
        require(std::abs(scalar)>0. && std::isfinite(std::abs(scalar)),
                "DynamicIsing CT-INT scalar outside floating-point range");
        for(int i:insertions) {
            MatType h=MatType::Zero(sites,sites);h(i,i)=.5i*pi;
            auto f=embed(gaussian(numberMajorana(h)),0,total);
            fs.push_back(f);fs.push_back(f);scalar*= -.5i;
        }
        if(e.kind==IsingVertex::SpinFlip) {
            fs.push_back(embed(spinProjector(c.spinPath[k]),sites,total));
            bool create=((c.spinPath[k]>>(2*e.index))&3u)==0;
            fs.push_back(embed(spinFlipRotation(e.index,create),sites,total));
        }
    }
    DataType out=scalar*traceGaussianFactors(fs);
    require(std::isfinite(std::abs(out)),"Nonfinite DynamicIsing Pfaffian weight");return out;
}
DataType DynamicIsing::fermionWordTrace(double beta,const std::vector<IsingEvent>& events) const {
    require(sites<=5,"DynamicIsing factorized Fock oracle supports sites<=5");
    int dim=1<<sites;MatType id=MatType::Identity(dim,dim),q=MatType::Zero(dim,dim);
    std::vector<MatType> c,dn;auto color=sublattices();
    for(int i=0;i<sites;++i) {c.push_back(annihilation(sites,i));dn.push_back(c.back().adjoint()*c.back()-.5*id);}
    for(auto [i,j]:edges) q+=hopping*(c[i].adjoint()*c[j]+c[j].adjoint()*c[i]);
    MatType p=id;double last=0.;
    for(const auto& e:events) {
        p=(matrixExp((e.time-last)*q)*p).eval();last=e.time;
        if(e.kind==IsingVertex::DensityBond) {
            auto [i,j]=edges[e.index];p=(-repulsion*dn[i]*dn[j]*p).eval();
        } else p=(-color[e.index]*transverse[e.index]*dn[e.index]*p).eval();
    }
    return (matrixExp((beta-last)*q)*p).trace();
}
DataType DynamicIsing::configurationReference(const IsingConfiguration& c) const {
    if(!validConfiguration(c)) return 0.;
    return fermionWordTrace(c.beta,c.events)*spinPathAmplitude(c);
}
DataType DynamicIsing::factorizedWordTrace(double beta,const std::vector<IsingEvent>& events) const {
    validate();if(!validEvents(beta,events)) return 0.;
    require(sites<=5,"DynamicIsing spin path sum supports sites<=5");
    double spin=0.;
    for(std::uint64_t v=0;v<(std::uint64_t(1)<<(2*sites));++v) {
        IsingConfiguration c{beta,events,{v}};
        for(const auto& e:events) {
            auto next=c.spinPath.back();
            if(e.kind==IsingVertex::SpinFlip) next^=std::uint64_t(3)<<(2*e.index);
            c.spinPath.push_back(next);
        }
        if(validConfiguration(c)) spin+=spinPathAmplitude(c);
    }
    return fermionWordTrace(beta,events)*spin;
}
MatType DynamicIsing::exactHamiltonian(bool staggered) const {
    validate();auto f=fullOperators(*this);auto color=sublattices();MatType h=f.h0;
    for(auto [i,j]:edges) h+=repulsion*f.dn[i]*f.dn[j];
    for(int i=0;i<sites;++i) h+=transverse[i]*(staggered ? color[i]:1)*f.dn[i]*f.sx[i];
    return h;
}
DataType DynamicIsing::exactWordTrace(double beta,const std::vector<IsingEvent>& events) const {
    validate();if(!validEvents(beta,events)) return 0.;
    auto f=fullOperators(*this);auto color=sublattices();
    MatType p=MatType::Identity(f.h0.rows(),f.h0.cols());double last=0.;
    for(const auto& e:events) {
        p=(matrixExp(-(e.time-last)*f.h0)*p).eval();last=e.time;
        if(e.kind==IsingVertex::DensityBond) {
            auto [i,j]=edges[e.index];p=(-repulsion*f.dn[i]*f.dn[j]*p).eval();
        } else p=(-transverse[e.index]*color[e.index]*f.dn[e.index]*f.sx[e.index]*p).eval();
    }
    return (matrixExp(-(beta-last)*f.h0)*p).trace();
}
}
