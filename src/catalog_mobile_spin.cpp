#include "catalog/mobile_spin.h"
#include <cmath>
#include <queue>

// Huffman/Chandrasekharan PRE94 043311 (2016), Section III:
// https://arxiv.org/html/1605.07420v2 . Independent real Jz is a direct
// extension: only the positive diagonal worldline factor depends on Jz.
namespace catalog {
namespace {
double spinZ(std::uint64_t v,int i) {
    return (double((v>>(2*i))&1u)+double((v>>(2*i+1))&1u)-1.)/2.;
}
Gaussian embed(const Gaussian& f,int offset,int total) {
    const int n=f.g.rows()/2;Gaussian out;
    out.g=MatType::Zero(2*total,2*total);
    out.g.block(2*offset,2*offset,2*n,2*n)=f.g;
    out.trace=std::ldexp(1.,total-n)*f.trace;return out;
}
MatType annihilation(int n,int i) {
    MatType out=MatType::Zero(1<<n,1<<n);
    for(unsigned v=0;v<(1u<<n);++v) if(v&(1u<<i))
        out(v^(1u<<i),v)=__builtin_popcount(v&((1u<<i)-1))%2 ? -1.:1.;
    return out;
}
struct FockOperators {
    MatType h0;
    std::vector<MatType> dn,sx,pairHop;
};
FockOperators fullOperators(const MobileSpin& m) {
    require(m.sites<=3,"MobileSpin full Fock oracle supports sites<=3");
    int n=3*m.sites,dim=1<<n;MatType id=MatType::Identity(dim,dim);
    std::vector<MatType> c,number,sz,creation;
    for(int i=0;i<n;++i) {c.push_back(annihilation(n,i));number.push_back(c.back().adjoint()*c.back());}
    FockOperators f;f.h0=MatType::Zero(dim,dim);
    for(int i=0;i<m.sites;++i) {
        int a=m.sites+2*i,b=a+1;
        f.dn.push_back(number[i]-.5*id);
        sz.push_back(.5*(number[a]+number[b]-id));
        creation.push_back(c[a].adjoint()*c[b].adjoint());
        f.sx.push_back(.5*(creation.back()+creation.back().adjoint()));
    }
    for(auto [i,j]:m.edges) {
        f.h0-=m.hopping*(c[i].adjoint()*c[j]+c[j].adjoint()*c[i]);
        f.h0+=m.longitudinal*sz[i]*sz[j];
        f.pairHop.push_back(.5*m.exchange*(creation[i]*creation[j].adjoint()+creation[j]*creation[i].adjoint()));
    }
    return f;
}
}
MobileSpin MobileSpin::square() {
    MobileSpin m;m.sites=4;m.edges={{0,1},{1,2},{2,3},{3,0}};
    m.transverse={.6,.4,.3,.7};return m;
}
DynamicIsing MobileSpin::projectorModel() const {
    DynamicIsing m;m.sites=sites;m.edges=edges;m.hopping=hopping;m.repulsion=repulsion;
    m.ising=-longitudinal;m.transverse=transverse;return m;
}
std::vector<int> MobileSpin::sublattices() const {return projectorModel().sublattices();}
void MobileSpin::validate() const {
    projectorModel().validate();
    require(std::isfinite(exchange) && exchange>=0.,"MobileSpin requires finite Jxy>=0");
}
bool MobileSpin::signFreeParameters() const {
    validate();if(exchange==0.) return true;
    std::vector<std::vector<int>> neighbors(sites);
    for(auto [i,j]:edges) {neighbors[i].push_back(j);neighbors[j].push_back(i);}
    std::vector<bool> seen(sites,false);std::queue<int> todo;
    for(int start=0;start<sites;++start) if(!seen[start]) {
        seen[start]=true;todo.push(start);int sign=0;
        while(!todo.empty()) {
            int i=todo.front();todo.pop();
            int local=(transverse[i]>0.)-(transverse[i]<0.);
            if(local && sign && local!=sign) return false;
            if(local) sign=local;
            for(int j:neighbors[i]) if(!seen[j]) {seen[j]=true;todo.push(j);}
        }
    }
    return true;
}
bool MobileSpin::validEvents(double beta,const std::vector<MobileEvent>& events) const {
    if(!std::isfinite(beta) || beta<=0.) return false;
    double last=0.;
    for(auto e:events) {
        if(!std::isfinite(e.time) || e.time<last || e.time>beta || e.index<0) return false;
        if(e.kind==MobileVertex::SpinFlip) {if(e.index>=sites) return false;}
        else if(e.kind==MobileVertex::DensityBond || e.kind==MobileVertex::PairHop) {
            if(e.index>=int(edges.size())) return false;
        } else return false;
        last=e.time;
    }
    return true;
}
MobileConfiguration MobileSpin::resolve(double beta,const std::vector<MobileEvent>& events,std::uint64_t initial) const {
    validate();MobileConfiguration c{beta,events,{initial}};
    if(!validEvents(beta,events)) return c;
    for(auto e:events) {
        auto next=c.spinPath.back();
        if(e.kind==MobileVertex::SpinFlip) next^=std::uint64_t(3)<<(2*e.index);
        else if(e.kind==MobileVertex::PairHop) {
            auto [i,j]=edges[e.index];next^=(std::uint64_t(3)<<(2*i))|(std::uint64_t(3)<<(2*j));
        }
        c.spinPath.push_back(next);
    }
    return c;
}
bool MobileSpin::validConfiguration(const MobileConfiguration& c) const {
    validate();
    if(!validEvents(c.beta,c.events) || c.spinPath.size()!=c.events.size()+1) return false;
    for(auto v:c.spinPath) if(v>=(std::uint64_t(1)<<(2*sites))) return false;
    if(c.spinPath.front()!=c.spinPath.back()) return false;
    for(size_t k=0;k<c.events.size();++k) {
        auto before=c.spinPath[k],after=c.spinPath[k+1],expected=before;auto e=c.events[k];
        if(e.kind==MobileVertex::SpinFlip) {
            auto pair=(before>>(2*e.index))&3u;if(pair!=0 && pair!=3) return false;
            expected^=std::uint64_t(3)<<(2*e.index);
        } else if(e.kind==MobileVertex::PairHop) {
            auto [i,j]=edges[e.index];auto a=(before>>(2*i))&3u,b=(before>>(2*j))&3u;
            if(!((a==0 && b==3)||(a==3 && b==0))) return false;
            expected^=(std::uint64_t(3)<<(2*i))|(std::uint64_t(3)<<(2*j));
        }
        if(after!=expected) return false;
    }
    return true;
}
Gaussian MobileSpin::spinProjector(std::uint64_t v) const {return projectorModel().spinProjector(v);}
Gaussian MobileSpin::spinFlipRotation(int site,bool create) const {return projectorModel().spinFlipRotation(site,create);}
double MobileSpin::spinEnergy(std::uint64_t v) const {
    double out=0.;for(auto [i,j]:edges) out+=longitudinal*spinZ(v,i)*spinZ(v,j);return out;
}
double MobileSpin::spinPathAmplitude(const MobileConfiguration& c) const {
    double out=1.,last=0.;
    for(size_t k=0;k<c.events.size();++k) {
        auto e=c.events[k];out*=std::exp(-(e.time-last)*spinEnergy(c.spinPath[k]));last=e.time;
        if(e.kind==MobileVertex::SpinFlip) out*=.5;
        else if(e.kind==MobileVertex::PairHop) out*=.5*exchange;
    }
    return out*std::exp(-(c.beta-last)*spinEnergy(c.spinPath.back()));
}
DataType MobileSpin::spinPathPfaffian(const MobileConfiguration& c) const {
    if(!validConfiguration(c)) return 0.;
    std::vector<Gaussian> fs{spinProjector(c.spinPath.front())};
    for(size_t k=0;k<c.events.size();++k) if(c.events[k].kind!=MobileVertex::DensityBond) {
        fs.push_back(spinProjector(c.spinPath[k]));
        for(int i=0;i<sites;++i) if(((c.spinPath[k]^c.spinPath[k+1])>>(2*i))&3u)
            fs.push_back(spinFlipRotation(i,((c.spinPath[k]>>(2*i))&3u)==0));
    }
    return traceGaussianFactors(fs);
}
DataType MobileSpin::configurationPfaffian(const MobileConfiguration& c) const {
    if(!validConfiguration(c)) return 0.;
    const int total=3*sites;auto color=sublattices();MatType hop=MatType::Zero(sites,sites);
    for(auto [i,j]:edges) hop(i,j)=hop(j,i)=hopping;
    std::vector<Gaussian> fs{embed(spinProjector(c.spinPath.front()),sites,total)};
    DataType scalar=1.;double last=0.,pi=std::acos(-1.);
    for(size_t k=0;k<=c.events.size();++k) {
        double time=k<c.events.size() ? c.events[k].time:c.beta,dt=time-last;last=time;
        scalar*=std::exp(-dt*spinEnergy(c.spinPath[k]));
        if(dt!=0.) fs.push_back(embed(gaussian(numberMajorana(dt*hop)),0,total));
        if(k==c.events.size()) break;
        auto e=c.events[k];std::vector<int> insertions;
        if(e.kind==MobileVertex::DensityBond) {
            if(repulsion==0.) return 0.;
            auto [i,j]=edges[e.index];insertions={i,j};scalar*= -repulsion;
        } else if(e.kind==MobileVertex::SpinFlip) {
            if(transverse[e.index]==0.) return 0.;
            insertions={e.index};scalar*= -.5*color[e.index]*transverse[e.index];
        } else {
            if(exchange==0.) return 0.;scalar*=.5*exchange;
        }
        require(std::abs(scalar)>0. && std::isfinite(std::abs(scalar)),"MobileSpin CT scalar outside floating-point range");
        for(int i:insertions) {
            MatType h=MatType::Zero(sites,sites);h(i,i)=.5i*pi;
            auto f=embed(gaussian(numberMajorana(h)),0,total);
            fs.push_back(f);fs.push_back(f);scalar*= -.5i;
        }
        if(e.kind!=MobileVertex::DensityBond) {
            fs.push_back(embed(spinProjector(c.spinPath[k]),sites,total));
            for(int i=0;i<sites;++i) if(((c.spinPath[k]^c.spinPath[k+1])>>(2*i))&3u)
                fs.push_back(embed(spinFlipRotation(i,((c.spinPath[k]>>(2*i))&3u)==0),sites,total));
        }
    }
    auto out=scalar*traceGaussianFactors(fs);
    require(std::isfinite(std::abs(out)),"Nonfinite MobileSpin Pfaffian weight");return out;
}
DataType MobileSpin::fermionWordTrace(double beta,const std::vector<MobileEvent>& events) const {
    require(sites<=5,"MobileSpin factorized Fock oracle supports sites<=5");
    int dim=1<<sites;MatType id=MatType::Identity(dim,dim),q=MatType::Zero(dim,dim);
    std::vector<MatType> c,dn;auto color=sublattices();
    for(int i=0;i<sites;++i) {c.push_back(annihilation(sites,i));dn.push_back(c.back().adjoint()*c.back()-.5*id);}
    for(auto [i,j]:edges) q+=hopping*(c[i].adjoint()*c[j]+c[j].adjoint()*c[i]);
    MatType p=id;double last=0.;
    for(auto e:events) {
        p=(matrixExp((e.time-last)*q)*p).eval();last=e.time;
        if(e.kind==MobileVertex::DensityBond) {
            auto [i,j]=edges[e.index];p=(-repulsion*dn[i]*dn[j]*p).eval();
        } else if(e.kind==MobileVertex::SpinFlip) p=(-color[e.index]*transverse[e.index]*dn[e.index]*p).eval();
    }
    return (matrixExp((beta-last)*q)*p).trace();
}
DataType MobileSpin::configurationReference(const MobileConfiguration& c) const {
    if(!validConfiguration(c)) return 0.;return fermionWordTrace(c.beta,c.events)*spinPathAmplitude(c);
}
DataType MobileSpin::factorizedWordTrace(double beta,const std::vector<MobileEvent>& events) const {
    validate();if(!validEvents(beta,events)) return 0.;
    require(sites<=5,"MobileSpin summed path oracle supports sites<=5");
    double spin=0.;
    for(std::uint64_t v=0;v<(std::uint64_t(1)<<(2*sites));++v) {
        auto c=resolve(beta,events,v);if(validConfiguration(c)) spin+=spinPathAmplitude(c);
    }
    return fermionWordTrace(beta,events)*spin;
}
MatType MobileSpin::exactHamiltonian(bool staggered) const {
    validate();auto f=fullOperators(*this);auto color=sublattices();MatType h=f.h0;
    for(auto [i,j]:edges) h+=repulsion*f.dn[i]*f.dn[j];
    for(const auto& p:f.pairHop) h+=(staggered ? -1.:1.)*p;
    for(int i=0;i<sites;++i) h+=transverse[i]*(staggered ? color[i]:1)*f.dn[i]*f.sx[i];
    return h;
}
DataType MobileSpin::exactWordTrace(double beta,const std::vector<MobileEvent>& events) const {
    validate();if(!validEvents(beta,events)) return 0.;
    auto f=fullOperators(*this);MatType p=MatType::Identity(f.h0.rows(),f.h0.cols());double last=0.;
    // Original basis: no staggered h; pair-hop vertices have the original minus
    // CT-INT coefficient. This independently checks both Marshall signs.
    for(auto e:events) {
        p=(matrixExp(-(e.time-last)*f.h0)*p).eval();last=e.time;
        if(e.kind==MobileVertex::DensityBond) {
            auto [i,j]=edges[e.index];p=(-repulsion*f.dn[i]*f.dn[j]*p).eval();
        } else if(e.kind==MobileVertex::SpinFlip) p=(-transverse[e.index]*f.dn[e.index]*f.sx[e.index]*p).eval();
        else p=(-f.pairHop[e.index]*p).eval();
    }
    return (matrixExp(-(beta-last)*f.h0)*p).trace();
}
}
