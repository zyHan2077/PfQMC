#include "catalog/pxp_models.h"
#include <algorithm>
#include <cmath>
// Turner et al., arXiv:1806.10933v2, Eqs.(2),(5), section II.2.
// Their discussion includes the full 2^L spin Hilbert space and its frozen
// clusters. Here all of those states and every physical parity copy remain.
namespace catalog {
namespace {
Gaussian stateProjector(int modes,std::uint64_t state) {
    Gaussian p;p.trace=1.;p.g=MatType::Zero(2*modes,2*modes);
    for(int i=0;i<modes;++i){p.g(2*i,2*i+1)=((state>>i)&1)?-1.i:1.i;p.g(2*i+1,2*i)=-p.g(2*i,2*i+1);}return p;
}
Gaussian downProjector(int modes,int site) {
    Gaussian p;p.trace=std::ldexp(1.,modes-1);p.g=MatType::Zero(2*modes,2*modes);
    // P_down=na fixes ONLY the a occupation; every other mode is unprojected.
    p.g(4*site,4*site+1)=-1.i;p.g(4*site+1,4*site)=1.i;return p;
}
MatType embed(const MatType& a,int site,int sites) {
    int radix=a.rows(),dim=1,place=1;for(int i=0;i<sites;++i)dim*=radix;for(int i=0;i<site;++i)place*=radix;
    MatType out=MatType::Zero(dim,dim);
    for(int s=0;s<dim;++s){int before=(s/place)%radix;for(int after=0;after<radix;++after)out(s+(after-before)*place,s)=a(after,before);}return out;
}
void appendFlip(const PXPModel& m,int site,std::vector<Gaussian>& factors,DataType& scalar) {
    auto half=gaussian(m.halfRotation(site));factors.push_back(half);factors.push_back(half);
    // X=-i exp(pi gamma_1 gamma_2/2). The -i is essential even though
    // every physical X matrix element is real; the Spin lift is NOT discarded.
    scalar *= -1.i;
    if(m.phaseFlux&&site==0) {
        auto flux=gaussian(m.halfRotation(2,true));factors.push_back(flux);factors.push_back(flux);
        // 1-2n_2=2na_2-1=+i gamma_0 gamma_1.
        scalar *= 1.i;
    }
}
}
PXPModel PXPModel::chain(int n,double amplitude) {
    PXPModel m;m.sites=n;m.omega.assign(n,amplitude);m.validate();return m;
}
void PXPModel::validate() const {
    require(sites>=2&&sites<=12&&int(omega.size())==sites,"PXP requires 2..12 sites and one amplitude/site");
    require(!phaseFlux||sites>=3,"PXP flux control needs at least three sites");
    for(double x:omega)require(std::isfinite(x),"Nonfinite PXP amplitude");
    require(std::isfinite(repulsion)&&std::isfinite(detuning),"Nonfinite PXP diagonal coupling");
}
bool PXPModel::signFreeParameters() const {validate();return !phaseFlux&&std::all_of(omega.begin(),omega.end(),[](double x){return x>=0.;});}
std::uint64_t PXPModel::occupation(std::uint64_t excitation,std::uint64_t odd,int n) {
    require(n>=2&&n<=12&&excitation<(std::uint64_t(1)<<n),"PXP spin state out of range");
    return JQModel::occupation(((std::uint64_t(1)<<n)-1)^excitation,odd,n);
}
std::uint64_t PXPModel::excitations(std::uint64_t state) const {
    validate();require(state<(std::uint64_t(1)<<modes()),"PXP occupation out of range");
    std::uint64_t out=0;for(int i=0;i<sites;++i)out|=(1-((state>>(2*i))&1))<<i;return out;
}
std::vector<int> PXPModel::enabledSites(std::uint64_t state) const {
    auto bits=excitations(state);std::vector<int> out;
    for(int i=0;i<sites;++i)if((i==0||!((bits>>(i-1))&1))&&(i+1==sites||!((bits>>(i+1))&1)))out.push_back(i);return out;
}
std::uint64_t PXPModel::flipped(std::uint64_t state,int site) const {
    auto allowed=enabledSites(state);require(std::find(allowed.begin(),allowed.end(),site)!=allowed.end(),"PXP flip is blocked");
    return state^(std::uint64_t(3)<<(2*site));
}
double PXPModel::diagonalEnergy(std::uint64_t state) const {
    auto bits=excitations(state);double e=detuning*__builtin_popcountll(bits);
    for(int i=0;i+1<sites;++i)if(((bits>>i)&3)==3)e+=repulsion;return e;
}
bool PXPModel::valid(const PXPPath& p) const {
    validate();if(!(p.beta>0.)||!std::isfinite(p.beta)||p.states.empty()||p.states.size()!=p.events.size()+1||p.states.front()!=p.states.back())return false;
    for(auto s:p.states)if(s>=(std::uint64_t(1)<<modes()))return false;double last=0.;
    for(size_t k=0;k<p.events.size();++k){auto e=p.events[k];if(!std::isfinite(e.time)||e.time<last||e.time>p.beta||e.site<0||e.site>=sites)return false;last=e.time;
        auto allowed=enabledSites(p.states[k]);if(std::find(allowed.begin(),allowed.end(),e.site)==allowed.end()||flipped(p.states[k],e.site)!=p.states[k+1])return false;}return true;
}
MatType PXPModel::halfRotation(int site,bool flux) const {
    validate();require(site>=0&&site<sites,"PXP rotation site out of range");
    MatType a=MatType::Zero(2*modes(),2*modes());int p=4*site+(flux?0:1),q=4*site+(flux?1:2);
    a(p,q)=std::acos(-1.)/2.;a(q,p)=-a(p,q);return a;
}
DataType PXPModel::weight(const PXPPath& p) const {
    if(!valid(p))return 0.;std::vector<Gaussian> factors;DataType scalar=1.;double last=0.;
    for(size_t k=0;k<p.events.size();++k){auto e=p.events[k];scalar*=omega[e.site]*std::exp(-(e.time-last)*diagonalEnergy(p.states[k]));last=e.time;
        factors.push_back(stateProjector(modes(),p.states[k]));appendFlip(*this,e.site,factors,scalar);}
    scalar*=std::exp(-(p.beta-last)*diagonalEnergy(p.states.back()));factors.push_back(stateProjector(modes(),p.states.back()));
    return scalar*traceGaussianFactors(factors);
}
DataType PXPModel::reference(const PXPPath& p) const {
    if(!valid(p))return 0.;double out=1.,last=0.;
    for(size_t k=0;k<p.events.size();++k){auto e=p.events[k];out*=omega[e.site]*std::exp(-(e.time-last)*diagonalEnergy(p.states[k]));last=e.time;
        if(phaseFlux&&e.site==0)out*=1-2*int((excitations(p.states[k])>>2)&1);}
    return out*std::exp(-(p.beta-last)*diagonalEnergy(p.states.back()));
}
DataType PXPModel::kineticWordWeight(const std::vector<int>& labels,double beta) const {
    validate();require(beta>=0&&std::isfinite(beta),"Invalid PXP word beta");
    if(labels.empty())return std::ldexp(1.,modes());std::vector<Gaussian> factors;DataType scalar=1.;int order=0;
    for(int site:labels){require(site>=0&&site<sites,"Invalid PXP word label");scalar*=beta*omega[site]/(++order);
        if(site>0)factors.push_back(downProjector(modes(),site-1));if(site+1<sites)factors.push_back(downProjector(modes(),site+1));
        appendFlip(*this,site,factors,scalar);}
    return scalar*traceGaussianFactors(factors);
}
MatType PXPModel::termMatrix(int site,bool full) const {
    validate();require(site>=0&&site<sites&&sites<=(full?4:8),"PXP matrix reference size/index bound");
    int radix=full?4:2,dim=1<<(sites*(full?2:1));MatType localN,localX;
    if(full){localN=.5*MatType::Identity(4,4)+JQModel::localSpin(2);localX=2.*JQModel::localSpin(0);}
    else {localN=MatType::Zero(2,2);localN(1,1)=1.;localX=MatType::Zero(2,2);localX(0,1)=localX(1,0)=1.;}
    MatType id=MatType::Identity(dim,dim),out=embed(localX,site,sites),empty=MatType::Identity(radix,radix)-localN;
    if(site>0)out=out*embed(empty,site-1,sites);if(site+1<sites)out=out*embed(empty,site+1,sites);
    if(phaseFlux&&site==0)out=out*(id-2.*embed(localN,2,sites));return out;
}
MatType PXPModel::exactHamiltonian(bool full) const {
    validate();require(sites<=(full?4:8),"PXP Hamiltonian reference size bound");
    int dim=1<<(sites*(full?2:1));MatType h=MatType::Zero(dim,dim),n;
    if(full)n=.5*MatType::Identity(4,4)+JQModel::localSpin(2);else{n=MatType::Zero(2,2);n(1,1)=1.;}
    for(int i=0;i<sites;++i){h-=omega[i]*termMatrix(i,full);auto ni=embed(n,i,sites);h+=detuning*ni;if(i+1<sites)h+=repulsion*ni*embed(n,i+1,sites);}return h;
}
}
