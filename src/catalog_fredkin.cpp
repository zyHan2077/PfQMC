#include "catalog/fredkin_models.h"
#include <algorithm>
#include <cmath>
namespace catalog {
namespace {
MatType fredkinPauli(int axis) {
    MatType s=MatType::Zero(2,2);if(axis==0)s(0,1)=s(1,0)=.5;else if(axis==1){s(0,1)=-.5i;s(1,0)=.5i;}else{s(0,0)=.5;s(1,1)=-.5;}return s;
}
MatType fredkinEmbed(const MatType& a,int site,int sites) {
    int radix=a.rows(),dim=1,place=1;for(int r=0;r<sites;++r)dim*=radix;for(int r=0;r<site;++r)place*=radix;
    MatType out=MatType::Zero(dim,dim);for(int s=0;s<dim;++s){int old=(s/place)%radix;for(int t=0;t<radix;++t)out(s+(t-old)*place,s)=a(t,old);}return out;
}
Gaussian fredkinProjector(int modes,std::uint64_t s) {
    Gaussian p;p.trace=1.;p.g=MatType::Zero(2*modes,2*modes);for(int i=0;i<modes;++i){p.g(2*i,2*i+1)=((s>>i)&1)?-1.i:1.i;p.g(2*i+1,2*i)=-p.g(2*i,2*i+1);}return p;
}
bool enabled(std::uint64_t s,const FredkinTerm& t) {
    return (((s>>(2*t.control))&1)==(t.up?0:1)) && (((s>>(2*t.first))&1)!=((s>>(2*t.second))&1));
}
std::uint64_t swapped(std::uint64_t s,const FredkinTerm& t){return s^(std::uint64_t(3)<<(2*t.first))^(std::uint64_t(3)<<(2*t.second));}
}
FredkinModel FredkinModel::chain(int n,double left,double right) {
    FredkinModel m;m.sites=n;m.fields.assign(n,0.);m.terms.clear();for(int i=0;i+2<n;++i)m.terms.push_back({i,i+1,i+2,true,left});for(int i=0;i+2<n;++i)m.terms.push_back({i+2,i,i+1,false,right});m.validate();return m;
}
void FredkinModel::validate() const {
    require(sites>=3 && sites<=12 && int(fields.size())==sites && !terms.empty(),"Fredkin requires 3..12 sites");
    require(std::isfinite(boundaryLeft) && std::isfinite(boundaryRight),"Nonfinite Fredkin boundary");for(auto x:fields)require(std::isfinite(x),"Nonfinite Fredkin field");
    for(const auto& t:terms)require(t.control>=0 && t.first>=0 && t.second>=0 && t.control<sites && t.first<sites && t.second<sites && t.control!=t.first && t.control!=t.second && t.first!=t.second && std::isfinite(t.coupling),"Invalid Fredkin conditional term");
}
bool FredkinModel::signFreeParameters() const {validate();for(const auto& t:terms)if(t.coupling<0)return false;return true;}
std::uint64_t FredkinModel::spinState(std::uint64_t s) const {std::uint64_t out=0;for(int r=0;r<sites;++r)out|=((s>>(2*r))&1)<<r;return out;}
std::vector<std::uint64_t> FredkinModel::nextStates(std::uint64_t s,int term) const {
    require(term>=-1 && term<int(terms.size()) && s<(std::uint64_t(1)<<modes()),"Invalid Fredkin state/term");std::vector<std::uint64_t> out;
    if(term>=0){if(enabled(s,terms[term]))out.push_back(swapped(s,terms[term]));return out;}
    for(const auto& t:terms)if(enabled(s,t)){auto after=swapped(s,t);if(std::find(out.begin(),out.end(),after)==out.end())out.push_back(after);}return out;
}
double FredkinModel::diagonalEnergy(std::uint64_t s) const {
    double e=boundaryLeft*double(s&1)+boundaryRight*(1.-double((s>>(2*(sites-1)))&1));
    for(int r=0;r<sites;++r)e+=fields[r]*(.5-double((s>>(2*r))&1));
    for(const auto& t:terms)if(enabled(s,t))e+=t.coupling/2.;return e;
}
double FredkinModel::transitionCoefficient(std::uint64_t before,std::uint64_t after,int term) const {
    require(term>=-1 && term<int(terms.size()),"Invalid Fredkin vertex label");double c=0.;
    if(term>=0)return enabled(before,terms[term]) && swapped(before,terms[term])==after ? terms[term].coupling/2.:0.;
    for(const auto& t:terms)if(enabled(before,t) && swapped(before,t)==after)c+=t.coupling/2.;return c;
}
bool FredkinModel::valid(const FredkinPath& p) const {
    validate();if(!(p.beta>0) || !std::isfinite(p.beta) || p.states.empty() || p.states.size()!=p.events.size()+1 || p.states.front()!=p.states.back())return false;
    for(auto s:p.states)if(s>=(std::uint64_t(1)<<modes()))return false;double last=0.;
    for(int k=0;k<int(p.events.size());++k){auto e=p.events[k];if(!std::isfinite(e.time) || e.time<last || e.time>p.beta || e.term<-1 || e.term>=int(terms.size()))return false;last=e.time;auto next=nextStates(p.states[k],e.term);if(std::find(next.begin(),next.end(),p.states[k+1])==next.end())return false;}return true;
}
MatType FredkinModel::vertexHalfRotation(std::uint64_t before,std::uint64_t after) const {
    MatType a=MatType::Zero(2*modes(),2*modes());int count=0;
    for(int r=0;r<sites;++r){auto diff=((before^after)>>(2*r))&3;if(!diff)continue;require(diff==3,"Fredkin transition changed a parity copy");++count;
        double direction=((before>>(2*r))&1)?1.:-1.;a(4*r,4*r+2)=-direction*.5*std::acos(-1.);a(4*r+2,4*r)=-a(4*r,4*r+2);}
    require(count==2,"Fredkin event must exchange two spins");return a;
}
DataType FredkinModel::weight(const FredkinPath& p) const {
    if(!valid(p))return 0.;double scalar=1.,last=0.;std::vector<Gaussian> factors;
    for(int k=0;k<int(p.events.size());++k){auto e=p.events[k];scalar*=transitionCoefficient(p.states[k],p.states[k+1],e.term)*std::exp(-(e.time-last)*diagonalEnergy(p.states[k]));last=e.time;
        factors.push_back(fredkinProjector(modes(),p.states[k]));auto half=gaussian(vertexHalfRotation(p.states[k],p.states[k+1]));factors.push_back(half);factors.push_back(half);}
    scalar*=std::exp(-(p.beta-last)*diagonalEnergy(p.states.back()));factors.push_back(fredkinProjector(modes(),p.states.back()));return scalar*traceGaussianFactors(factors);
}
DataType FredkinModel::reference(const FredkinPath& p) const {
    if(!valid(p))return 0.;double w=1.,last=0.;for(int k=0;k<int(p.events.size());++k){auto e=p.events[k];w*=transitionCoefficient(p.states[k],p.states[k+1],e.term)*std::exp(-(e.time-last)*diagonalEnergy(p.states[k]));last=e.time;}return w*std::exp(-(p.beta-last)*diagonalEnergy(p.states.back()));
}
MatType FredkinModel::termMatrix(int term,bool full) const {
    validate();require(term>=0 && term<int(terms.size()) && sites<=(full?4:8),"Fredkin oracle size/index bound");int dim=1<<(sites*(full?2:1));auto id=MatType::Identity(dim,dim).eval();auto t=terms[term];MatType p=.25*id;
    for(int axis=0;axis<3;++axis){auto s=full?JQModel::localSpin(axis):fredkinPauli(axis);p-=fredkinEmbed(s,t.first,sites)*fredkinEmbed(s,t.second,sites);}
    auto z=full?JQModel::localSpin(2):fredkinPauli(2);MatType control=.5*id+(t.up?1.:-1.)*fredkinEmbed(z,t.control,sites);return control*p;
}
MatType FredkinModel::exactHamiltonian(bool full) const {
    validate();require(sites<=(full?4:8),"Fredkin Hamiltonian oracle size bound");int dim=1<<(sites*(full?2:1));auto id=MatType::Identity(dim,dim).eval();MatType h=MatType::Zero(dim,dim);
    for(int t=0;t<int(terms.size());++t)h+=terms[t].coupling*termMatrix(t,full);
    auto z=full?JQModel::localSpin(2):fredkinPauli(2);h+=boundaryLeft*(.5*id-fredkinEmbed(z,0,sites))+boundaryRight*(.5*id+fredkinEmbed(z,sites-1,sites));
    for(int r=0;r<sites;++r)h+=fields[r]*fredkinEmbed(z,r,sites);return h;
}
}
