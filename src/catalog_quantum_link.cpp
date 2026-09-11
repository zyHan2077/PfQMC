#include "catalog/quantum_link.h"
#include <algorithm>
#include <cmath>
namespace catalog {
Gaussian linkOccupationProjector(int n,std::uint64_t s) {
 require(n>0 && n<63 && s<(std::uint64_t(1)<<n),"Link projector dimension");
 Gaussian p;p.trace=1.;p.g=MatType::Zero(2*n,2*n);
 for(int i=0;i<n;++i){p.g(2*i,2*i+1)=(s&(std::uint64_t(1)<<i))?-1.i:1.i;p.g(2*i+1,2*i)=-p.g(2*i,2*i+1);}
 return p;
}
std::vector<std::pair<int,int>> QuantumLink::edges() const {
 require(sites>=2 && sites<=8 && (!periodic || sites>=3),"Quantum link: path/cycle with 2..8 sites, cycle>=3");
 std::vector<std::pair<int,int>> e;for(int i=0;i+1<sites;++i)e.emplace_back(i,i+1);if(periodic)e.emplace_back(sites-1,0);return e;
}
void QuantumLink::validate() const {
 (void)edges();require(std::isfinite(hopping) && hopping>=0 && std::isfinite(repulsion) && std::isfinite(electricDensity) && std::isfinite(chemical),"Quantum link finite real parameters; t>=0");
}
double QuantumLink::diagonalEnergy(std::uint64_t s) const {
 double energy=0.;for(int i=0;i<sites;++i)energy+=chemical*(double((s>>i)&1)-.5);
 auto e=edges();for(int b=0;b<int(e.size());++b){auto [i,j]=e[b];double ni=(s>>i)&1,nj=(s>>j)&1,d=int((s>>(sites+2*b))&1)+int((s>>(sites+2*b+1))&1)-1;energy+=repulsion*(ni-.5)*(nj-.5)+.5*electricDensity*d*(nj-ni);}
 return energy;
}
std::vector<std::uint64_t> QuantumLink::nextStates(std::uint64_t s) const {
 auto e=edges();require(s<(std::uint64_t(1)<<modes()),"Quantum link occupation outside Fock space");std::vector<std::uint64_t> out;
 for(int b=0;b<int(e.size());++b){auto [i,j]=e[b];int ni=(s>>i)&1,nj=(s>>j)&1,pair=(s>>(sites+2*b))&3;
  if((ni==1 && nj==0 && pair==0)||(ni==0 && nj==1 && pair==3))out.push_back(s^(std::uint64_t(1)<<i)^(std::uint64_t(1)<<j)^(std::uint64_t(3)<<(sites+2*b)));
 }return out;
}
bool QuantumLink::valid(const LinkPath& p) const {
 validate();if(!(p.beta>0) || !std::isfinite(p.beta) || p.states.empty() || p.times.size()+1!=p.states.size() || p.states.front()!=p.states.back())return false;
 double last=0.;for(auto s:p.states)if(s>=(std::uint64_t(1)<<modes()))return false;
 for(int k=0;k<int(p.times.size());++k){if(!std::isfinite(p.times[k]) || p.times[k]<last || p.times[k]>p.beta)return false;last=p.times[k];auto next=nextStates(p.states[k]);if(std::find(next.begin(),next.end(),p.states[k+1])==next.end())return false;}return true;
}
MatType QuantumLink::hopRotation(std::uint64_t before,std::uint64_t after) const {
 int n=modes();auto e=edges();auto next=nextStates(before);require(std::find(next.begin(),next.end(),after)!=next.end(),"Invalid quantum-link hop");
 for(int b=0;b<int(e.size());++b)if(((before^after)>>(sites+2*b))&3) {
  auto [i,j]=e[b];bool right=(before>>i)&1;int source=right?i:j,dest=right?j:i;
  MatType k=MatType::Zero(n,n);k(dest,source)=1.;k(source,dest)=-1.;
  auto pair=pairMajorana(n,sites+2*b,sites+2*b+1,1.,false);
  return .5*std::acos(-1.)*(numberMajorana(k)+(right?1.:-1.)*(pair-pair.adjoint()));
 }throw std::runtime_error("Quantum link did not change its link");
}
DataType QuantumLink::weight(const LinkPath& p) const {
 if(!valid(p))return 0.;std::vector<Gaussian> factors;double scalar=1.,last=0.;int n=modes();
 for(int k=0;k<int(p.times.size());++k){scalar*=hopping*std::exp(-(p.times[k]-last)*diagonalEnergy(p.states[k]));last=p.times[k];factors.push_back(linkOccupationProjector(n,p.states[k]));factors.push_back(gaussian(hopRotation(p.states[k],p.states[k+1])));}
 scalar*=std::exp(-(p.beta-last)*diagonalEnergy(p.states.back()));factors.push_back(linkOccupationProjector(n,p.states.back()));return scalar*traceGaussianFactors(factors);
}
DataType QuantumLink::reference(const LinkPath& p) const {
 if(!valid(p))return 0.;double scalar=1.,last=0.;
 for(int k=0;k<int(p.times.size());++k){auto diff=(p.states[k]^p.states[k+1])&((std::uint64_t(1)<<sites)-1);int a=__builtin_ctzll(diff),b=63-__builtin_clzll(diff);auto between=((std::uint64_t(1)<<b)-1)^((std::uint64_t(1)<<(a+1))-1);int sign=__builtin_popcountll(p.states[k]&between)%2?-1:1;scalar*=sign*hopping*std::exp(-(p.times[k]-last)*diagonalEnergy(p.states[k]));last=p.times[k];}
 return scalar*std::exp(-(p.beta-last)*diagonalEnergy(p.states.back()));
}
MatType QuantumLink::exactHamiltonian() const {
 validate();int n=modes();require(n<=9,"Quantum link full Fock oracle requires <=9 modes");int dim=1<<n;MatType id=MatType::Identity(dim,dim);std::vector<MatType> annih,number;
 for(int i=0;i<n;++i){MatType a=MatType::Zero(dim,dim);for(unsigned s=0;s<unsigned(dim);++s)if(s&(1u<<i))a(s^(1u<<i),s)=__builtin_popcount(s&((1u<<i)-1))%2?-1.:1.;annih.push_back(a);number.push_back(a.adjoint()*a);}
 MatType h=MatType::Zero(dim,dim);for(int i=0;i<sites;++i)h+=chemical*(number[i]-.5*id);
 auto e=edges();for(int b=0;b<int(e.size());++b){auto [i,j]=e[b];int a=sites+2*b;MatType B=annih[a].adjoint()*annih[a+1].adjoint();MatType move=annih[j].adjoint()*B*annih[i];h-=hopping*(move+move.adjoint());h+=repulsion*(number[i]-.5*id)*(number[j]-.5*id);h+=.5*electricDensity*(number[a]+number[a+1]-id)*(number[j]-number[i]);}
 return h;
}
}
