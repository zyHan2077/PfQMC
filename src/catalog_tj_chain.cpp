#include "catalog/tj_chain.h"
#include <algorithm>
#include <cmath>
namespace catalog {
namespace {int occ(std::uint64_t s,int i){return (s>>(2*i))&3;}int count(int a){return __builtin_popcount(unsigned(a));}
int hopSign(std::uint64_t s,int from,int to){int sign=(__builtin_popcountll(s&((std::uint64_t(1)<<from)-1))%2)?-1:1;s^=std::uint64_t(1)<<from;return sign*((__builtin_popcountll(s&((std::uint64_t(1)<<to)-1))%2)?-1:1);}
}
void TJChain::validate()const{require(sites>=2&&sites<=15&&(!periodic||sites>=3),"TJ size must be 2..15 (ring >=3)");for(double x:{hopping,exchange,doublon,density,chemical})require(std::isfinite(x),"TJ finite real parameters required");}
bool TJChain::signFreeParameters()const{validate();return !periodic&&hopping>=0&&exchange>=0;}
std::vector<std::pair<int,int>> TJChain::edges()const{validate();std::vector<std::pair<int,int>> out;for(int i=0;i+1<sites;++i)out.emplace_back(i,i+1);if(periodic)out.emplace_back(sites-1,0);return out;}
std::vector<std::uint64_t> TJChain::nextStates(std::uint64_t s)const{validate();require(s<(std::uint64_t(1)<<modes()),"TJ state outside Fock space");std::vector<std::uint64_t> out;for(auto [i,j]:edges()){int a=occ(s,i),b=occ(s,j);if((a==0&&(b==1||b==2))||(b==0&&(a==1||a==2))||(a+b==3&&a*b==2))out.push_back(s^(std::uint64_t(a^b)<<(2*i))^(std::uint64_t(a^b)<<(2*j)));}return out;}
bool TJChain::valid(const LinkPath& p)const{validate();if(!std::isfinite(p.beta)||p.beta<=0||p.states.empty()||p.states.front()!=p.states.back()||p.states.size()!=p.times.size()+1)return false;for(auto s:p.states)if(s>=(std::uint64_t(1)<<modes()))return false;double last=0.;for(int k=0;k<int(p.times.size());++k){if(!std::isfinite(p.times[k])||p.times[k]<last||p.times[k]>p.beta)return false;last=p.times[k];auto ns=nextStates(p.states[k]);if(std::find(ns.begin(),ns.end(),p.states[k+1])==ns.end())return false;}return true;}
int TJChain::squeezedPhase(std::uint64_t s)const{validate();require(s<(std::uint64_t(1)<<modes()),"TJ state outside Fock space");int rank=0,exponent=0;for(int i=0;i<sites;++i){int a=occ(s,i);if(a==3){rank=0;continue;}if(a){++rank;if(a==2)exponent+=rank;}}return exponent%2?-1:1;}
double TJChain::diagonalEnergy(std::uint64_t s)const{validate();require(s<(std::uint64_t(1)<<modes()),"TJ state outside Fock space");double e=0.;for(int i=0;i<sites;++i)e+=doublon*(occ(s,i)==3)-chemical*count(occ(s,i));for(auto [i,j]:edges()){int a=occ(s,i),b=occ(s,j),za=(a==1)-(a==2),zb=(b==1)-(b==2);e+=exchange*(za*zb-count(a)*count(b))/4.+density*count(a)*count(b);}return e;}
double TJChain::vertexScalar(std::uint64_t a,std::uint64_t b)const{auto ns=nextStates(a);require(std::find(ns.begin(),ns.end(),b)!=ns.end(),"TJ invalid transition");return __builtin_popcountll(a^b)==2?hopping:-exchange/2.;}
MatType TJChain::rotation(std::uint64_t a,std::uint64_t b)const{(void)vertexScalar(a,b);MatType k=MatType::Zero(modes(),modes());auto add=[&](int from,int to){k(to,from)+=M_PI/2;k(from,to)-=M_PI/2;};if(__builtin_popcountll(a^b)==2){int from=__builtin_ctzll(a&~b),to=__builtin_ctzll(b&~a);add(from,to);}else for(int i=0;i<sites;++i)if(occ(a,i)!=occ(b,i))add(2*i+(occ(a,i)==2),2*i+(occ(b,i)==2));return numberMajorana(k);}
DataType TJChain::weight(const LinkPath& p)const{if(!valid(p))return 0.;std::vector<Gaussian> factors;double z=1.,last=0.;for(int k=0;k<int(p.times.size());++k){z*=vertexScalar(p.states[k],p.states[k+1])*std::exp(-(p.times[k]-last)*diagonalEnergy(p.states[k]));last=p.times[k];factors.push_back(linkOccupationProjector(modes(),p.states[k]));factors.push_back(gaussian(rotation(p.states[k],p.states[k+1])));}factors.push_back(linkOccupationProjector(modes(),p.states.back()));return z*std::exp(-(p.beta-last)*diagonalEnergy(p.states.back()))*traceGaussianFactors(factors);}
DataType TJChain::reference(const LinkPath& p)const{if(!valid(p))return 0.;double z=1.,last=0.;for(int k=0;k<int(p.times.size());++k){auto a=p.states[k],b=p.states[k+1];int phase=1;if(__builtin_popcountll(a^b)==2)phase=hopSign(a,__builtin_ctzll(a&~b),__builtin_ctzll(b&~a));z*=phase*vertexScalar(a,b)*std::exp(-(p.times[k]-last)*diagonalEnergy(a));last=p.times[k];}return z*std::exp(-(p.beta-last)*diagonalEnergy(p.states.back()));}
MatType TJChain::exactHamiltonian()const{
 validate();require(modes()<=10,"TJ full Fock oracle limited to ten modes");int dim=1<<modes();MatType h=MatType::Zero(dim,dim);
 auto apply=[&](unsigned s,const std::vector<std::pair<int,bool>>& ops,double c){unsigned q=s;for(auto it=ops.rbegin();it!=ops.rend();++it){int m=it->first;bool create=it->second;if(bool((q>>m)&1)==create)return;if(__builtin_popcount(q&((1u<<m)-1))%2)c=-c;q^=1u<<m;}h(q,s)+=c;};
 for(unsigned s=0;s<unsigned(dim);++s){double e=0.;for(int i=0;i<sites;++i){int up=(s>>(2*i))&1,dn=(s>>(2*i+1))&1;e+=doublon*up*dn-chemical*(up+dn);}for(auto [i,j]:edges()){
  int ui=(s>>(2*i))&1,di=(s>>(2*i+1))&1,uj=(s>>(2*j))&1,dj=(s>>(2*j+1))&1;
  e+=exchange*((ui-di)*(uj-dj)-(ui+di)*(uj+dj))/4.+density*(ui+di)*(uj+dj);
  for(int spin=0;spin<2;++spin)if(!((s>>(2*i+1-spin))&1)&&!((s>>(2*j+1-spin))&1)){
   apply(s,{{2*i+spin,true},{2*j+spin,false}},-hopping);apply(s,{{2*j+spin,true},{2*i+spin,false}},-hopping);
  }
  apply(s,{{2*i,true},{2*i+1,false},{2*j+1,true},{2*j,false}},exchange/2.);
  apply(s,{{2*i+1,true},{2*i,false},{2*j,true},{2*j+1,false}},exchange/2.);
 }h(s,s)+=e;}return h;
}
}
