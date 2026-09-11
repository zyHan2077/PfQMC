#include "catalog/color_projector.h"
#include <algorithm>
#include <cmath>
#include <set>

namespace catalog {
namespace {
Gaussian colorOccupationProjector(int n,std::uint64_t s) {
    Gaussian p;p.trace=1.;p.g=MatType::Zero(2*n,2*n);
    for(int i=0;i<n;++i){p.g(2*i,2*i+1)=((s>>i)&1)?-1.i:1.i;p.g(2*i+1,2*i)=-p.g(2*i,2*i+1);}
    return p;
}
MatType localBilinear(int dest,int source) {
    MatType a=MatType::Zero(8,8);
    for(unsigned s=0;s<8;++s)if((s>>source)&1) {
        unsigned t=s^(1u<<source);double sign=__builtin_popcount(s&((1u<<source)-1))%2?-1.:1.;
        if(!((t>>dest)&1)){sign*=__builtin_popcount(t&((1u<<dest)-1))%2?-1.:1.;a(t^(1u<<dest),s)=sign;}
    }return a;
}
MatType product(const MatType& a,const MatType& b) {
    MatType out(64,64);
    for(int s=0;s<64;++s)for(int t=0;t<64;++t)out(t,s)=a(t&7,s&7)*b(t>>3,s>>3);
    return out;
}
int statePhase(std::uint64_t s,int sites) {
    int sign=1;for(int r=0;r<sites;++r)sign*=ColorProjector::occupationPhase((s>>(3*r))&7);return sign;
}
}
void ColorProjector::validate() const {
    require(sites>=2 && sites<=8,"Color projector supports 2..8 cells");
    std::set<std::pair<int,int>> found;
    for(auto e:edges){require(e.first>=0 && e.second>=0 && e.first<sites && e.second<sites && e.first!=e.second,"Color edge outside graph");if(e.first>e.second)std::swap(e.first,e.second);require(found.insert(e).second,"Duplicate color edge");}
    require(!edges.empty() && std::isfinite(singlet) && std::isfinite(permutation) && std::isfinite(charging) && std::isfinite(chemical) && std::isfinite(density) && std::isfinite(activeDensity),"Nonfinite color parameters");
}
bool ColorProjector::signFreeParameters() const {validate();return singlet>=0 && permutation>=0;}
int ColorProjector::color(unsigned s) {
    require(s<8,"Color local occupation outside Fock space");
    static const int colors[8]={-1,0,1,2,2,1,0,-1};return colors[s];
}
int ColorProjector::occupationPhase(unsigned s){require(s<8,"Color local occupation outside Fock space");return s==5?-1:1;}
unsigned ColorProjector::occupation(int n,int c){require((n==1 || n==2) && c>=0 && c<3,"Invalid active color");return n==1?(1u<<c):(7u^(1u<<c));}
double ColorProjector::diagonalEnergy(std::uint64_t s) const {
    double energy=0.;std::vector<double> q(sites);
    for(int r=0;r<sites;++r){q[r]=__builtin_popcountll((s>>(3*r))&7)-1.5;energy+=charging*q[r]*q[r]+chemical*q[r];}
    for(auto e:edges){double ni=q[e.first]+1.5,nj=q[e.second]+1.5;energy+=density*q[e.first]*q[e.second]+activeDensity*ni*(3-ni)*nj*(3-nj)/4.;}return energy;
}
double ColorProjector::coefficient(ColorVertex kind) const {
    require(kind==ColorVertex::Singlet || kind==ColorVertex::Permutation,"Unknown color vertex");
    return kind==ColorVertex::Singlet?singlet/3.:permutation;
}
std::vector<std::uint64_t> ColorProjector::nextStates(std::uint64_t before,int edge,ColorVertex kind) const {
    require(edge>=0 && edge<int(edges.size()) && before<(std::uint64_t(1)<<modes()),"Invalid color transition index");(void)coefficient(kind);
    auto [i,j]=edges[edge];unsigned si=(before>>(3*i))&7,sj=(before>>(3*j))&7;
    int ci=color(si),cj=color(sj);if(ci<0 || cj<0 || (kind==ColorVertex::Singlet && ci!=cj))return {};
    auto make=[&](int di,int dj){auto mask=(std::uint64_t(7)<<(3*i))|(std::uint64_t(7)<<(3*j));return (before&~mask)|(std::uint64_t(occupation(__builtin_popcount(si),di))<<(3*i))|(std::uint64_t(occupation(__builtin_popcount(sj),dj))<<(3*j));};
    if(kind==ColorVertex::Permutation)return {make(cj,ci)};
    return {make(0,0),make(1,1),make(2,2)};
}
bool ColorProjector::valid(const ColorPath& p) const {
    validate();if(!(p.beta>0) || !std::isfinite(p.beta) || p.states.empty() || p.states.size()!=p.events.size()+1 || p.states.front()!=p.states.back())return false;
    for(auto s:p.states)if(s>=(std::uint64_t(1)<<modes()))return false;
    double last=0.;for(int k=0;k<int(p.events.size());++k){auto e=p.events[k];if(!std::isfinite(e.time) || e.time<last || e.time>p.beta || e.edge<0 || e.edge>=int(edges.size()) || (e.kind!=ColorVertex::Singlet && e.kind!=ColorVertex::Permutation))return false;last=e.time;auto next=nextStates(p.states[k],e.edge,e.kind);if(std::find(next.begin(),next.end(),p.states[k+1])==next.end())return false;}return true;
}
MatType ColorProjector::vertexRotation(std::uint64_t before,std::uint64_t after) const {
    MatType k=MatType::Zero(modes(),modes());
    for(int r=0;r<sites;++r){unsigned a=(before>>(3*r))&7,b=(after>>(3*r))&7;if(a==b)continue;
        int n=__builtin_popcount(a);require(color(a)>=0 && color(b)>=0 && __builtin_popcount(b)==n && __builtin_popcount(a^b)==2,"Invalid color rotation");
        int source=__builtin_ctz(a&~b),dest=__builtin_ctz(b&~a);double direction=n==1?1.:-1.;
        k(3*r+dest,3*r+source)=direction;k(3*r+source,3*r+dest)=-direction;
    }return .5*std::acos(-1.)*numberMajorana(k);
}
DataType ColorProjector::weight(const ColorPath& p) const {
    if(!valid(p))return 0.;std::vector<Gaussian> factors;double scalar=1.,last=0.;
    for(int k=0;k<int(p.events.size());++k){auto e=p.events[k];scalar*=coefficient(e.kind)*std::exp(-(e.time-last)*diagonalEnergy(p.states[k]));last=e.time;factors.push_back(colorOccupationProjector(modes(),p.states[k]));factors.push_back(gaussian(vertexRotation(p.states[k],p.states[k+1])));}
    scalar*=std::exp(-(p.beta-last)*diagonalEnergy(p.states.back()));factors.push_back(colorOccupationProjector(modes(),p.states.back()));return scalar*traceGaussianFactors(factors);
}
DataType ColorProjector::reference(const ColorPath& p) const {
    if(!valid(p))return 0.;double scalar=1.,last=0.;
    for(int k=0;k<int(p.events.size());++k){auto e=p.events[k];scalar*=coefficient(e.kind)*statePhase(p.states[k],sites)*statePhase(p.states[k+1],sites)*std::exp(-(e.time-last)*diagonalEnergy(p.states[k]));last=e.time;}
    return scalar*std::exp(-(p.beta-last)*diagonalEnergy(p.states.back()));
}
MatType ColorProjector::localSpin(int axis) {
    require(axis>=0 && axis<3,"Unknown SO3 spin axis");int b=(axis+1)%3,c=(axis+2)%3;return -1.i*(localBilinear(b,c)-localBilinear(c,b));
}
MatType ColorProjector::localActive() {
    MatType a=MatType::Zero(8,8);for(unsigned s=0;s<8;++s){int n=__builtin_popcount(s);a(s,s)=n*(3-n)/2.;}return a;
}
MatType ColorProjector::localBondVertex(ColorVertex kind) {
    require(kind==ColorVertex::Singlet || kind==ColorVertex::Permutation,"Unknown color vertex");
    MatType dot=MatType::Zero(64,64);for(int a=0;a<3;++a)dot+=product(localSpin(a),localSpin(a));
    MatType aa=product(localActive(),localActive()),square=(dot*dot).eval();
    if(kind==ColorVertex::Singlet)return (square-aa)/3.;return square+dot-aa;
}
MatType ColorProjector::vertexMatrix(int edge,ColorVertex kind) const {
    validate();require(sites<=3 && edge>=0 && edge<int(edges.size()),"Color full Fock vertex supports <=3 cells");
    auto local=localBondVertex(kind);auto [i,j]=edges[edge];int dim=1<<modes();MatType out=MatType::Zero(dim,dim);
    std::uint64_t mask=(std::uint64_t(7)<<(3*i))|(std::uint64_t(7)<<(3*j));
    for(unsigned s=0;s<unsigned(dim);++s){int u=((s>>(3*i))&7)|(((s>>(3*j))&7)<<3);
        for(unsigned v=0;v<64;++v)if(std::abs(local(v,u))>0.){auto t=(s&~mask)|((std::uint64_t(v)&7)<<(3*i))|((std::uint64_t(v)>>3)<<(3*j));out(t,s)=local(v,u);}}
    return out;
}
MatType ColorProjector::exactHamiltonian() const {
    validate();require(sites<=3,"Color full Fock Hamiltonian supports <=3 cells");int dim=1<<modes();MatType h=MatType::Zero(dim,dim);
    for(int s=0;s<dim;++s)h(s,s)=diagonalEnergy(s);
    for(int e=0;e<int(edges.size());++e)h-=singlet*vertexMatrix(e,ColorVertex::Singlet)+permutation*vertexMatrix(e,ColorVertex::Permutation);
    return h;
}
}
