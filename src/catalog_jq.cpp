#include "catalog/jq_models.h"
#include <algorithm>
#include <cmath>
#include <queue>
#include <set>
namespace catalog {
namespace {
Gaussian jqProjector(int modes,std::uint64_t s) {
    Gaussian p;p.trace=1.;p.g=MatType::Zero(2*modes,2*modes);
    for(int i=0;i<modes;++i){p.g(2*i,2*i+1)=((s>>i)&1)?-1.i:1.i;p.g(2*i+1,2*i)=-p.g(2*i,2*i+1);}return p;
}
MatType embed(const MatType& a,int site,int sites,int radix) {
    int dim=1;for(int i=0;i<sites;++i)dim*=radix;int place=1;for(int i=0;i<site;++i)place*=radix;
    MatType out=MatType::Zero(dim,dim);for(int s=0;s<dim;++s){int old=(s/place)%radix;for(int t=0;t<radix;++t)out(s+(t-old)*place,s)=a(t,old);}return out;
}
MatType pauli(int axis) {
    MatType s=MatType::Zero(2,2);if(axis==0)s(0,1)=s(1,0)=.5;
    else if(axis==1){s(0,1)=-.5i;s(1,0)=.5i;}else{s(0,0)=.5;s(1,1)=-.5;}return s;
}
}
JQModel JQModel::rectangle(int lx,int ly,int order,double j,double q) {
    require(lx>=2 && ly>=2 && lx*ly<=12 && order>=2 && order<=3,"JQ rectangle/order outside reference limits");
    JQModel m;m.sites=lx*ly;m.fields.assign(m.sites,0.);m.terms.clear();auto id=[&](int x,int y){return x+lx*y;};
    for(int y=0;y<ly;++y)for(int x=0;x<lx;++x){if(x+1<lx)m.terms.push_back({{{id(x,y),id(x+1,y)}},j});if(y+1<ly)m.terms.push_back({{{id(x,y),id(x,y+1)}},j});}
    for(int y=0;y+order<=ly;++y)for(int x=0;x+1<lx;++x){JQTerm t{{},q};for(int a=0;a<order;++a)t.bonds.push_back({id(x,y+a),id(x+1,y+a)});m.terms.push_back(t);}
    for(int y=0;y+1<ly;++y)for(int x=0;x+order<=lx;++x){JQTerm t{{},q};for(int a=0;a<order;++a)t.bonds.push_back({id(x+a,y),id(x+a,y+1)});m.terms.push_back(t);}
    m.validate();return m;
}
void JQModel::validate() const {
    require(sites>=2 && sites<=12 && int(fields.size())==sites && !terms.empty(),"JQ model shape outside reference limits");for(double h:fields)require(std::isfinite(h),"Nonfinite JQ field");
    for(const auto& t:terms){require(!t.bonds.empty() && t.bonds.size()<=3 && std::isfinite(t.coupling),"Invalid JQ term");std::set<int> touched;
        for(auto [i,j]:t.bonds)require(i>=0 && j>=0 && i<sites && j<sites && i!=j && touched.insert(i).second && touched.insert(j).second,"JQ product bonds must be disjoint and inside graph");}
}
std::vector<int> JQModel::sublattices() const {
    validate();std::vector<std::vector<int>> graph(sites);for(const auto& t:terms)for(auto [i,j]:t.bonds){graph[i].push_back(j);graph[j].push_back(i);}
    std::vector<int> color(sites,0);for(int root=0;root<sites;++root)if(!color[root]){std::queue<int> q;q.push(root);color[root]=1;
        while(!q.empty()){int i=q.front();q.pop();for(int j:graph[i]){if(!color[j]){color[j]=-color[i];q.push(j);}else if(color[j]==color[i])return {};}}}return color;
}
bool JQModel::signFreeParameters() const {if(sublattices().empty())return false;for(const auto& t:terms)if(t.coupling<0)return false;return true;}
double JQModel::diagonalEnergy(std::uint64_t s) const {double e=0.;for(int r=0;r<sites;++r)e+=fields[r]*(.5-double((s>>(2*r))&1));return e;}
std::uint64_t JQModel::spinState(std::uint64_t s) const {std::uint64_t out=0;for(int r=0;r<sites;++r)out|=((s>>(2*r))&1)<<r;return out;}
std::uint64_t JQModel::occupation(std::uint64_t spin,std::uint64_t odd,int sites) {
    require(sites>=1 && sites<=12 && spin<(std::uint64_t(1)<<sites) && odd<(std::uint64_t(1)<<sites),"JQ spin/parity state outside space");std::uint64_t out=0;
    for(int r=0;r<sites;++r){auto a=(spin>>r)&1,b=a^((odd>>r)&1);out|=(a|(b<<1))<<(2*r);}return out;
}
std::vector<std::uint64_t> JQModel::nextStates(std::uint64_t s,int term) const {
    require(term>=0 && term<int(terms.size()) && s<(std::uint64_t(1)<<modes()),"JQ transition index outside model");const auto& bonds=terms[term].bonds;
    for(auto [i,j]:bonds)if(((s>>(2*i))&1)==((s>>(2*j))&1))return {};
    std::vector<std::uint64_t> out;for(unsigned choice=0;choice<(1u<<bonds.size());++choice){auto t=s;for(unsigned b=0;b<bonds.size();++b)if((choice>>b)&1){auto [i,j]=bonds[b];t^=(std::uint64_t(3)<<(2*i))|(std::uint64_t(3)<<(2*j));}out.push_back(t);}return out;
}
bool JQModel::valid(const JQPath& p) const {
    validate();if(!(p.beta>0) || !std::isfinite(p.beta) || p.states.empty() || p.states.size()!=p.events.size()+1 || p.states.front()!=p.states.back())return false;
    for(auto s:p.states)if(s>=(std::uint64_t(1)<<modes()))return false;double last=0.;
    for(int k=0;k<int(p.events.size());++k){auto e=p.events[k];if(!std::isfinite(e.time) || e.time<last || e.time>p.beta || e.term<0 || e.term>=int(terms.size()))return false;last=e.time;auto next=nextStates(p.states[k],e.term);if(std::find(next.begin(),next.end(),p.states[k+1])==next.end())return false;}return true;
}
MatType JQModel::vertexHalfRotation(std::uint64_t before,std::uint64_t after,int term) const {
    auto next=nextStates(before,term);require(std::find(next.begin(),next.end(),after)!=next.end(),"Invalid JQ resolved vertex");MatType a=MatType::Zero(2*modes(),2*modes());
    // Each swapped singlet bond has physical matrix element -1/2. The first
    // spin rotation carries that minus sign; it is not replaced by abs(weight).
    for(auto [i,j]:terms[term].bonds)if(((before^after)>>(2*i))&1)for(int r:{i,j}){
        double direction=((before>>(2*r))&1)?1.:-1.;if(r==i)direction=-direction;
        a(4*r,4*r+2)=-direction*.5*std::acos(-1.);a(4*r+2,4*r)=-a(4*r,4*r+2);
    }return a;
}
DataType JQModel::weight(const JQPath& p) const {
    if(!valid(p))return 0.;std::vector<Gaussian> factors;double scalar=1.,last=0.;
    for(int k=0;k<int(p.events.size());++k){auto e=p.events[k];scalar*=std::ldexp(terms[e.term].coupling,-int(terms[e.term].bonds.size()))*std::exp(-(e.time-last)*diagonalEnergy(p.states[k]));last=e.time;
        factors.push_back(jqProjector(modes(),p.states[k]));auto half=gaussian(vertexHalfRotation(p.states[k],p.states[k+1],e.term));factors.push_back(half);factors.push_back(half);}
    scalar*=std::exp(-(p.beta-last)*diagonalEnergy(p.states.back()));factors.push_back(jqProjector(modes(),p.states.back()));return scalar*traceGaussianFactors(factors);
}
DataType JQModel::reference(const JQPath& p) const {
    if(!valid(p))return 0.;double w=1.,last=0.;for(int k=0;k<int(p.events.size());++k){auto e=p.events[k];int flipped=__builtin_popcountll(spinState(p.states[k])^spinState(p.states[k+1]))/2;
        w*=(flipped%2?-1.:1.)*std::ldexp(terms[e.term].coupling,-int(terms[e.term].bonds.size()))*std::exp(-(e.time-last)*diagonalEnergy(p.states[k]));last=e.time;}return w*std::exp(-(p.beta-last)*diagonalEnergy(p.states.back()));
}
MatType JQModel::localSpin(int axis) {
    require(axis>=0 && axis<3,"Invalid JQ spin axis");std::vector<MatType> g;
    for(int i=0;i<2;++i){MatType a=MatType::Zero(4,4);for(unsigned s=0;s<4;++s)if((s>>i)&1)a(s^(1u<<i),s)=__builtin_popcount(s&((1u<<i)-1))%2?-1.:1.;g.push_back(a+a.adjoint());g.push_back(1.i*(a.adjoint()-a));}
    return -.5i*g[(axis+1)%3]*g[(axis+2)%3];
}
MatType JQModel::vertexMatrix(int term) const {
    validate();require(sites<=4 && term>=0 && term<int(terms.size()),"JQ full Fock supports <=4 sites");int dim=1<<modes();auto id=MatType::Identity(dim,dim).eval();MatType out=id;
    for(auto [i,j]:terms[term].bonds){MatType p=.25*id;for(int axis=0;axis<3;++axis){auto s=localSpin(axis);p-=embed(s,i,sites,4)*embed(s,j,sites,4);}out=(out*p).eval();}return out;
}
MatType JQModel::spinVertexMatrix(int term) const {
    validate();require(sites<=8 && term>=0 && term<int(terms.size()),"JQ Pauli oracle supports <=8 spins");int dim=1<<sites;auto id=MatType::Identity(dim,dim).eval();MatType out=id;
    for(auto [i,j]:terms[term].bonds){MatType p=.25*id;for(int axis=0;axis<3;++axis){auto s=pauli(axis);p-=embed(s,i,sites,2)*embed(s,j,sites,2);}out=(out*p).eval();}return out;
}
MatType JQModel::exactHamiltonian() const {validate();require(sites<=4,"JQ Fock size bound");MatType h=MatType::Zero(1<<modes(),1<<modes());for(int s=0;s<h.rows();++s)h(s,s)=diagonalEnergy(s);for(int t=0;t<int(terms.size());++t)h-=terms[t].coupling*vertexMatrix(t);return h;}
MatType JQModel::exactSpinHamiltonian() const {validate();require(sites<=8,"JQ spin size bound");MatType h=MatType::Zero(1<<sites,1<<sites);for(int s=0;s<h.rows();++s)h(s,s)=diagonalEnergy(occupation(s,0,sites));for(int t=0;t<int(terms.size());++t)h-=terms[t].coupling*spinVertexMatrix(t);return h;}
}
