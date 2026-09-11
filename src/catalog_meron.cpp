#include "catalog/meron_models.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <queue>
#include <set>
// Cluster resummation: Huffman and Chandrasekharan, arXiv:1605.07420v2,
// Section IV, Eq. (22). E-gated pair-fermion embedding retains all odd blockers:
// untouched sites supply four states instead of two, a positive extra factor.
// The entire H also has a rotated pair-worldline/conditional-Wei formulation;
// no claim is made that meron grouping is its unique sign-free representation.
namespace catalog {
namespace {
constexpr double pi=3.1415926535897932384626433832795;
MatType annihilation(int n,int i) {
    MatType out=MatType::Zero(1<<n,1<<n);
    for(unsigned s=0;s<(1u<<n);++s)if(s&(1u<<i))out(s^(1u<<i),s)=__builtin_popcount(s&((1u<<i)-1))%2?-1.:1.;
    return out;
}
MatType hopMatrix(const MeronModel& m) {
    MatType h=MatType::Zero(m.sites,m.sites);for(auto [i,j]:m.edges)h(i,j)=h(j,i)=m.hopping;return h;
}
Gaussian embed(const Gaussian& f,int offset,int total) {
    Gaussian out;out.g=MatType::Zero(2*total,2*total);out.g.block(2*offset,2*offset,f.g.rows(),f.g.cols())=f.g;
    out.trace=std::ldexp(1.,total-f.g.rows()/2)*f.trace;return out;
}
Gaussian projector(int sites,std::uint64_t occupation) {
    Gaussian f;f.trace=1.;f.g=MatType::Zero(4*sites,4*sites);
    for(int i=0;i<2*sites;++i){f.g(2*i,2*i+1)=DataType(0,(occupation&(std::uint64_t(1)<<i))?-1.:1.);f.g(2*i+1,2*i)=-f.g(2*i,2*i+1);}return f;
}
Gaussian rotation(int sites,int site,bool create) {
    auto a=pairMajorana(2*sites,2*site,2*site+1,1.,false);
    return gaussian((create?.5*pi:-.5*pi)*MatType(a-a.adjoint()));
}
std::vector<std::pair<std::uint64_t,double>> transitions(const MeronModel& m,const MeronEvent& event,std::uint64_t before) {
    if(event.kind==MeronKind::DensityBond)return {{before,1.}};
    if(event.kind==MeronKind::Monomer) {
        auto pair=(before>>(2*event.index))&3u;if(pair!=0 && pair!=3)return {};
        return {{before,.5},{before^(std::uint64_t(3)<<(2*event.index)),.5}};
    }
    auto [i,j]=m.edges[event.index];auto a=(before>>(2*i))&3u,b=(before>>(2*j))&3u;
    if((a!=0 && a!=3)||(b!=0 && b!=3)||a==b)return {};
    return {{before,.5},{before^(std::uint64_t(3)<<(2*i))^(std::uint64_t(3)<<(2*j)),-.5}};
}
DataType fermionReference(const MeronModel& m,double beta,const std::vector<MeronEvent>& events) {
    int n=m.sites;require(n<=6,"Meron fermion Fock oracle limited to six sites");int dim=1<<n;
    std::vector<MatType> c,dn;MatType q=MatType::Zero(dim,dim),id=MatType::Identity(dim,dim);
    for(int i=0;i<n;++i){c.push_back(annihilation(n,i));dn.push_back(c.back().adjoint()*c.back()-.5*id);}
    for(auto [i,j]:m.edges)q+=m.hopping*(c[i].adjoint()*c[j]+c[j].adjoint()*c[i]);
    auto sigma=m.sublattices();MatType product=id;double last=0.;
    for(const auto& e:events) {
        product=matrixExp((e.time-last)*q)*product;last=e.time;
        if(e.kind==MeronKind::DensityBond){auto [i,j]=m.edges[e.index];product=-dn[i]*dn[j]*product;}
        else if(e.kind==MeronKind::Monomer)product=double(sigma[e.index])*dn[e.index]*product;
    }
    return (matrixExp((beta-last)*q)*product).trace();
}
double vertexScalar(const MeronModel& m,const std::vector<MeronEvent>& events) {
    double scalar=1.;for(const auto& e:events)scalar*=e.kind==MeronKind::DensityBond?m.repulsion:e.kind==MeronKind::SingletBond?m.exchange:m.coupling;return scalar;
}
struct Physical {
    MatType h0;
    std::vector<MatType> density,singlet,monomer;
};
Physical physical(const MeronModel& m) {
    require(m.sites<=2,"Meron full 3*sites-mode Fock oracle limited to two sites");
    int n=3*m.sites,dim=1<<n;MatType id=MatType::Identity(dim,dim);
    std::vector<MatType> c,number,dn,E,sz,plus,p;
    for(int i=0;i<n;++i){c.push_back(annihilation(n,i));number.push_back(c.back().adjoint()*c.back());}
    for(int i=0;i<m.sites;++i) {
        int a=m.sites+2*i,b=a+1;dn.push_back(number[i]-.5*id);
        E.push_back(id-number[a]-number[b]+2.*number[a]*number[b]);
        sz.push_back(.5*(number[a]+number[b]-id));plus.push_back(c[a].adjoint()*c[b].adjoint());
        p.push_back(.5*(E.back()+plus.back()+plus.back().adjoint()));
    }
    Physical out;out.h0=MatType::Zero(dim,dim);auto sigma=m.sublattices();
    for(auto [i,j]:m.edges) {
        out.h0-=m.hopping*(c[i].adjoint()*c[j]+c[j].adjoint()*c[i]);
        out.density.push_back(m.repulsion*dn[i]*dn[j]);
        out.singlet.push_back(-m.exchange*(.25*E[i]*E[j]-sz[i]*sz[j]-.5*(plus[i]*plus[j].adjoint()+plus[i].adjoint()*plus[j])));
    }
    for(int i=0;i<m.sites;++i)out.monomer.push_back(-m.coupling*sigma[i]*dn[i]*p[i]);
    return out;
}
}
MeronModel MeronModel::square() {MeronModel m;m.sites=4;m.edges={{0,1},{1,2},{2,3},{3,0}};return m;}
std::vector<int> MeronModel::sublattices() const {
    require(sites>0 && sites<=20,"Meron occupation-mask bound is 20 sites");
    std::vector<std::vector<int>> adjacent(sites);std::set<std::pair<int,int>> seen;
    for(auto [i,j]:edges){require(i>=0&&j>=0&&i<sites&&j<sites&&i!=j,"Invalid meron edge");require(seen.insert(std::minmax(i,j)).second,"Duplicate meron edge");adjacent[i].push_back(j);adjacent[j].push_back(i);}
    std::vector<int> color(sites,0);std::queue<int> todo;
    for(int root=0;root<sites;++root)if(!color[root]){color[root]=1;todo.push(root);while(!todo.empty()){int i=todo.front();todo.pop();for(int j:adjacent[i]){if(!color[j]){color[j]=-color[i];todo.push(j);}require(color[j]==-color[i],"Meron graph must be bipartite");}}}
    return color;
}
void MeronModel::validate() const {
    (void)sublattices();require(std::isfinite(hopping)&&std::isfinite(repulsion)&&std::isfinite(exchange)&&std::isfinite(coupling)&&repulsion>=0&&exchange>=0&&coupling>=0,"Meron requires real t and nonnegative V,J,h");
}
bool MeronModel::validEvents(double beta,const std::vector<MeronEvent>& events) const {
    if(!std::isfinite(beta)||beta<=0)return false;double last=0.;
    for(const auto& e:events){if(!std::isfinite(e.time)||e.time<last||e.time>beta||e.index<0)return false;last=e.time;if(e.kind==MeronKind::Monomer){if(e.index>=sites)return false;}else if(e.kind==MeronKind::DensityBond||e.kind==MeronKind::SingletBond){if(e.index>=int(edges.size()))return false;}else return false;}return true;
}
MeronClusters MeronModel::clusters(const std::vector<MeronEvent>& events) const {
    validate();auto color=sublattices();std::vector<int> parent,nodeSite,ends,first(sites,-1),last(sites,-1);
    auto node=[&](int site,bool end){int i=parent.size();parent.push_back(i);nodeSite.push_back(site);ends.push_back(end?1:0);return i;};
    std::function<int(int)> find=[&](int i){return parent[i]==i?i:parent[i]=find(parent[i]);};
    auto unite=[&](int i,int j){parent[find(i)]=find(j);};
    auto legs=[&](int site,bool end){int in=node(site,end),out=node(site,end);if(last[site]>=0)unite(last[site],in);else first[site]=in;last[site]=out;return std::make_pair(in,out);};
    MeronClusters result;
    for(const auto& e:events) {
        if(e.kind==MeronKind::DensityBond)continue;++result.spinVertices;
        if(e.kind==MeronKind::Monomer){require(e.index>=0&&e.index<sites,"Invalid cluster monomer");legs(e.index,true);}
        else {require(e.index>=0&&e.index<int(edges.size()),"Invalid cluster bond");auto [i,j]=edges[e.index];auto a=legs(i,false),b=legs(j,false);unite(a.first,b.first);unite(a.second,b.second);}
    }
    for(int i=0;i<sites;++i){if(first[i]<0)++result.untouched;else unite(first[i],last[i]);}
    std::map<int,std::vector<int>> components;
    for(int i=0;i<int(parent.size());++i){int root=find(i);components[root];if(ends[i])components[root].push_back(color[nodeSite[i]]);}
    result.clusters=components.size();
    for(const auto& entry:components){const auto& endpoint=entry.second;require(endpoint.empty()||endpoint.size()==2,"Invalid meron line topology");if(endpoint.size()==2&&endpoint[0]!=endpoint[1])++result.merons;}
    return result;
}
double MeronModel::clusterTrace(const std::vector<MeronEvent>& events) const {
    auto c=clusters(events);if(c.merons)return 0.;
    return std::ldexp(1.,c.clusters-c.spinVertices+2*c.untouched);
}
DataType MeronModel::fermionPfaffian(double beta,const std::vector<MeronEvent>& events) const {
    validate();require(validEvents(beta,events),"Invalid meron event chronology");
    auto sigma=sublattices();auto h=numberMajorana(hopMatrix(*this));std::vector<MatType> pulses;DataType scalar=1.;double last=0.;
    auto density=[&](int i){MatType d=MatType::Zero(sites,sites);d(i,i)=1.;auto a=numberMajorana(d);pulses.push_back(DataType(0,pi/2)*a);pulses.push_back(DataType(0,pi/2)*a);scalar*=DataType(0,-.5);};
    for(const auto& e:events){pulses.push_back((e.time-last)*h);last=e.time;if(e.kind==MeronKind::DensityBond){auto [i,j]=edges[e.index];density(i);density(j);scalar*=-1.;}else if(e.kind==MeronKind::Monomer){density(e.index);scalar*=sigma[e.index];}}
    pulses.push_back((beta-last)*h);return scalar*traceProduct(pulses);
}
DataType MeronModel::groupedWeight(double beta,const std::vector<MeronEvent>& events) const {
    return vertexScalar(*this,events)*clusterTrace(events)*fermionPfaffian(beta,events);
}
DataType MeronModel::referenceWeight(double beta,const std::vector<MeronEvent>& events) const {
    validate();require(validEvents(beta,events),"Invalid meron reference events");require(sites<=6,"Meron spin Fock oracle limited to six active spins");
    int dim=1<<sites;MatType product=MatType::Identity(dim,dim);std::vector<bool> touched(sites,false);
    for(const auto& e:events) {
        if(e.kind==MeronKind::DensityBond)continue;MatType v=MatType::Zero(dim,dim);
        if(e.kind==MeronKind::Monomer){touched[e.index]=true;for(int s=0;s<dim;++s){v(s,s)=.5;v(s^(1<<e.index),s)=.5;}}
        else {auto [i,j]=edges[e.index];touched[i]=touched[j]=true;for(int s=0;s<dim;++s)if(((s>>i)&1)!=((s>>j)&1)){v(s,s)=.5;v(s^(1<<i)^(1<<j),s)=-.5;}}
        product=v*product;
    }
    int untouched=std::count(touched.begin(),touched.end(),false);
    return vertexScalar(*this,events)*std::ldexp(1.,untouched)*product.trace()*fermionReference(*this,beta,events);
}
bool MeronModel::validPath(double beta,const std::vector<MeronEvent>& events,const std::vector<std::uint64_t>& path) const {
    validate();if(!validEvents(beta,events)||path.size()!=events.size()+1||path.front()!=path.back())return false;
    for(auto v:path)if(v>=(std::uint64_t(1)<<(2*sites)))return false;
    for(int k=0;k<int(events.size());++k){bool found=false;for(auto [after,amplitude]:transitions(*this,events[k],path[k]))if(after==path[k+1])found=true;if(!found)return false;}return true;
}
DataType MeronModel::resolvedPfaffian(double beta,const std::vector<MeronEvent>& events,const std::vector<std::uint64_t>& path) const {
    if(!validPath(beta,events,path))return 0.;int total=3*sites;auto sigma=sublattices();auto a=numberMajorana(hopMatrix(*this));
    std::vector<Gaussian> factors{embed(projector(sites,path.front()),sites,total)};DataType scalar=1.;double last=0.;
    auto density=[&](int i){MatType d=MatType::Zero(sites,sites);d(i,i)=1.;auto g=gaussian(DataType(0,pi/2)*numberMajorana(d));factors.push_back(embed(g,0,total));factors.push_back(embed(g,0,total));scalar*=DataType(0,-.5);};
    for(int k=0;k<int(events.size());++k) {
        const auto& e=events[k];factors.push_back(embed(gaussian((e.time-last)*a),0,total));last=e.time;
        if(e.kind==MeronKind::DensityBond){auto [i,j]=edges[e.index];density(i);density(j);scalar*=-repulsion;continue;}
        if(e.kind==MeronKind::Monomer){density(e.index);scalar*=coupling*sigma[e.index];}else scalar*=exchange;
        for(auto [after,amplitude]:transitions(*this,e,path[k]))if(after==path[k+1])scalar*=amplitude;
        factors.push_back(embed(projector(sites,path[k]),sites,total));
        for(int i=0;i<sites;++i)if(((path[k]>>(2*i))&3u)!=((path[k+1]>>(2*i))&3u))factors.push_back(embed(rotation(sites,i,((path[k]>>(2*i))&3u)==0),sites,total));
    }
    factors.push_back(embed(gaussian((beta-last)*a),0,total));return scalar*traceGaussianFactors(factors);
}
DataType MeronModel::resolvedReference(double beta,const std::vector<MeronEvent>& events,const std::vector<std::uint64_t>& path) const {
    if(!validPath(beta,events,path))return 0.;double spin=1.;for(int k=0;k<int(events.size());++k)for(auto [after,amplitude]:transitions(*this,events[k],path[k]))if(after==path[k+1])spin*=amplitude;
    return vertexScalar(*this,events)*spin*fermionReference(*this,beta,events);
}
DataType MeronModel::exhaustivePfaffian(double beta,const std::vector<MeronEvent>& events,int* negativePaths) const {
    require(sites<=2&&events.size()<=5,"Exhaustive resolved Pfaffian oracle limited to two sites/five events");require(validEvents(beta,events),"Invalid exhaustive events");
    if(negativePaths)*negativePaths=0;DataType sum=0.;std::vector<std::uint64_t> path;
    std::function<void(int)> visit=[&](int k){if(k==int(events.size())){if(path.back()==path.front()){auto z=resolvedPfaffian(beta,events,path);sum+=z;if(negativePaths&&z.real()<-1e-12)++*negativePaths;}return;}for(auto [after,amplitude]:transitions(*this,events[k],path.back())){path.push_back(after);visit(k+1);path.pop_back();}};
    for(std::uint64_t first=0;first<(std::uint64_t(1)<<(2*sites));++first){path={first};visit(0);}return sum;
}
MatType MeronModel::exactHamiltonian() const {
    validate();auto p=physical(*this);MatType h=p.h0;for(const auto* group:{&p.density,&p.singlet,&p.monomer})for(const auto& v:*group)h+=v;return h;
}
DataType MeronModel::exactWordTrace(double beta,const std::vector<MeronEvent>& events) const {
    require(validEvents(beta,events),"Invalid physical meron events");auto p=physical(*this);MatType product=MatType::Identity(p.h0.rows(),p.h0.cols());double last=0.;
    for(const auto& e:events){product=matrixExp(-(e.time-last)*p.h0)*product;last=e.time;const auto& v=e.kind==MeronKind::DensityBond?p.density[e.index]:e.kind==MeronKind::SingletBond?p.singlet[e.index]:p.monomer[e.index];product=-v*product;}
    return (matrixExp(-(beta-last)*p.h0)*product).trace();
}
}
