#include "catalog/color_projector.h"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <queue>
#include <random>
using namespace catalog;

struct Records {
    std::ostream& out;int rows=0,positive=0,negative=0,zero=0;double error=0.,phase=0.;
    explicit Records(std::ostream& s):out(s){out<<std::setprecision(17)<<"case,index,modes,weight_real,weight_imag,reference_real,reference_imag,expected,relative_error,J,K,charging,chemical,density,active_density,beta,states,events\n";}
    void add(const std::string& name,int index,const ColorProjector& m,const ColorPath& p,DataType z,DataType ref,int expected=1) {
        double e=std::abs(z-ref)/std::max(1.,std::abs(ref));require(std::isfinite(e) && e<2e-8,"Color signed Pfaffian/reference disagreement");error=std::max(error,e);
        if(!expected){require(std::abs(z)<1e-11,"Expected zero color path");++zero;}
        else {require(std::abs(z)>1e-14,"Unresolved color path weight");double d=std::abs(z/std::abs(z)-double(expected));require(d<2e-8,"Wrong color signed phase");phase=std::max(phase,d);if(expected>0)++positive;else ++negative;}
        out<<name<<','<<index<<','<<m.modes()<<','<<z.real()<<','<<z.imag()<<','<<ref.real()<<','<<ref.imag()<<','<<expected<<','<<e<<','<<m.singlet<<','<<m.permutation<<','<<m.charging<<','<<m.chemical<<','<<m.density<<','<<m.activeDensity<<','<<p.beta<<',';
        for(auto s:p.states)out<<s<<':';out<<',';for(auto e0:p.events)out<<e0.time<<':'<<e0.edge<<':'<<(e0.kind==ColorVertex::Singlet?'P':'S')<<';';out<<'\n';++rows;
    }
};
static int phase(std::uint64_t s,int sites){int p=1;for(int r=0;r<sites;++r)p*=ColorProjector::occupationPhase((s>>(3*r))&7);return p;}
static void setTimes(ColorPath& p) {for(int k=0;k<int(p.events.size());++k)p.events[k].time=p.beta*(k+1)/(p.events.size()+1.);}
static DataType physicalWord(const ColorProjector& m,const ColorPath& p,const std::vector<MatType>& v) {
    DataType z=1.;double last=0.;
    for(int k=0;k<int(p.events.size());++k){auto e=p.events[k];double c=e.kind==ColorVertex::Singlet?m.singlet:m.permutation;z*=std::exp(-(e.time-last)*m.diagonalEnergy(p.states[k]))*c*v[2*e.edge+(e.kind==ColorVertex::Permutation)](p.states[k+1],p.states[k]);last=e.time;}
    return z*std::exp(-(p.beta-last)*m.diagonalEnergy(p.states.back()));
}
static std::vector<MatType> physicalVertices(const ColorProjector& m) {
    std::vector<MatType> v;for(int e=0;e<int(m.edges.size());++e)for(auto k:{ColorVertex::Singlet,ColorVertex::Permutation})v.push_back(m.vertexMatrix(e,k));return v;
}
static void localControls() {
    auto a=ColorProjector::localActive();MatType casimir=MatType::Zero(8,8);
    for(int k=0;k<3;++k){auto s=ColorProjector::localSpin(k);require((s-s.adjoint()).norm()<1e-13,"SO3 spin not Hermitian");casimir+=s*s;}
    require((casimir-2.*a).norm()<1e-13,"Full Fock spin0+spin1 decomposition failed");
    auto p=ColorProjector::localBondVertex(ColorVertex::Singlet),swap=ColorProjector::localBondVertex(ColorVertex::Permutation);
    require((p*p-p).norm()<1e-13 && (p*swap-p).norm()<1e-13,"Color projector/permutation algebra failed");
    MatType square=swap*swap;for(int s=0;s<64;++s)for(int t=0;t<64;++t)require(std::abs(square(t,s)-(t==s?a(s&7,s&7)*a(s>>3,s>>3):DataType(0.)))<1e-13,"Color swap square/blocked space failed");
    // This full-Fock coefficient is genuinely eight-fermion order.
    auto coeff=[&](const MatType& v){DataType sum=0.;for(unsigned bits=0;bits<16;++bits){unsigned state=(bits&3)|((bits>>2)<<3);sum+=(__builtin_popcount(bits)%2?-1.:1.)*v(state,state);}return sum;};
    require(std::abs(coeff(p)-1.)<1e-13 && std::abs(coeff(swap)-3.)<1e-13,"Local four-density coefficient disappeared");
    // R008's positive Gaussian vertex is a positive product across its two cells.
    // The N_i=N_j=1 block here is a rank-one entangled, non-Slater projector.
    // This excludes its unshifted projector from that positive product cone;
    // it does NOT exclude arbitrary constant shifts or arbitrary new expansions.
    MatType block(9,9),pt(9,9);
    for(int x=0;x<3;++x)for(int y=0;y<3;++y)for(int u=0;u<3;++u)for(int v=0;v<3;++v)block(3*x+y,3*u+v)=p((1<<x)|((1<<y)<<3),(1<<u)|((1<<v)<<3));
    for(int x=0;x<3;++x)for(int y=0;y<3;++y)for(int u=0;u<3;++u)for(int v=0;v<3;++v)pt(3*x+y,3*u+v)=block(3*x+v,3*u+y);
    Eigen::SelfAdjointEigenSolver<MatType> spectrum(block),partialTranspose(pt);
    require((block*block-block).norm()<1e-13 && std::abs(block.trace()-1.)<1e-13 && std::abs(partialTranspose.eigenvalues()(0)+1./3.)<1e-13,"Entangled singlet block certificate failed");
    MatType pair=MatType::Zero(6,6);for(int c=0;c<3;++c){pair(c,3+c)=1./std::sqrt(3.);pair(3+c,c)=-1./std::sqrt(3.);}
    require(Eigen::FullPivLU<MatType>(pair).rank()==6,"Singlet unexpectedly has Slater rank one");
    std::cerr<<"PASS local controls: eight-fermion coefficient(P,swap)=(1,3), singlet partial-transpose min=-1/3, pair rank=6\n";
}
static void controls(Records& r) {
    localControls();ColorProjector m;m.sites=2;m.edges={{0,1}};
    auto v=physicalVertices(m);
    // All four particle/hole spin1 combinations, with an odd-length nonpalindrome.
    for(int ni=1;ni<=2;++ni)for(int nj=1;nj<=2;++nj){ColorPath p;
        for(int c:{0,1,2,0})p.states.push_back(ColorProjector::occupation(ni,c)|(std::uint64_t(ColorProjector::occupation(nj,c))<<3));
        p.events={{0.,0,ColorVertex::Singlet},{0.,0,ColorVertex::Singlet},{0.,0,ColorVertex::Singlet}};setTimes(p);
        auto ref=physicalWord(m,p,v);r.add("nonpalindrome_three_color_cycle",2*(ni-1)+nj,m,p,m.weight(p),ref);
        if(ni==1 && nj==2){ColorProjector bad=m;bad.singlet=-m.singlet;require(!bad.signFreeParameters(),"Negative J incorrectly admitted");r.add("negative_singlet_three_color_cycle",0,bad,p,bad.weight(p),physicalWord(bad,p,v),-1);}
    }
    ColorPath blocked;blocked.states={0,0};blocked.events={{.3,0,ColorVertex::Singlet}};r.add("vacuum_blocker",0,m,blocked,m.weight(blocked),0.,0);
    blocked.states={63,63};r.add("filled_blocker",0,m,blocked,m.weight(blocked),0.,0);
    ColorPath invalid;invalid.states={9,18};invalid.events={{.3,0,ColorVertex::Singlet}};r.add("nonclosed",0,m,invalid,m.weight(invalid),0.,0);
    ColorProjector z=m;z.singlet=0.;ColorPath zero;zero.states={9,18,9};zero.events={{.2,0,ColorVertex::Singlet},{.4,0,ColorVertex::Singlet}};r.add("zero_coupling",0,z,zero,z.weight(zero),0.,0);
    // Every physical transition, including mixed N=1/N=2 blocks, is checked
    // against independent canonical-Fock spin polynomials and the local gauge.
    for(int kind=0;kind<2;++kind)for(unsigned s=0;s<64;++s){auto k=kind?ColorVertex::Permutation:ColorVertex::Singlet;auto next=m.nextStates(s,0,k);
        for(unsigned t=0;t<64;++t){bool allowed=std::find(next.begin(),next.end(),t)!=next.end();DataType expected=allowed?double(phase(s,2)*phase(t,2))*(kind?1.:1./3.):0.;require(std::abs(v[kind](t,s)-expected)<1e-13,"Independent full Fock vertex/gauge certificate failed");}}
    // Complete 3-vertex SSE coefficient. Includes every diagonal vertex and
    // both singlet/permutation choices; no path is chosen by its sign.
    DataType sum=0.;int index=0;
    for(unsigned s=0;s<64;++s)for(auto k0:{ColorVertex::Singlet,ColorVertex::Permutation})for(auto t:m.nextStates(s,0,k0))for(auto k1:{ColorVertex::Singlet,ColorVertex::Permutation})for(auto u:m.nextStates(t,0,k1))for(auto k2:{ColorVertex::Singlet,ColorVertex::Permutation}) {
        auto next=m.nextStates(u,0,k2);if(std::find(next.begin(),next.end(),s)==next.end())continue;
        ColorPath p;p.states={s,t,u,s};p.events={{.2,0,k0},{.4,0,k1},{.6,0,k2}};auto w=m.weight(p);r.add("all_three_vertex_words",index++,m,p,w,physicalWord(m,p,v));sum+=w;
    }
    MatType vertex=m.singlet*v[0]+m.permutation*v[1],thermal=MatType::Zero(64,64);for(int s=0;s<64;++s)thermal(s,s)=std::exp(-.8*m.diagonalEnergy(s));
    require(std::abs(sum-(thermal*vertex*vertex*vertex).trace())<2e-10,"Complete SSE order-three trace failed");
    auto h=m.exactHamiltonian();auto thermalExact=matrixExp(-.8*h),interaction=matrixExp(.8*vertex);
    require((thermalExact-thermal*interaction).norm()<2e-10,"Independent physical Hamiltonian exponential failed");
    // The genuine negative biquadratic model is included without dropping the
    // active-sector-dependent constant in the full physical Fock space.
    ColorProjector pure=m;pure.permutation=pure.charging=pure.chemical=pure.density=0.;pure.activeDensity=-pure.singlet/3.;
    MatType dot=MatType::Zero(64,64);for(int axis=0;axis<3;++axis){auto s=ColorProjector::localSpin(axis);for(int u=0;u<64;++u)for(int t=0;t<64;++t)dot(t,u)+=s(t&7,u&7)*s(t>>3,u>>3);}
    require((pure.exactHamiltonian()+pure.singlet/3.*dot*dot).norm()<1e-12,"Pure negative biquadratic embedding dropped active constant");
    ColorPath bp;bp.states={9,18,36,9};bp.events={{.2,0,ColorVertex::Singlet},{.4,0,ColorVertex::Singlet},{.6,0,ColorVertex::Singlet}};r.add("pure_negative_biquadratic",0,pure,bp,pure.weight(bp),physicalWord(pure,bp,v));
    std::cerr<<"PASS full Fock controls: all 8192 local matrix entries and "<<index<<" resolved SSE words\n";
}
struct Step {std::uint64_t to;int edge;ColorVertex kind;};
static std::vector<Step> steps(const ColorProjector& m,std::uint64_t s) {
    std::vector<Step> out;for(int e=0;e<int(m.edges.size());++e)for(auto k:{ColorVertex::Singlet,ColorVertex::Permutation})for(auto t:m.nextStates(s,e,k))out.push_back({t,e,k});return out;
}
static ColorPath randomClosed(const ColorProjector& m,std::mt19937_64& rng,int index) {
    ColorPath p;p.beta=.6+.03*(index%7);std::uint64_t start=0;
    for(int site=0;site<m.sites;++site){unsigned s=index%4?ColorProjector::occupation(1+rng()%2,rng()%3):rng()%8;start|=std::uint64_t(s)<<(3*site);}p.states={start};
    for(int k=0;k<3+index%4;++k){auto next=steps(m,p.states.back());if(next.empty())break;auto a=next[rng()%next.size()];p.events.push_back({0.,a.edge,a.kind});p.states.push_back(a.to);}
    // Close by a shortest path in the configuration graph, not a forced reverse.
    std::uint64_t end=p.states.back();std::queue<std::uint64_t> q;q.push(end);std::map<std::uint64_t,std::pair<std::uint64_t,Step>> prev;prev[end]={end,{end,0,ColorVertex::Singlet}};
    while(!q.empty() && !prev.count(start)){auto s=q.front();q.pop();for(auto a:steps(m,s))if(!prev.count(a.to)){prev[a.to]={s,a};q.push(a.to);}}
    require(prev.count(start),"Random color walk cannot return");std::vector<Step> close;
    for(auto s=start;s!=end;s=prev.at(s).first)close.push_back(prev.at(s).second);
    std::reverse(close.begin(),close.end());for(auto a:close){p.events.push_back({0.,a.edge,a.kind});p.states.push_back(a.to);}setTimes(p);return p;
}
int main(int argc,char** argv) {
    try {
        mkl_set_num_threads(1);mkl_set_dynamic(0);int samples=32;unsigned seed=240602926;std::string csv;
        for(int i=1;i<argc;++i){std::string key=argv[i];require(i+1<argc,"Missing color CLI value");std::string value=argv[++i];if(key=="--samples")samples=std::stoi(value);else if(key=="--seed")seed=std::stoul(value);else if(key=="--csv")csv=value;else throw std::runtime_error("Unknown color option");}
        require(samples>0,"Positive color samples required");std::ofstream file;if(!csv.empty()){file.open(csv);require(file.good(),"Cannot open color CSV");}Records r(csv.empty()?std::cout:file);controls(r);std::mt19937_64 rng(seed);
        for(int sites:{2,3}){ColorProjector m;if(sites==2){m.sites=2;m.edges={{0,1}};}require(m.signFreeParameters(),"Default color parameters not positive");auto v=physicalVertices(m);auto h=m.exactHamiltonian();require((h-h.adjoint()).norm()<1e-11,"Physical color Hamiltonian not Hermitian");
            int moved=0;for(int index=0;index<samples;++index){auto p=randomClosed(m,rng,index);for(int k=0;k<int(p.events.size());++k)if(p.states[k]!=p.states[k+1]){++moved;break;}
                auto z=m.weight(p),ref=physicalWord(m,p,v);require(std::abs(ref-m.reference(p))<1e-10,"Color independent path oracle disagreement");r.add(sites==3?"triangle_random_paths":"two_cell_random_paths",index,m,p,z,ref);}
            require(moved>0,"Random color ensemble never changed a color");std::cerr<<"PASS color "<<sites<<" cells moved="<<moved<<'\n';
        }
        std::cerr<<"PASS color_projector rows="<<r.rows<<" positive="<<r.positive<<" zero="<<r.zero<<" negative="<<r.negative<<" max_relative_error="<<r.error<<" max_phase_error="<<r.phase<<" seed="<<seed<<'\n';
    }catch(const std::exception& e){std::cerr<<"FAIL color_projector: "<<e.what()<<'\n';return 1;}return 0;
}
