#include "catalog/jq_models.h"
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
    explicit Records(std::ostream& s):out(s){out<<std::setprecision(17)<<"case,index,modes,weight_real,weight_imag,reference_real,reference_imag,expected,relative_error,beta,terms,fields,states,events\n";}
    void add(const std::string& name,int index,const JQModel& m,const JQPath& p,DataType z,DataType ref,int expected=1) {
        double e=std::abs(z-ref)/std::max(1.,std::abs(ref));require(std::isfinite(e) && e<3e-8,"JQ signed Pfaffian/reference mismatch");error=std::max(error,e);
        if(!expected){require(std::abs(z)<1e-11,"Expected zero JQ path");++zero;}else{require(std::abs(z)>1e-15,"Unresolved JQ weight");double d=std::abs(z/std::abs(z)-double(expected));require(d<3e-8,"Wrong JQ signed phase");phase=std::max(phase,d);if(expected>0)++positive;else ++negative;}
        out<<name<<','<<index<<','<<m.modes()<<','<<z.real()<<','<<z.imag()<<','<<ref.real()<<','<<ref.imag()<<','<<expected<<','<<e<<','<<p.beta<<',';
        for(const auto& t:m.terms){out<<t.coupling<<':';for(auto [i,j]:t.bonds)out<<i<<'-'<<j<<'/';out<<';';}out<<',';for(auto h:m.fields)out<<h<<':';out<<',';for(auto s:p.states)out<<s<<':';out<<',';for(auto e0:p.events)out<<e0.time<<':'<<e0.term<<';';out<<'\n';++rows;
    }
};
static void times(JQPath& p){for(int k=0;k<int(p.events.size());++k)p.events[k].time=p.beta*(k+1)/(p.events.size()+1.);}
static std::vector<MatType> vertices(const JQModel& m,bool full) {std::vector<MatType> v;for(int t=0;t<int(m.terms.size());++t)v.push_back(full?m.vertexMatrix(t):m.spinVertexMatrix(t));return v;}
static DataType physicalWord(const JQModel& m,const JQPath& p,const std::vector<MatType>& v,bool full) {
    DataType z=1.;double last=0.;for(int k=0;k<int(p.events.size());++k){auto e=p.events[k];auto s=full?p.states[k]:m.spinState(p.states[k]),t=full?p.states[k+1]:m.spinState(p.states[k+1]);z*=m.terms[e.term].coupling*v[e.term](t,s)*std::exp(-(e.time-last)*m.diagonalEnergy(p.states[k]));last=e.time;}return z*std::exp(-(p.beta-last)*m.diagonalEnergy(p.states.back()));
}
static void controls(Records& r) {
    MatType casimir=MatType::Zero(4,4);for(int a=0;a<3;++a){auto s=JQModel::localSpin(a);require((s-s.adjoint()).norm()<1e-13,"JQ physical spin not Hermitian");casimir+=s*s;}
    require((casimir-.75*MatType::Identity(4,4)).norm()<1e-13,"JQ has unwanted spin-zero blockers");
    MatType plus=JQModel::localSpin(0)+1.i*JQModel::localSpin(1);for(unsigned p=0;p<2;++p)require(std::abs(plus(JQModel::occupation(0,p,1),JQModel::occupation(1,p,1))-1.)<1e-13,"JQ local parity-copy basis incorrect");
    JQModel m;auto f=vertices(m,true),s=vertices(m,false);auto sub=m.sublattices();require(!sub.empty() && m.signFreeParameters(),"JQ square must be sign free");
    // Independent full physical Fock vertex matrices versus Pauli spin matrices,
    // all local gauge-copy sectors and all allowed/forbidden transitions.
    int checked=0;
    for(int k=0;k<int(f.size());++k)for(unsigned before=0;before<256;++before){auto next=m.nextStates(before,k);
        for(unsigned after=0;after<256;++after){bool paritySame=true;for(int site=0;site<4;++site)if((__builtin_popcount((before>>(2*site))&3)%2)!=(__builtin_popcount((after>>(2*site))&3)%2))paritySame=false;
            DataType expected=paritySame?s[k](m.spinState(after),m.spinState(before)):0.;require(std::abs(f[k](after,before)-expected)<1e-13,"JQ Fock/Pauli copy identity failed");
            bool allowed=std::find(next.begin(),next.end(),after)!=next.end();require((std::abs(expected)>1e-13)==allowed,"JQ transition support mismatch");
            if(allowed){int flips=0;for(int site=0;site<4;++site)if(sub[site]>0 && ((before^after)>>(2*site)&1))++flips;require(std::abs(expected-(flips%2?-1.:1.)*std::ldexp(1.,-int(m.terms[k].bonds.size())))<1e-13,"JQ Marshall gauge certificate failed");}++checked;
        }}
    auto q=s[4],other=s[5];require((q*q-q).norm()<1e-13,"JQ product is not a projector");
    require(std::abs((q*other-other*q).norm()-std::sqrt(3./8.))<1e-13,"JQ noncommuting-bond certificate failed");
    // A conditional two-spin swap distinguishes Q2 from ordinary pair hopping
    // plus a pure four-spin XY ring. The unchanged remote bond is essential.
    require(std::abs(q(6,5)+.25)<1e-13 && q(2,1)==DataType(0.) && std::abs(s[0](6,5)+.5)<1e-13 && std::abs(s[0](2,1)+.5)<1e-13,"JQ conditional-swap certificate failed");
    require(std::abs(q(10,5)-.25)<1e-13,"JQ four-site matrix element disappeared");
    DataType density4=0.;for(unsigned bits=0;bits<16;++bits){unsigned physical=0;for(int site=0;site<4;++site)physical|=((bits>>site)&1)<<(2*site);density4+=(__builtin_popcount(bits)%2?-1.:1.)*f[4](physical,physical);}
    require(std::abs(density4-1.)<1e-13,"JQ eight-fermion four-density coefficient failed");
    // Odd-length, nonpalindromic Q2 path in EVERY parity-copy sector.
    for(unsigned parity=0;parity<16;++parity){JQPath p;for(unsigned state:{5u,6u,10u,5u})p.states.push_back(JQModel::occupation(state,parity,4));p.events={{.2,4},{.4,4},{.6,4}};
        r.add("nonpalindrome_Q2_all_parities",parity,m,p,m.weight(p),physicalWord(m,p,f,true));
        if(parity==3){JQModel bad=m;bad.terms[4].coupling=-bad.terms[4].coupling;require(!bad.signFreeParameters(),"Negative Q admitted");r.add("negative_Q_odd_word",0,bad,p,bad.weight(p),physicalWord(bad,p,f,true),-1);}}
    // Positive couplings on a triangle: the physical spin exchange gives a
    // genuine negative signed Gaussian trace, testing bipartiteness itself.
    JQModel triangle;triangle.sites=3;triangle.fields.assign(3,0.);triangle.terms={{{{0,1}},1.},{{{1,2}},1.},{{{2,0}},1.}};require(!triangle.signFreeParameters(),"Triangle incorrectly admitted");
    JQPath negative;negative.states={3,12,48,3};negative.events={{.2,0},{.4,1},{.6,2}};auto tv=vertices(triangle,true);r.add("triangle_signed_exchange",0,triangle,negative,triangle.weight(negative),physicalWord(triangle,negative,tv,true),-1);
    JQPath blocked;blocked.states={0,0};blocked.events={{.3,4}};r.add("parallel_spin_annihilated",0,m,blocked,m.weight(blocked),0.,0);
    JQPath invalid;invalid.states={51,60};invalid.events={{.3,4}};r.add("nonclosed",0,m,invalid,m.weight(invalid),0.,0);
    JQModel zero=m;zero.terms[4].coupling=0.;JQPath zp;zp.states={51,60,51};zp.events={{.2,4},{.4,4}};r.add("zero_Q",0,zero,zp,zero.weight(zp),0.,0);
    // Complete order-two coefficient in the entire 256-dimensional Fock space.
    DataType sum=0.;int words=0;for(unsigned before=0;before<256;++before)for(int a=0;a<int(m.terms.size());++a)for(auto middle:m.nextStates(before,a))for(int b=0;b<int(m.terms.size());++b){auto next=m.nextStates(middle,b);if(std::find(next.begin(),next.end(),before)==next.end())continue;
        JQPath p;p.states={before,middle,before};p.events={{.2,a},{.6,b}};auto z=m.weight(p);r.add("all_order_two_Fock_words",words++,m,p,z,physicalWord(m,p,f,true));sum+=z;}
    MatType h=m.exactHamiltonian(),hs=m.exactSpinHamiltonian();require(std::abs(sum-(h*h).trace())<3e-9,"JQ complete SSE trace coefficient failed");
    Eigen::SelfAdjointEigenSolver<MatType> ef(h),es(hs);for(int i=0;i<256;++i)require(std::abs(ef.eigenvalues()(i)-es.eigenvalues()(i/16))<2e-11,"JQ full Fock spectrum not 16 spin copies");
    require(std::abs(matrixExp(-.4*h).trace()-16.*matrixExp(-.4*hs).trace())<3e-9,"JQ full partition trace lost copy multiplicity");
    std::cerr<<"PASS JQ controls matrix_entries="<<checked<<" complete_words="<<words<<" Q_commutator=sqrt(3/8) eight_fermion_coefficient=1 full_Fock_copies=16\n";
}
struct Step {std::uint64_t to;int term;};
static std::vector<Step> steps(const JQModel& m,std::uint64_t s){std::vector<Step> out;for(int t=0;t<int(m.terms.size());++t)for(auto u:m.nextStates(s,t))out.push_back({u,t});return out;}
static JQPath randomClosed(const JQModel& m,std::mt19937_64& rng,int index) {
    JQPath p;p.beta=.6+.03*(index%9);auto mask=(std::uint64_t(1)<<m.sites)-1;auto start=JQModel::occupation(rng()&mask,rng()&mask,m.sites);p.states={start};
    for(int k=0;k<3+index%4;++k){auto next=steps(m,p.states.back());if(next.empty())break;auto a=next[rng()%next.size()];p.events.push_back({0.,a.term});p.states.push_back(a.to);}
    auto end=p.states.back();std::queue<std::uint64_t> q;q.push(end);std::map<std::uint64_t,std::pair<std::uint64_t,Step>> prev;prev[end]={end,{end,0}};
    while(!q.empty() && !prev.count(start)){auto s=q.front();q.pop();for(auto a:steps(m,s))if(!prev.count(a.to)){prev[a.to]={s,a};q.push(a.to);}}
    require(prev.count(start),"JQ random path did not close");std::vector<Step> close;for(auto s=start;s!=end;s=prev.at(s).first)close.push_back(prev.at(s).second);std::reverse(close.begin(),close.end());
    for(auto a:close){p.events.push_back({0.,a.term});p.states.push_back(a.to);}times(p);return p;
}
int main(int argc,char** argv) {
    try {
        mkl_set_num_threads(1);mkl_set_dynamic(0);int samples=32;unsigned seed=611343;std::string csv;
        for(int i=1;i<argc;++i){std::string key=argv[i];require(i+1<argc,"Missing JQ CLI argument");std::string value=argv[++i];if(key=="--samples")samples=std::stoi(value);else if(key=="--seed")seed=std::stoul(value);else if(key=="--csv")csv=value;else throw std::runtime_error("Unknown JQ option");}
        require(samples>0,"Positive JQ samples required");std::ofstream file;if(!csv.empty()){file.open(csv);require(file.good(),"Cannot open JQ CSV");}Records r(csv.empty()?std::cout:file);controls(r);std::mt19937_64 rng(seed);
        for(int order:{2,3}){JQModel m=JQModel::rectangle(order,2,order);for(int site=0;site<m.sites;++site)m.fields[site]=.03*(site-1);require(m.signFreeParameters(),"JQ rectangular builder not sign free");auto sv=vertices(m,false);
            int qterm=-1;for(int t=0;t<int(m.terms.size());++t)if(int(m.terms[t].bonds.size())==order)qterm=t;require(qterm>=0,"JQ builder omitted multi-bond term");
            // Explicit simultaneous multi-bond flips exercise Q3, not merely J.
            std::uint64_t spin=0;for(auto [i,j]:m.terms[qterm].bonds)spin|=std::uint64_t(1)<<i;
            auto initial=JQModel::occupation(spin,5,m.sites);auto next=m.nextStates(initial,qterm);JQPath p;p.states={initial,next.back(),initial};p.events={{.2,qterm},{.6,qterm}};
            r.add(order==2?"explicit_Q2":"explicit_Q3",0,m,p,m.weight(p),physicalWord(m,p,sv,false));
            int moved=0;for(int index=0;index<samples;++index){auto path=randomClosed(m,rng,index);for(int k=0;k<int(path.events.size());++k)if(path.states[k]!=path.states[k+1]){++moved;break;}
                auto ref=physicalWord(m,path,sv,false);require(std::abs(ref-m.reference(path))<1e-11,"JQ independent Pauli word disagreement");r.add(order==2?"Q2_random":"Q3_random",index,m,path,m.weight(path),ref);}
            require(moved>0,"JQ random ensemble never exchanged spins");std::cerr<<"PASS JQ order="<<order<<" sites="<<m.sites<<" moved="<<moved<<'\n';
        }
        std::cerr<<"PASS jq rows="<<r.rows<<" positive="<<r.positive<<" zero="<<r.zero<<" negative="<<r.negative<<" max_relative_error="<<r.error<<" max_phase_error="<<r.phase<<" seed="<<seed<<'\n';
    }catch(const std::exception& e){std::cerr<<"FAIL jq: "<<e.what()<<'\n';return 1;}return 0;
}
