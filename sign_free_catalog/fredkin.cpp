#include "catalog/fredkin_models.h"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <queue>
#include <random>
using namespace catalog;
struct Records {
    std::ostream& out;int rows=0,positive=0,negative=0,zero=0;double error=0.,phase=0.;
    explicit Records(std::ostream& s):out(s){out<<std::setprecision(17)<<"case,index,modes,weight_real,weight_imag,reference_real,reference_imag,expected,relative_error,beta,boundary_left,boundary_right,terms,fields,states,events\n";}
    void add(const std::string& name,int index,const FredkinModel& m,const FredkinPath& p,DataType z,DataType ref,int expected=1) {
        double e=std::abs(z-ref)/std::max(1.,std::abs(ref));require(std::isfinite(e) && e<3e-8,"Fredkin signed Pfaffian/reference mismatch");error=std::max(error,e);
        if(!expected){require(std::abs(z)<1e-11,"Expected zero Fredkin path");++zero;}else{require(std::abs(z)>1e-15,"Unresolved Fredkin weight");double d=std::abs(z/std::abs(z)-double(expected));require(d<3e-8,"Wrong Fredkin signed phase");phase=std::max(phase,d);if(expected>0)++positive;else ++negative;}
        out<<name<<','<<index<<','<<m.modes()<<','<<z.real()<<','<<z.imag()<<','<<ref.real()<<','<<ref.imag()<<','<<expected<<','<<e<<','<<p.beta<<','<<m.boundaryLeft<<','<<m.boundaryRight<<',';
        for(auto t:m.terms)out<<t.control<<':'<<t.first<<':'<<t.second<<':'<<t.up<<':'<<t.coupling<<';';out<<',';for(auto h:m.fields)out<<h<<':';out<<',';for(auto s:p.states)out<<s<<':';out<<',';for(auto e0:p.events)out<<e0.time<<':'<<e0.term<<';';out<<'\n';++rows;
    }
};
struct Oracle {MatType h;std::vector<MatType> terms;bool full;};
static Oracle oracle(const FredkinModel& m,bool full) {Oracle o;o.full=full;o.h=m.exactHamiltonian(full);for(int t=0;t<int(m.terms.size());++t)o.terms.push_back(m.termMatrix(t,full));return o;}
static DataType physicalWord(const FredkinModel& m,const FredkinPath& p,const Oracle& o) {
    DataType w=1.;double last=0.;for(int k=0;k<int(p.events.size());++k){auto e=p.events[k];auto s=o.full?p.states[k]:m.spinState(p.states[k]),t=o.full?p.states[k+1]:m.spinState(p.states[k+1]);
        DataType c=e.term<0?-o.h(t,s):-m.terms[e.term].coupling*o.terms[e.term](t,s);w*=c*std::exp(-(e.time-last)*o.h(s,s).real());last=e.time;}
    auto s=o.full?p.states.back():m.spinState(p.states.back());return w*std::exp(-(p.beta-last)*o.h(s,s).real());
}
static void times(FredkinPath& p){for(int k=0;k<int(p.events.size());++k)p.events[k].time=p.beta*(k+1)/(p.events.size()+1.);}
static void controls(Records& r) {
    FredkinModel m;auto f=oracle(m,true),s=oracle(m,false);require(m.signFreeParameters(),"Default Fredkin chain outside positive domain");int entries=0;
    for(unsigned before=0;before<256;++before){require(std::abs(f.h(before,before)-m.diagonalEnergy(before))<1e-13,"Fredkin dropped diagonal part of conditional projectors");
        for(unsigned after=0;after<256;++after){bool same=true;for(int site=0;site<4;++site)if((__builtin_popcount((before>>(2*site))&3)%2)!=(__builtin_popcount((after>>(2*site))&3)%2))same=false;
            for(int t=0;t<int(m.terms.size());++t){DataType expected=same?s.terms[t](m.spinState(after),m.spinState(before)):0.;require(std::abs(f.terms[t](after,before)-expected)<1e-13,"Fredkin full-Fock parity copy identity failed");if(before!=after)require(std::abs(-m.terms[t].coupling*expected-m.transitionCoefficient(before,after,t))<1e-13,"Fredkin transition coefficient differs from physical operator");++entries;}
            if(before!=after){require(f.h(after,before).real()<=1e-13,"Positive Fredkin couplings have a positive off-diagonal H entry");require(std::abs(-f.h(after,before)-m.transitionCoefficient(before,after,-1))<1e-13,"Fredkin combined transition coefficient failed");}
        }}
    // Genuine six-fermion term, including the equal-coupling model. Its purely
    // diagonal cubic-density coefficient cancels between left and right terms,
    // so that would be the wrong certificate. Z0 X1 X2 has coefficient -1/8.
    auto three=FredkinModel::chain(3);auto threeH=three.exactHamiltonian(false);MatType word=MatType::Zero(8,8);
    for(unsigned state=0;state<8;++state)word(state^6,state)=(state&1)?-1.:1.;
    require(std::abs((threeH*word).trace()/8.+.125)<1e-13,"Fredkin irreducible three-spin coefficient failed");
    require(std::abs(threeH(4,2)+.5)<1e-13 && threeH(5,3)==DataType(0.),"Fredkin single-control conditional swap failed");
    require((s.terms[0]*s.terms[1]-s.terms[1]*s.terms[0]).norm()>.1,"Fredkin terms unexpectedly form commuting stabilizers");
    Eigen::SelfAdjointEigenSolver<MatType> ef(f.h),es(s.h);for(int i=0;i<256;++i)require(std::abs(ef.eigenvalues()(i)-es.eigenvalues()(i/16))<2e-11,"Fredkin Fock spectrum not 16 copies");
    require(es.eigenvalues()(0)>-1e-12 && es.eigenvalues()(1)>1e-5,"Fredkin physical spin ground-state multiplicity incorrect");cVecType dyck=cVecType::Zero(16);dyck(10)=dyck(12)=1./std::sqrt(2.);require((s.h*dyck).norm()<1e-13,"Normalized Fredkin Dyck state has wrong energy");
    require(std::abs(matrixExp(-.4*f.h).trace()-16.*matrixExp(-.4*s.h).trace())<3e-9,"Fredkin trace-copy multiplicity lost");
    // Every local parity sector participates; no spin-zero state is removed.
    for(unsigned parity=0;parity<16;++parity){FredkinPath p;p.states={JQModel::occupation(12,parity,4),JQModel::occupation(10,parity,4),JQModel::occupation(12,parity,4)};p.events={{.2,0},{.6,3}};r.add("all_parity_copies",parity,m,p,m.weight(p),physicalWord(m,p,f));}
    FredkinPath frozen;frozen.states={JQModel::occupation(15,7,4)};r.add("non_Dyck_state_retained",0,m,frozen,m.weight(frozen),physicalWord(m,frozen,f));
    FredkinPath blocked;blocked.states={0,15,0};blocked.events={{.2,0},{.6,0}};r.add("inactive_control_or_pair",0,m,blocked,m.weight(blocked),0.,0);
    FredkinPath invalid;invalid.states={51,60};invalid.events={{.3,0}};r.add("nonclosed",0,m,invalid,m.weight(invalid),0.,0);
    auto z=m;z.terms[0].coupling=0.;FredkinPath zp;zp.states={JQModel::occupation(12,0,4),JQModel::occupation(10,0,4),JQModel::occupation(12,0,4)};zp.events={{.2,0},{.6,0}};r.add("zero_coupling",0,z,zp,z.weight(zp),0.,0);
    // Exact complete order-two CT path sum at the stated times. D contains all
    // diagonal bulk and boundary terms, and parallel labels are all retained.
    DataType sum=0.;int count=0;for(unsigned before=0;before<256;++before)for(int a=0;a<int(m.terms.size());++a)for(auto middle:m.nextStates(before,a))for(int b=0;b<int(m.terms.size());++b){auto next=m.nextStates(middle,b);if(std::find(next.begin(),next.end(),before)==next.end())continue;
        FredkinPath p;p.states={before,middle,before};p.events={{.2,a},{.55,b}};auto w=m.weight(p);r.add("all_two_vertex_Fock_words",count++,m,p,w,physicalWord(m,p,f));sum+=w;}
    MatType transfer=-f.h,d0=MatType::Zero(256,256),d1=d0,d2=d0;for(int i=0;i<256;++i){double e=f.h(i,i).real();transfer(i,i)=0.;d0(i,i)=std::exp(-.2*e);d1(i,i)=std::exp(-.35*e);d2(i,i)=std::exp(-.25*e);}
    require(std::abs(sum-(d2*transfer*d1*transfer*d0).trace())<2e-10,"Fredkin complete CT word sum failed");
    // A four-step loop remains negative AFTER combining all parallel physical
    // hopping terms, so it cannot be repaired by merely relabeling vertices.
    auto bad=FredkinModel::chain(6);bad.terms[2].coupling=-2.;require(!bad.signFreeParameters(),"Negative Fredkin coupling admitted");auto bo=oracle(bad,false);FredkinPath bp;
    for(unsigned state:{50u,42u,44u,52u,50u})bp.states.push_back(JQModel::occupation(state,9,6));bp.events={{.1,-1},{.3,-1},{.5,-1},{.7,-1}};
    const double coeff[4]={-.5,1.,.5,.5};for(int k=0;k<4;++k)require(std::abs(bad.transitionCoefficient(bp.states[k],bp.states[k+1],-1)-coeff[k])<1e-13,"Fredkin combined negative loop changed");
    auto negative=bad.weight(bp);r.add("negative_combined_physical_cycle",0,bad,bp,negative,physicalWord(bad,bp,bo),-1);
    DataType resolved=0.;int words=0;auto resolvedPath=bp;
    std::function<void(int)> visit=[&](int at){if(at==4){auto ref=physicalWord(bad,resolvedPath,bo),w=bad.weight(resolvedPath);r.add("negative_domain_resolved_word",words++,bad,resolvedPath,w,ref,ref.real()>0?1:-1);resolved+=w;return;}
        for(int t=0;t<int(bad.terms.size());++t)if(bad.transitionCoefficient(bp.states[at],bp.states[at+1],t)!=0.){resolvedPath.events[at].term=t;visit(at+1);}};visit(0);
    require(std::abs(resolved-negative)<1e-12,"Combining negative-domain terms lost their signs");
    auto good=FredkinModel::chain(6);auto go=oracle(good,false);r.add("positive_nonpalindrome_cycle",0,good,bp,good.weight(bp),physicalWord(good,bp,go));
    std::cerr<<"PASS Fredkin controls matrix_entries="<<entries<<" complete_words="<<count<<" six_fermion_ZXX=-1/8 full_Fock_copies=16 negative_combined_cycle_product=-1/8\n";
}
struct Step {std::uint64_t to;int term;};
static std::vector<Step> steps(const FredkinModel& m,std::uint64_t s,bool combined){std::vector<Step> out;if(combined){for(auto t:m.nextStates(s,-1))out.push_back({t,-1});}else for(int k=0;k<int(m.terms.size());++k)for(auto t:m.nextStates(s,k))out.push_back({t,k});return out;}
static FredkinPath randomClosed(const FredkinModel& m,std::mt19937_64& rng,int index) {
    FredkinPath p;p.beta=.6+.03*(index%8);auto mask=(std::uint64_t(1)<<m.sites)-1;auto start=JQModel::occupation(rng()&mask,rng()&mask,m.sites);p.states={start};bool combined=index%2;
    for(int k=0;k<3+index%4;++k){auto next=steps(m,p.states.back(),combined);if(next.empty())break;auto a=next[rng()%next.size()];p.events.push_back({0.,a.term});p.states.push_back(a.to);}
    auto end=p.states.back();std::queue<std::uint64_t> q;q.push(end);std::map<std::uint64_t,std::pair<std::uint64_t,Step>> prev;prev[end]={end,{end,0}};
    while(!q.empty() && !prev.count(start)){auto s=q.front();q.pop();for(auto a:steps(m,s,combined))if(!prev.count(a.to)){prev[a.to]={s,a};q.push(a.to);}}
    require(prev.count(start),"Fredkin random walk cannot close");std::vector<Step> close;for(auto s=start;s!=end;s=prev.at(s).first)close.push_back(prev.at(s).second);std::reverse(close.begin(),close.end());for(auto a:close){p.events.push_back({0.,a.term});p.states.push_back(a.to);}times(p);return p;
}
int main(int argc,char** argv) {
    try {
        mkl_set_num_threads(1);mkl_set_dynamic(0);int samples=32;unsigned seed=160503842;std::string csv;
        for(int i=1;i<argc;++i){std::string key=argv[i];require(i+1<argc,"Missing Fredkin CLI value");std::string value=argv[++i];if(key=="--samples")samples=std::stoi(value);else if(key=="--seed")seed=std::stoul(value);else if(key=="--csv")csv=value;else throw std::runtime_error("Unknown Fredkin option");}
        require(samples>0,"Positive Fredkin samples required");std::ofstream file;if(!csv.empty()){file.open(csv);require(file.good(),"Cannot open Fredkin CSV");}Records r(csv.empty()?std::cout:file);controls(r);std::mt19937_64 rng(seed);
        for(int sites:{4,6}){auto m=FredkinModel::chain(sites,.8,1.1);for(int i=0;i<sites;++i)m.fields[i]=.03*(i-1);auto o=oracle(m,sites<=4);int moved=0;
            for(int index=0;index<samples;++index){auto p=randomClosed(m,rng,index);if(!p.events.empty())++moved;auto ref=physicalWord(m,p,o);require(std::abs(ref-m.reference(p))<1e-11,"Fredkin independent physical word mismatch");r.add(sites==4?"four_spin_random":"six_spin_random",index,m,p,m.weight(p),ref);}require(moved>0,"Fredkin ensemble never moved");std::cerr<<"PASS Fredkin sites="<<sites<<" moved="<<moved<<'\n';}
        std::cerr<<"PASS fredkin rows="<<r.rows<<" positive="<<r.positive<<" zero="<<r.zero<<" negative="<<r.negative<<" max_relative_error="<<r.error<<" max_phase_error="<<r.phase<<" seed="<<seed<<'\n';
    }catch(const std::exception& e){std::cerr<<"FAIL fredkin: "<<e.what()<<'\n';return 1;}return 0;
}
