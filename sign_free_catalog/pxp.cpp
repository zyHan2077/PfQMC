#include "catalog/pxp_models.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <queue>
#include <random>
using namespace catalog;
namespace {
struct Records {
    std::ostream& out;int rows=0,pos=0,zero=0,neg=0;double error=0.,phase=0.;
    Records(std::ostream& o):out(o){out<<"case,sample,modes,real_pf,imag_pf,real_reference,imag_reference,relative_error,expected,word\n"<<std::setprecision(17);}
    void add(std::string name,int sample,int modes,DataType w,DataType ref,int expected,const std::vector<int>& word={}) {
        double e=std::abs(w-ref)/std::max(1.,std::abs(ref));error=std::max(error,e);
        require(std::isfinite(e)&&e<2e-9,"PXP signed Pf/reference mismatch: "+name);
        if(expected==0){require(ref==DataType(0.)&&std::abs(w)<1e-10,"Nonzero PXP zero control");++zero;}
        else {require(expected*ref.real()>0.&&std::abs(ref.imag())<1e-10*std::abs(ref),"Wrong PXP reference sign");
            require(std::abs(w)>0.,"PXP weight vanished");double p=std::abs(w/std::abs(w)-double(expected));require(p<1e-8,"Wrong signed PXP phase");
            if(expected>0){++pos;phase=std::max(phase,p);}else ++neg;}
        ++rows;out<<name<<','<<sample<<','<<modes<<','<<w.real()<<','<<w.imag()<<','<<ref.real()<<','<<ref.imag()<<','<<e<<','<<expected<<",\"";
        for(int i:word)out<<i<<';';out<<"\"\n";
    }
};
PXPPath path(const PXPModel& m,std::vector<std::uint64_t> excitations,std::uint64_t parity,double beta=.8) {
    PXPPath p;p.beta=beta;for(auto s:excitations)p.states.push_back(PXPModel::occupation(s,parity,m.sites));
    for(size_t k=1;k<p.states.size();++k){auto diff=excitations[k-1]^excitations[k];require(diff&&!(diff&(diff-1)),"Path must flip one spin");
        p.events.push_back({beta*double(k)/p.states.size(),__builtin_ctzll(diff)});}return p;
}
DataType physical(const MatType& h,const PXPPath& p) {
    DataType w=1.;double last=0.;for(size_t k=0;k<p.events.size();++k){auto e=p.events[k];w*=-h(p.states[k+1],p.states[k])*std::exp(-(e.time-last)*h(p.states[k],p.states[k]).real());last=e.time;}
    return w*std::exp(-(p.beta-last)*h(p.states.back(),p.states.back()).real());
}
struct WordOracle {
    std::vector<std::vector<int>> next;std::vector<std::vector<DataType>> amplitude;int dim;
    WordOracle(const PXPModel& m):dim(1<<m.modes()) {
        // Matrices come from canonical-Fock Majoranas, not from enabledSites,
        // occupation flips, Gaussian factors, or the path-reference function.
        for(int i=0;i<m.sites;++i){auto v=m.termMatrix(i);std::vector<int> to(dim,-1);std::vector<DataType> a(dim,0.);
            for(int s=0;s<dim;++s)for(int t=0;t<dim;++t)if(v(t,s)!=DataType(0.)){require(to[s]<0,"PXP oracle vertex is not a partial permutation");to[s]=t;a[s]=v(t,s);}
            next.push_back(to);amplitude.push_back(a);}
    }
    DataType weight(const PXPModel& m,const std::vector<int>& word,double beta) const {
        DataType tr=0.;for(int start=0;start<dim;++start){int state=start;DataType a=1.;
            for(int i:word){a*=amplitude[i][state];state=next[i][state];if(state<0)break;}if(state==start)tr+=a;}
        int order=0;for(int i:word)tr*=beta*m.omega[i]/(++order);return tr;
    }
};
void fixed(Records& r) {
    PXPModel m;auto h=m.exactHamiltonian(),spin=m.exactHamiltonian(false);
    require((h-h.adjoint()).norm()<1e-12,"PXP physical H is not Hermitian");
    Eigen::SelfAdjointEigenSolver<MatType> ef(h),es(spin);double copies=0.;
    for(int i=0;i<h.rows();++i)copies=std::max(copies,std::abs(ef.eigenvalues()(i)-es.eigenvalues()(i/16)));
    require(copies<2e-11,"PXP full Fock is not sixteen spin copies");
    int row=0;
    for(int s=0;s<h.rows();++s) {
        require(std::abs(h(s,s)-m.diagonalEnergy(s))<1e-12,"PXP diagonal polynomial mismatch");
        auto allowed=m.enabledSites(s);
        for(int t=0;t<h.rows();++t)if(s!=t){bool expected=false;for(int i:allowed)expected|=m.flipped(s,i)==std::uint64_t(t);
            require(expected==(h(t,s)!=DataType(0.)),"PXP full polynomial hopping support mismatch");
            if(expected){auto p=path(m,{m.excitations(s),m.excitations(t),m.excitations(s)},0);
                p.states={std::uint64_t(s),std::uint64_t(t),std::uint64_t(s)};
                require(std::abs(m.reference(p)-physical(h,p))<1e-11,"PXP local scalar/full Fock mismatch");
                r.add("all_two_vertex_full_Fock",row++,m.modes(),m.weight(p),physical(h,p),1);}}
    }
    auto flat=PXPModel::chain(4,1.);flat.repulsion=flat.detuning=0.;auto hf=flat.exactHamiltonian();
    for(int parity=0;parity<16;++parity){auto p=path(flat,{0,1,5,4,0},parity);
        r.add("nonpalindromic_four_flips_all_parities",parity,flat.modes(),flat.weight(p),physical(hf,p),1);}
    // The all-excited state is retained in full Hilbert space, even though all
    // flips are frozen. Its empty CT word has positive diagonal Boltzmann weight.
    auto frozen=path(m,{15},5);r.add("full_Hilbert_frozen_adjacent_excitations",0,m.modes(),m.weight(frozen),physical(h,frozen),1);
    auto part=path(m,{3,11,3},9);r.add("dynamic_region_next_to_frozen_pair",0,m.modes(),m.weight(part),physical(h,part),1);
    auto blocked=path(m,{3,2,3},0);r.add("zero_neighbor_blocked_flip",0,m.modes(),m.weight(blocked),0.,0);
    auto disabled=flat;disabled.omega[0]=0.;auto p=path(disabled,{0,1,0},0);r.add("zero_amplitude",0,disabled.modes(),disabled.weight(p),0.,0);
    auto open=path(flat,{0,1},0);r.add("zero_nonclosed_path",0,flat.modes(),flat.weight(open),0.,0);
    auto flux=flat;flux.phaseFlux=true;require(!flux.signFreeParameters(),"Conditional flux control admitted to positive domain");
    auto bad=path(flux,{0,1,5,4,0},10);auto hb=flux.exactHamiltonian();
    r.add("negative_conditional_flux_closed_loop",0,flux.modes(),flux.weight(bad),physical(hb,bad),-1);
    WordOracle fluxOracle(flux);std::vector<int> fluxWord{0,2,0,2};
    r.add("negative_complete_flux_word",0,flux.modes(),flux.kineticWordWeight(fluxWord,.8),fluxOracle.weight(flux,fluxWord,.8),-1,fluxWord);
    // Changing a constant site's Omega sign is only a local Z gauge. A closed
    // loop remains positive; it is deliberately NOT used as a negative control.
    auto gauge=flat;gauge.omega[0]=-1.;auto gp=path(gauge,{0,1,5,4,0},3);
    r.add("constant_negative_Omega_gauge_control",0,gauge.modes(),gauge.weight(gp),physical(gauge.exactHamiltonian(),gp),1);
    // Two X flips: the Gaussian rotations alone have Spin lift -I, while
    // the two retained -i scalars restore the physical positive X^2 trace.
    auto half=flat.halfRotation(0);std::vector<MatType> lift{half,half,half,half};
    r.add("negative_raw_Gaussian_Spin_lift",0,flat.modes(),traceProduct(lift,true),-std::ldexp(1.,flat.modes()),-1);
    r.add("zero_single_kinetic_vertex_trace",0,flat.modes(),flat.kineticWordWeight({0},.8),0.,0,{0});
    // Full kinetic CT expansion, all ordered words through order four.
    WordOracle oracle(flat);std::vector<DataType> sums(5,0.);int words=0;
    for(int order=0;order<=4;++order){int count=1<<(2*order);
        for(int w=0;w<count;++w){std::vector<int> labels(order);int x=w;for(int k=0;k<order;++k){labels[k]=x%4;x/=4;}
            auto ref=oracle.weight(flat,labels,.8),pf=flat.kineticWordWeight(labels,.8);sums[order]+=pf;
            r.add("complete_kinetic_word",words++,flat.modes(),pf,ref,ref.real()>0?1:0,labels);}}
    MatType power=MatType::Identity(16,16),kinetic=-flat.exactHamiltonian(false);double coefficient=1.;
    for(int order=0;order<=4;++order){if(order){power=kinetic*power;coefficient*=.8/order;}
        require(std::abs(sums[order]-16.*coefficient*power.trace())<1e-9,"PXP complete word sum disagrees with physical H power");}
    // Full interacting H: exact diagonal propagation and off-diagonal vertices
    // converge to the exact partition function as the time grid is refined.
    MatType d=spin.diagonal().asDiagonal(),k=d-spin;double exactZ=0.;for(double e:es.eigenvalues())exactZ+=std::exp(-.8*e);
    double previous=1e99;for(int slices:{32,64,128}){MatType step=(MatType::Identity(16,16)+(.8/slices)*k)*matrixExp(-(.8/slices)*d),prod=MatType::Identity(16,16);
        for(int j=0;j<slices;++j)prod=step*prod;double error=std::abs(prod.trace().real()-exactZ);require(error<previous,"PXP full-H time slicing failed to converge");previous=error;}
    MatType n=MatType::Zero(16,16);for(int s=0;s<16;++s)n(s,s)=__builtin_popcount(unsigned(s));
    double numberComm=(spin*n-n*spin).norm(),vertexComm=(flat.termMatrix(0,false)*flat.termMatrix(1,false)-flat.termMatrix(1,false)*flat.termMatrix(0,false)).norm();
    require(numberComm>1.&&vertexComm>1.,"PXP lost dynamic number changes or noncommuting vertices");
    std::cerr<<"PXP physical_copy_error="<<copies<<" number_commutator="<<numberComm<<" adjacent_vertex_commutator="<<vertexComm<<'\n';
}
std::vector<std::uint64_t> close(const PXPModel& m,std::uint64_t from,std::uint64_t to) {
    int dim=1<<m.sites;std::vector<int> previous(dim,-1);std::queue<int> q;q.push(from);previous[from]=from;
    while(!q.empty()&&previous[to]<0){int state=q.front();q.pop();auto occ=PXPModel::occupation(state,0,m.sites);
        for(int i:m.enabledSites(occ)){int next=state^(1<<i);if(previous[next]<0){previous[next]=state;q.push(next);}}}
    require(previous[to]>=0,"PXP closure path missing");std::vector<std::uint64_t> suffix;
    for(int x=to;x!=int(from);x=previous[x])suffix.push_back(x);std::reverse(suffix.begin(),suffix.end());return suffix;
}
void randomPaths(Records& r,int samples,std::mt19937_64& rng) {
    std::uniform_real_distribution<double> amp(.3,1.2),diag(-.5,.5),time(.01,.79);int nonpal=0,multi=0;
    for(int n:{4,6})for(int sample=0;sample<samples;++sample){auto m=PXPModel::chain(n);for(double& x:m.omega)x=amp(rng);m.repulsion=diag(rng);m.detuning=diag(rng);
        std::vector<std::uint64_t> walk{rng()%(1u<<n)};for(int k=0;k<4+sample%5;++k){auto allowed=m.enabledSites(PXPModel::occupation(walk.back(),0,n));if(allowed.empty())break;walk.push_back(walk.back()^(1u<<allowed[rng()%allowed.size()]));}
        auto suffix=close(m,walk.back(),walk.front());walk.insert(walk.end(),suffix.begin(),suffix.end());auto p=path(m,walk,rng()%(1u<<n));
        std::vector<double> times(p.events.size());for(double& t:times)t=time(rng);std::sort(times.begin(),times.end());for(size_t k=0;k<times.size();++k)p.events[k].time=times[k];
        nonpal+=!std::equal(walk.begin(),walk.end(),walk.rbegin());unsigned updates=0;for(auto e:p.events)updates|=1u<<e.site;multi+=__builtin_popcount(updates)>1;
        auto ref=n==4?physical(m.exactHamiltonian(),p):m.reference(p);r.add("random_full_Hilbert_closed_path",sample,m.modes(),m.weight(p),ref,1);}
    if(samples>=32)require(nonpal>0&&multi>0,"PXP ensemble missed nonretracing multispin paths");std::cerr<<"PXP random_nonpal="<<nonpal<<" random_multispin="<<multi<<'\n';
}
}
int main(int argc,char** argv){try {
    int samples=32;unsigned seed=731;std::string csv;for(int i=1;i<argc;++i){std::string arg=argv[i];require(i+1<argc,"Missing PXP argument value");std::string value=argv[++i];
        if(arg=="--samples")samples=std::stoi(value);else if(arg=="--seed")seed=std::stoul(value);else if(arg=="--csv")csv=value;else throw std::runtime_error("Unknown PXP argument");}
    require(samples>0,"PXP samples must be positive");mkl_set_num_threads(1);mkl_set_dynamic(0);std::ofstream file;if(!csv.empty()){file.open(csv);require(file.good(),"Cannot open PXP CSV");}
    Records r(csv.empty()?std::cout:file);std::mt19937_64 rng(seed);fixed(r);randomPaths(r,samples,rng);r.out.flush();require(r.out.good(),"PXP CSV write failed");
    std::cerr<<std::setprecision(17)<<"PASS PXP rows="<<r.rows<<" positive="<<r.pos<<" zero="<<r.zero<<" negative_controls="<<r.neg<<" max_relative_error="<<r.error
        <<" max_phase_error="<<r.phase<<" samples="<<samples<<" seed="<<seed<<'\n';return 0;
}catch(const std::exception& e){std::cerr<<"FAIL PXP: "<<e.what()<<'\n';return 1;}}
