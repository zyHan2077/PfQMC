#include "catalog/dynamic_ising.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>

using namespace catalog;
namespace {
IsingConfiguration resolve(double beta,const std::vector<IsingEvent>& events,std::uint64_t initial) {
    IsingConfiguration c{beta,events,{initial}};
    for(const auto& e:events) {
        auto next=c.spinPath.back();
        if(e.kind==IsingVertex::SpinFlip) next^=std::uint64_t(3)<<(2*e.index);
        c.spinPath.push_back(next);
    }
    return c;
}
struct Recorder {
    std::ostream& out;
    std::uint64_t seed;
    int rows=0,positive=0,zero=0,negative=0;
    double maxError=0.,maxPhase=0.;
    Recorder(std::ostream& out,std::uint64_t seed):out(out),seed(seed) {
        out<<"case,sample,seed,sites,events,weight_real,weight_imag,reference_real,reference_imag,sign_real,sign_imag,expected,scaled_error,beta,event_word,spin_path,hopping,repulsion,ising,transverse,edges\n";
        out<<std::setprecision(17);
    }
    void add(const std::string& name,int sample,int sites,DataType w,DataType ref,int sign,
             const IsingConfiguration* c=nullptr,const DynamicIsing* model=nullptr) {
        require(std::isfinite(std::abs(w)) && std::isfinite(std::abs(ref)),"Nonfinite result: "+name);
        double error=std::abs(w-ref)/std::max(1.,std::abs(ref));
        require(error<2e-9,"Pfaffian/reference disagreement: "+name);
        maxError=std::max(maxError,error);
        DataType phase=std::abs(w)>0. ? w/std::abs(w):DataType(std::numeric_limits<double>::quiet_NaN(),0.);
        if(sign==0) {
            require(std::abs(w)<2e-10 && std::abs(ref)<2e-10,"Zero control has nonzero weight: "+name);++zero;
        } else {
            require(sign*ref.real()>0. && std::abs(ref.imag())<1e-10*std::max(1.,std::abs(ref)),"Invalid signed reference: "+name);
            require(std::abs(w)>1e-14,"Unresolved tiny signed weight: "+name);
            require(std::abs(phase-double(sign))<1e-8,"Incorrect measured sign: "+name);
            if(sign>0) {++positive;maxPhase=std::max(maxPhase,std::abs(phase-1.));} else ++negative;
        }
        ++rows;
        out<<name<<','<<sample<<','<<seed<<','<<sites<<','<<(c ? c->events.size():0)<<','
           <<w.real()<<','<<w.imag()<<','<<ref.real()<<','<<ref.imag()<<','<<phase.real()<<','<<phase.imag()<<','
           <<sign<<','<<error<<','<<(c ? c->beta:0.)<<",\"";
        if(c) for(auto e:c->events) out<<e.time<<':'<<(e.kind==IsingVertex::SpinFlip ? 's':'f')<<':'<<e.index<<';';
        out<<"\",\"";if(c) for(auto v:c->spinPath) out<<v<<';';out<<"\",";
        if(model) {
            out<<model->hopping<<','<<model->repulsion<<','<<model->ising<<",\"";
            for(double h:model->transverse) out<<h<<';';out<<"\",\"";
            for(auto [i,j]:model->edges) out<<i<<':'<<j<<';';out<<'"';
        } else out<<",,,\"\",\"\"";
        out<<'\n';
    }
};
void mustThrow(const std::function<void()>& fn,const std::string& what) {
    bool threw=false;try {fn();} catch(const std::exception&) {threw=true;}require(threw,"Missing guard: "+what);
}
void controls(Recorder& r) {
    DynamicIsing m;auto s=IsingVertex::SpinFlip;
    std::vector<IsingEvent> ev{{.3,s,0},{.7,s,0}};auto c=resolve(1.2,ev,0);
    r.add("spin_projector_closed_loop",0,2,m.spinPathPfaffian(c),1.,1,&c,&m);
    // Deliberately reverse the first pair rotation: the physical matrix element
    // acquires a minus sign. This probes the signed Pfaffian, not abs(weight).
    auto wrong=traceGaussianFactors({m.spinProjector(0),m.spinFlipRotation(0,false),
                                    m.spinProjector(3),m.spinFlipRotation(0,false)});
    r.add("negative_reversed_pair_rotation",0,2,wrong,-1.,-1);
    MatType h=MatType::Zero(2,2);h(0,1)=h(1,0)=.8;
    MatType n0=MatType::Zero(2,2),n1=n0;n0(0,0)=1.;n1(1,1)=1.;
    const double pi=std::acos(-1.);
    std::vector<MatType> pulses{numberMajorana(.3*h),numberMajorana(.5i*pi*n0),numberMajorana(.5i*pi*n0),
        numberMajorana(.4*h),numberMajorana(.5i*pi*n1),numberMajorana(.5i*pi*n1),numberMajorana(.5*h)};
    // This fermion-only density word omits the staggered factor on purpose.
    // Its corresponding Ising path is nonclosed and has total physical weight 0.
    DataType prefactor=.6*.4*DataType(-.5i)*DataType(-.5i);
    r.add("negative_unstaggered_density_subtrace",0,2,prefactor*traceProduct(pulses,true),
          prefactor*exactTrace(pulses),-1);
    auto zero=[&](const std::string& name,const IsingConfiguration& x) {
        r.add(name,0,2,m.configurationPfaffian(x),m.configurationReference(x),0,&x,&m);
    };
    zero("zero_blocked_single_occupancy",resolve(1.2,ev,1));
    zero("zero_nonclosed_path",resolve(1.2,{{.3,s,0},{.7,s,1}},0));
    auto bad=c;bad.events[0].time=-.1;zero("zero_negative_time",bad);
    bad=c;bad.events[0].time=.9;zero("zero_unsorted_time",bad);
    bad=c;bad.events[1].time=1.3;zero("zero_time_after_beta",bad);
    bad=c;bad.events[0].time=std::numeric_limits<double>::quiet_NaN();zero("zero_nonfinite_time",bad);
    bad=c;bad.spinPath.pop_back();zero("zero_path_size",bad);
    bad=c;bad.events[0].index=2;zero("zero_event_site",bad);
    bad=c;bad.spinPath[0]=16;zero("zero_occupation_outside_space",bad);
    bad=c;bad.beta=0.;zero("zero_nonpositive_beta",bad);
    bad=resolve(1.2,{{.3,IsingVertex::DensityBond,0}},0);bad.spinPath.back()=3;
    zero("zero_density_changes_spin",bad);
    r.add("zero_odd_spin_word",0,2,m.factorizedWordTrace(1.2,{{.3,s,0}}),m.exactWordTrace(1.2,{{.3,s,0}}),0);
    auto empty=resolve(1.2,{},1);r.add("blocked_spectator_empty_word",0,2,m.configurationPfaffian(empty),m.configurationReference(empty),1,&empty,&m);
    auto equal=resolve(1.2,{{.3,s,0},{.3,s,0}},0);
    r.add("equal_time_ordered_vertices",0,2,m.configurationPfaffian(equal),m.configurationReference(equal),1,&equal,&m);
    DynamicIsing absent=m;absent.transverse[0]=0.;
    r.add("zero_absent_transverse_coupling",0,2,absent.configurationPfaffian(c),absent.configurationReference(c),0,&c,&absent);
    absent=m;absent.repulsion=0.;auto noV=resolve(1.2,{{.3,IsingVertex::DensityBond,0}},0);
    r.add("zero_absent_density_coupling",0,2,absent.configurationPfaffian(noV),absent.configurationReference(noV),0,&noV,&absent);
    DynamicIsing badGraph=m;badGraph.sites=3;badGraph.transverse={.4,.5,.6};badGraph.edges={{0,1},{1,2},{2,0}};
    mustThrow([&]{badGraph.validate();},"odd cycle");badGraph=m;badGraph.edges.push_back({1,0});
    mustThrow([&]{badGraph.validate();},"duplicate edge");badGraph=m;badGraph.edges={{0,0}};
    mustThrow([&]{badGraph.validate();},"self edge");badGraph=m;badGraph.repulsion=-.1;
    mustThrow([&]{badGraph.validate();},"negative V");
    auto physical=m.exactHamiltonian(),rotated=m.exactHamiltonian(true);
    require((physical-physical.adjoint()).norm()<1e-12,"Non-Hermitian physical H");
    // U=exp(i*pi*S_z on the negative sublattice), built independently by occupations.
    cVecType u=cVecType::Ones(physical.rows());
    for(int v=0;v<physical.rows();++v) {
        double sz=(((v>>4)&1)+((v>>5)&1)-1.)/2.;u(v)=std::exp(1.i*pi*sz);
    }
    MatType transformed=u.conjugate().asDiagonal()*physical*u.asDiagonal();
    require((transformed-rotated).norm()<1e-12,"Staggered basis rotation changed the physical Hamiltonian");
    DynamicIsing disconnected;disconnected.sites=5;disconnected.edges={{0,2},{3,4}};
    disconnected.transverse={.6,.3,-.4,.7,.5};auto colors=disconnected.sublattices();
    require(colors[0]==-colors[2] && colors[3]==-colors[4],"BFS confused labels with bipartite colors");
    auto d=resolve(1.2,{{.3,s,2},{.7,s,2}},0);
    r.add("disconnected_nonparity_labels",0,5,disconnected.configurationPfaffian(d),disconnected.configurationReference(d),1,&d,&disconnected);
    DynamicIsing star=DynamicIsing::square();star.edges={{0,1},{0,2},{0,3}};
    auto x=resolve(1.2,{{.2,IsingVertex::DensityBond,1},{.4,s,2},{.8,s,2}},0);
    r.add("bipartite_branched_graph",0,4,star.configurationPfaffian(x),star.configurationReference(x),1,&x,&star);
}
void fixedWords(Recorder& r) {
    auto f=IsingVertex::DensityBond,s=IsingVertex::SpinFlip;
    for(auto m:{DynamicIsing{},DynamicIsing::square()}) {
        std::vector<IsingEvent> ev=m.sites==2 ? std::vector<IsingEvent>{{.13,f,0},{.27,s,0},{.49,s,1},{.78,s,0},{1.03,s,1}}:
            std::vector<IsingEvent>{{.08,f,0},{.23,s,0},{.41,s,1},{.63,f,2},{.82,s,0},{1.01,s,1}};
        DataType sum=0.;int paths=0;
        for(std::uint64_t v=0;v<(std::uint64_t(1)<<(2*m.sites));++v) {
            auto c=resolve(1.2,ev,v);if(!m.validConfiguration(c)) continue;
            auto w=m.configurationPfaffian(c),ref=m.configurationReference(c);sum+=w;
            r.add("fixed_resolved_path",paths++,m.sites,w,ref,1,&c,&m);
        }
        require(paths==(m.sites==2 ? 4:64),"Unexpected number of resolved spin paths");
        auto word=IsingConfiguration{1.2,ev,{}};
        r.add("fixed_sum_vs_factorized",0,m.sites,sum,m.factorizedWordTrace(1.2,ev),1,&word,&m);
        if(m.sites==2) r.add("fixed_sum_vs_full_Fock",0,2,sum,m.exactWordTrace(1.2,ev),1,&word,&m);
    }
}
void randomWords(Recorder& r,int samples,std::mt19937_64& rng) {
    std::uniform_real_distribution<double> time(.04,1.16),coupling(.4,.9);
    for(auto m:{DynamicIsing{},DynamicIsing::square()}) for(int sample=0;sample<samples;++sample) {
        for(double& h:m.transverse) h=coupling(rng)*(rng()%2 ? 1.:-1.);
        m.hopping=coupling(rng);m.repulsion=coupling(rng);m.ising=(rng()%2 ? 1.:-1.)*coupling(rng);
        std::vector<IsingEvent> ev;
        int bonds=rng()%3,pairs=1+rng()%2;
        for(int k=0;k<bonds;++k) ev.push_back({time(rng),IsingVertex::DensityBond,int(rng()%m.edges.size())});
        for(int k=0;k<pairs;++k) {
            int site=rng()%m.sites;ev.push_back({time(rng),IsingVertex::SpinFlip,site});ev.push_back({time(rng),IsingVertex::SpinFlip,site});
        }
        std::stable_sort(ev.begin(),ev.end(),[](auto a,auto b){return a.time<b.time;});
        std::uint64_t v=0;for(int i=0;i<m.sites;++i) if(rng()%2) v|=std::uint64_t(3)<<(2*i);
        auto c=resolve(1.2,ev,v);
        r.add("random_resolved_path",sample,m.sites,m.configurationPfaffian(c),m.configurationReference(c),1,&c,&m);
        if(sample<2) r.add("random_spin_Pfaffian",sample,m.sites,m.spinPathPfaffian(c),1.,1,&c,&m);
    }
}
}
int main(int argc,char** argv) {
    try {
        int samples=8;std::uint64_t seed=160507420;std::string csv;
        for(int i=1;i<argc;++i) {
            std::string a=argv[i];
            if(a=="--help") {
                std::cout<<"sign_free_dynamic_ising --samples N --seed N --csv PATH\n"
                    <<"Signed CT-INT weight verifier for the local three-orbital dynamic Ising-fermion model.\n"
                    <<"Includes projectors, all expansion scalars, Fock references, negative and zero controls.\n"
                    <<"This executable is not a production Monte Carlo updater.\n";return 0;
            }
            require(i+1<argc,"Missing argument value");std::string value=argv[++i];
            if(a=="--samples") samples=std::stoi(value);
            else if(a=="--seed") seed=std::stoull(value);
            else if(a=="--csv") csv=value;
            else throw std::runtime_error("Unknown argument: "+a);
        }
        require(samples>0,"samples must be positive");
        mkl_set_num_threads(1);mkl_set_dynamic(0);
        std::ofstream file;if(!csv.empty()) {file.open(csv);require(file.good(),"Cannot open CSV");}
        Recorder r(csv.empty() ? std::cout:file,seed);std::mt19937_64 rng(seed);
        controls(r);fixedWords(r);randomWords(r,samples,rng);
        r.out.flush();require(r.out.good(),"Failed to write validation output");
        std::cerr<<std::setprecision(15)<<"PASS DynamicIsing rows="<<r.rows<<" positive="<<r.positive<<" zero="<<r.zero
            <<" negative_controls="<<r.negative<<" max_scaled_error="<<r.maxError<<" max_positive_phase_error="<<r.maxPhase
            <<" samples_per_size="<<samples<<" seed="<<seed<<"\n";
    } catch(const std::exception& e) {std::cerr<<"FAIL DynamicIsing: "<<e.what()<<"\n";return 1;}
    return 0;
}
