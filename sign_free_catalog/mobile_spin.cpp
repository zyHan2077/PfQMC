#include "catalog/mobile_spin.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>

using namespace catalog;
namespace {
constexpr auto density=MobileVertex::DensityBond,flip=MobileVertex::SpinFlip,hop=MobileVertex::PairHop;
struct Recorder {
    std::ostream& out;std::uint64_t seed;
    int rows=0,positive=0,zero=0,negative=0;double maxError=0.,maxPhase=0.;
    Recorder(std::ostream& output,std::uint64_t seed):out(output),seed(seed) {
        out<<"case,sample,seed,sites,events,weight_real,weight_imag,reference_real,reference_imag,sign_real,sign_imag,expected,scaled_error,beta,event_word,spin_path,hopping,repulsion,longitudinal,exchange,transverse,edges,sign_free_parameters\n"<<std::setprecision(17);
    }
    void add(const std::string& name,int sample,const MobileSpin& m,DataType w,DataType ref,int expected,
             const MobileConfiguration* c=nullptr) {
        require(std::isfinite(std::abs(w)) && std::isfinite(std::abs(ref)),"Nonfinite weight: "+name);
        double error=std::abs(w-ref)/std::max(1.,std::abs(ref));maxError=std::max(maxError,error);
        require(error<2e-9,"Pfaffian/reference disagreement: "+name);
        DataType phase=std::abs(w)>0. ? w/std::abs(w):DataType(std::numeric_limits<double>::quiet_NaN(),0.);
        if(expected==0) {
            require(std::abs(w)<2e-10 && std::abs(ref)<2e-10,"Zero control nonzero: "+name);++zero;
        } else {
            require(expected*ref.real()>0. && std::abs(ref.imag())<1e-10*std::max(1.,std::abs(ref)),"Incorrect reference sign: "+name);
            require(std::abs(w)>1e-14,"Unresolved tiny signed weight: "+name);
            require(std::abs(phase-double(expected))<1e-8,"Incorrect measured phase: "+name);
            if(expected>0) {++positive;maxPhase=std::max(maxPhase,std::abs(phase-1.));} else ++negative;
        }
        ++rows;out<<name<<','<<sample<<','<<seed<<','<<m.sites<<','<<(c ? c->events.size():0)<<','
           <<w.real()<<','<<w.imag()<<','<<ref.real()<<','<<ref.imag()<<','<<phase.real()<<','<<phase.imag()<<','
           <<expected<<','<<error<<','<<(c ? c->beta:0.)<<",\"";
        if(c) for(auto e:c->events) out<<e.time<<':'<<(e.kind==flip ? 's':e.kind==hop ? 'p':'f')<<':'<<e.index<<';';
        out<<"\",\"";if(c) for(auto v:c->spinPath) out<<v<<';';
        out<<"\","<<m.hopping<<','<<m.repulsion<<','<<m.longitudinal<<','<<m.exchange<<",\"";
        for(double h:m.transverse) out<<h<<';';out<<"\",\"";
        for(auto [i,j]:m.edges) out<<i<<':'<<j<<';';out<<"\","<<m.signFreeParameters()<<'\n';
    }
    void configuration(const std::string& name,int sample,const MobileSpin& m,const MobileConfiguration& c,int expected=1) {
        add(name,sample,m,m.configurationPfaffian(c),m.configurationReference(c),expected,&c);
    }
};
void mustThrow(const std::function<void()>& fn,const std::string& name) {
    bool threw=false;try {fn();} catch(const std::exception&) {threw=true;}
    require(threw,"Missing parameter guard: "+name);
}
void fixedWords(Recorder& r) {
    for(auto m:{MobileSpin{},MobileSpin::square()}) {
        std::vector<MobileEvent> events=m.sites==2 ? std::vector<MobileEvent>{{.2,flip,0},{.5,hop,0},{.9,flip,1}}:
            std::vector<MobileEvent>{{.1,density,2},{.2,flip,0},{.45,hop,0},{.65,hop,1},{1.,flip,2}};
        DataType sum=0.;int count=0;
        for(std::uint64_t v=0;v<(std::uint64_t(1)<<(2*m.sites));++v) {
            auto c=m.resolve(1.2,events,v);if(!m.validConfiguration(c)) continue;
            auto w=m.configurationPfaffian(c);sum+=w;
            r.add("fixed_mobile_closed_path",count,m,w,m.configurationReference(c),1,&c);
            if(count<2) r.add("fixed_spin_Pfaffian",count,m,m.spinPathPfaffian(c),1.,1,&c);
            ++count;
        }
        require(count==(m.sites==2 ? 2:8),"Unexpected closed mobile path count");
        MobileConfiguration word{1.2,events,{}};
        r.add("fixed_sum_vs_factorized",0,m,sum,m.factorizedWordTrace(1.2,events),1,&word);
        if(m.sites==2) {
            r.add("fixed_sum_vs_original_6mode_Fock",0,m,sum,m.exactWordTrace(1.2,events),1,&word);
            require(std::abs(sum-0.0003586743649135766)<1e-12,"Independent Python fixed value mismatch");
        }
    }
}
void controls(Recorder& r) {
    MobileSpin m;std::vector<MobileEvent> ev{{.2,flip,0},{.5,hop,0},{.9,flip,1}};
    auto c=m.resolve(1.2,ev,0);
    // Each site has ONE Sx vertex, allowed because the pair moves across the bond.
    require(m.validConfiguration(c) && c.spinPath==std::vector<std::uint64_t>({0,3,12,0}),"Mobile pair did not close");
    auto badSign=m;badSign.transverse[1]*=-1.;require(!badSign.signFreeParameters(),"Mixed h sign not identified");
    r.configuration("negative_relative_h_sign",0,badSign,c,-1);
    auto word=MobileConfiguration{1.2,ev,{}};
    r.add("negative_h_sum_vs_original_6mode_Fock",0,badSign,badSign.factorizedWordTrace(1.2,ev),badSign.exactWordTrace(1.2,ev),-1,&word);
    auto both=m;for(double& h:both.transverse) h=-h;require(both.signFreeParameters(),"Both negative h should be allowed");
    r.configuration("both_h_negative",0,both,c);
    auto wrong=traceGaussianFactors({m.spinProjector(0),m.spinFlipRotation(0,false),m.spinProjector(3),m.spinFlipRotation(0,false)});
    r.add("negative_reversed_pair_rotation",0,m,wrong,-1.,-1);
    auto zero=[&](const std::string& name,const MobileConfiguration& x){r.configuration(name,0,m,x,0);};
    zero("zero_blocked_source",m.resolve(1.2,ev,1));
    zero("zero_blocked_pair_hop",m.resolve(1.2,ev,4));
    zero("zero_hop_empty_pair",m.resolve(1.2,{{.3,hop,0},{.8,hop,0}},0));
    zero("zero_hop_full_pair",m.resolve(1.2,{{.3,hop,0},{.8,hop,0}},15));
    zero("zero_nonclosed",m.resolve(1.2,{{.2,flip,0},{.9,flip,1}},0));
    auto bad=c;bad.events[0].time=-.1;zero("zero_negative_time",bad);
    bad=c;bad.events[0].time=.8;zero("zero_unsorted_time",bad);
    bad=c;bad.events.back().time=1.3;zero("zero_after_beta",bad);
    bad=c;bad.events[0].time=std::numeric_limits<double>::quiet_NaN();zero("zero_nonfinite_time",bad);
    bad=c;bad.events[1].index=1;zero("zero_event_edge",bad);
    bad=c;bad.events[0].index=2;zero("zero_event_site",bad);
    bad=c;bad.spinPath.pop_back();zero("zero_path_length",bad);
    bad=c;bad.spinPath[1]=15;zero("zero_wrong_transition",bad);
    bad=c;bad.spinPath[0]=16;zero("zero_occupation_bound",bad);
    bad=c;bad.beta=0.;zero("zero_nonpositive_beta",bad);
    bad=m.resolve(1.2,{{.3,density,0}},0);bad.spinPath.back()=3;zero("zero_density_changes_spin",bad);
    auto absent=m;absent.exchange=0.;r.configuration("zero_absent_exchange",0,absent,c,0);
    absent=m;absent.transverse[0]=0.;r.configuration("zero_absent_source",0,absent,c,0);
    absent=m;absent.repulsion=0.;r.configuration("zero_absent_density",0,absent,m.resolve(1.2,{{.3,density,0}},0),0);
    auto empty=m.resolve(1.2,{},1);r.configuration("spin_zero_spectator_empty_word",0,m,empty);
    auto equal=m.resolve(1.2,{{.3,flip,0},{.3,hop,0},{.3,flip,1}},0);
    r.configuration("equal_time_vertex_order",0,m,equal);
    auto cancelled=m.resolve(1.2,{{.2,flip,0},{.5,hop,0},{.8,flip,1}},0);
    r.configuration("zero_half_period_density_trace",0,m,cancelled,0);
    auto physical=m.exactHamiltonian(),rotated=m.exactHamiltonian(true);
    require((physical-physical.adjoint()).norm()<1e-12,"Physical H is not Hermitian");
    cVecType u=cVecType::Ones(physical.rows());double pi=std::acos(-1.);auto color=m.sublattices();
    for(int v=0;v<physical.rows();++v) for(int i=0;i<m.sites;++i) if(color[i]<0) {
        double sz=(((v>>(m.sites+2*i))&1)+((v>>(m.sites+2*i+1))&1)-1.)/2.;u(v)*=std::exp(1.i*pi*sz);
    }
    require((u.conjugate().asDiagonal()*physical*u.asDiagonal()-rotated).norm()<1e-12,"Hamiltonian Marshall rotation mismatch");
    MobileSpin disconnected;disconnected.sites=5;disconnected.edges={{0,2},{3,4}};
    disconnected.transverse={.6,.2,.4,-.5,-.7};
    require(disconnected.signFreeParameters(),"Independent components may have opposite h signs");
    auto colors=disconnected.sublattices();require(colors[0]==-colors[2],"Used numeric label parity instead of BFS");
    r.configuration("disconnected_nonparity_graph",0,disconnected,disconnected.resolve(1.2,{{.2,flip,0},{.5,hop,0},{.9,flip,2}},0));
    auto star=MobileSpin::square();star.edges={{0,1},{0,2},{0,3}};
    r.configuration("branched_nonparity_graph",0,star,star.resolve(1.2,{{.2,flip,0},{.5,hop,1},{.9,flip,2}},0));
    MobileSpin invalid=m;invalid.sites=3;invalid.transverse={.3,.4,.5};invalid.edges={{0,1},{1,2},{2,0}};
    mustThrow([&]{invalid.validate();},"odd cycle");invalid=m;invalid.edges={{0,0}};
    mustThrow([&]{invalid.validate();},"self edge");invalid=m;invalid.edges.push_back({1,0});
    mustThrow([&]{invalid.validate();},"duplicate edge");invalid=m;invalid.exchange=-.1;
    mustThrow([&]{invalid.validate();},"negative exchange");invalid=m;invalid.repulsion=-.1;
    mustThrow([&]{invalid.validate();},"negative V");
    // Verify the Ising overlap explicitly, rather than counting Jxy=0 again.
    MobileSpin limit=m;limit.exchange=0.;limit.transverse={.6,-.4};require(limit.signFreeParameters(),"Ising limit h signs should be unrestricted");
    DynamicIsing ising;ising.ising=-limit.longitudinal;ising.transverse=limit.transverse;
    require((limit.exactHamiltonian()-ising.exactHamiltonian()).norm()<1e-12,"Ising-limit Hamiltonian mismatch");
    std::vector<MobileEvent> lev{{.13,density,0},{.27,flip,0},{.49,flip,1},{.78,flip,0},{1.03,flip,1}};
    std::vector<IsingEvent> iev;for(auto e:lev) iev.push_back({e.time,e.kind==flip ? IsingVertex::SpinFlip:IsingVertex::DensityBond,e.index});
    auto lc=limit.resolve(1.2,lev,0);IsingConfiguration ic{1.2,iev,lc.spinPath};
    r.add("Ising_limit_signed_overlap",0,limit,limit.configurationPfaffian(lc),ising.configurationPfaffian(ic),1,&lc);
}
void randomWords(Recorder& r,int samples,std::mt19937_64& rng) {
    std::uniform_real_distribution<double> coupling(.4,.9),time(.04,1.16);
    for(auto m:{MobileSpin{},MobileSpin::square()}) for(int sample=0;sample<samples;++sample) {
        int sign=rng()%2 ? 1:-1;for(double& h:m.transverse) h=sign*coupling(rng);
        m.hopping=coupling(rng);m.repulsion=coupling(rng);m.exchange=coupling(rng);
        m.longitudinal=(rng()%2 ? -1.:1.)*coupling(rng);
        int current=rng()%m.sites;std::vector<MobileEvent> events{{0.,flip,current}};
        int steps=1+rng()%3;
        for(int k=0;k<steps;++k) {
            std::vector<int> incident;for(int e=0;e<int(m.edges.size());++e)
                if(m.edges[e].first==current || m.edges[e].second==current) incident.push_back(e);
            int e=incident[rng()%incident.size()];events.push_back({0.,hop,e});
            auto [i,j]=m.edges[e];current=(current==i ? j:i);
        }
        events.push_back({0.,flip,current});
        std::vector<double> times(events.size());for(double& t:times) t=time(rng);std::sort(times.begin(),times.end());
        for(size_t k=0;k<events.size();++k) events[k].time=times[k];
        if(rng()%2) events.push_back({time(rng),density,int(rng()%m.edges.size())});
        std::stable_sort(events.begin(),events.end(),[](auto a,auto b){return a.time<b.time;});
        auto initial=rng()%2 ? std::uint64_t(0):(std::uint64_t(1)<<(2*m.sites))-1;
        auto c=m.resolve(1.2,events,initial);require(m.validConfiguration(c),"Generated mobile walk is invalid");
        require(m.signFreeParameters(),"Generated parameters outside proof");r.configuration("random_mobile_path",sample,m,c);
        if(sample<2) {
            r.add("random_spin_Pfaffian",sample,m,m.spinPathPfaffian(c),1.,1,&c);
            if(m.sites==2) r.add("random_sum_vs_original_6mode_Fock",sample,m,m.factorizedWordTrace(1.2,events),m.exactWordTrace(1.2,events),1,&c);
        }
    }
}
}
int main(int argc,char** argv) {
    try {
        int samples=8;std::uint64_t seed=160507421;std::string csv;
        for(int i=1;i<argc;++i) {
            std::string arg=argv[i];if(arg=="--help") {
                std::cout<<"sign_free_mobile_spin --samples N --seed N --csv PATH\n"
                    <<"Signed CT-INT configuration verifier for mobile auxiliary pairs coupled to fermions.\n"
                    <<"Includes full 3N-mode Pfaffians, original Fock oracles, negative and zero controls.\n"
                    <<"This executable is not a production Monte Carlo updater.\n";return 0;
            }
            require(i+1<argc,"Missing argument value");std::string value=argv[++i];
            if(arg=="--samples") samples=std::stoi(value);
            else if(arg=="--seed") seed=std::stoull(value);
            else if(arg=="--csv") csv=value;
            else throw std::runtime_error("Unknown argument: "+arg);
        }
        require(samples>0,"samples must be positive");mkl_set_num_threads(1);mkl_set_dynamic(0);
        std::ofstream file;if(!csv.empty()) {file.open(csv);require(file.good(),"Cannot open CSV");}
        Recorder r(csv.empty() ? std::cout:file,seed);std::mt19937_64 rng(seed);
        fixedWords(r);controls(r);randomWords(r,samples,rng);
        r.out.flush();require(r.out.good(),"Failed to write validation output");
        std::cerr<<std::setprecision(15)<<"PASS MobileSpin rows="<<r.rows<<" positive="<<r.positive<<" zero="<<r.zero
            <<" negative_controls="<<r.negative<<" max_scaled_error="<<r.maxError<<" max_positive_phase_error="<<r.maxPhase
            <<" samples_per_size="<<samples<<" seed="<<seed<<'\n';
    } catch(const std::exception& e) {std::cerr<<"FAIL MobileSpin: "<<e.what()<<'\n';return 1;}
    return 0;
}
