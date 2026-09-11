#include "catalog/zigzag_models.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <queue>
#include <random>

using namespace catalog;
namespace {
struct Edge {int to,term;double coefficient;};
struct Graph {
    std::vector<std::vector<Edge>> adjacent;
    std::vector<int> phases,parents,negativeCycle;
    int positiveEdges=0,negativeEdges=0;
};
std::vector<int> decode(int v,int cells) {
    std::vector<int> out(cells);for(int i=cells-1;i>=0;--i) {out[i]=v%4;v/=4;}return out;
}
int encode(const std::vector<int>& state) {int out=0;for(int s:state) out=4*out+s;return out;}
Graph makeGraph(const ZigzagLadder& m,bool dimer=true) {
    require(m.rungs<=4,"Validation graph is restricted to rungs<=4");int dim=1<<(2*m.cells());
    Graph g;g.adjacent.resize(dim);auto terms=m.terms(dimer);
    for(int v=0;v<dim;++v) {
        auto state=decode(v,m.cells());
        for(int e=0;e<int(terms.size());++e) {
            const auto& t=terms[e];int before=0;for(int i:t.cells) before=4*before+state[i];
            for(int after=0;after<t.matrix.rows();++after) if(after!=before && std::abs(t.matrix(after,before))>1e-12) {
                require(std::abs(t.matrix(after,before).imag())<1e-12,"Complex local zigzag vertex");
                auto next=state;int value=after;for(int k=int(t.cells.size())-1;k>=0;--k) {next[t.cells[k]]=value%4;value/=4;}
                double coefficient=-t.matrix(after,before).real();g.adjacent[v].push_back({encode(next),e,coefficient});
                if(coefficient>0.) ++g.positiveEdges;else ++g.negativeEdges;
            }
        }
    }
    g.phases.assign(dim,0);g.parents.assign(dim,-1);std::queue<int> todo;
    for(int root=0;root<dim;++root) if(!g.phases[root]) {
        g.phases[root]=1;todo.push(root);
        while(!todo.empty()) {
            int x=todo.front();todo.pop();
            for(auto e:g.adjacent[x]) {
                int phase=g.phases[x]*(e.coefficient>0. ? 1:-1);
                if(!g.phases[e.to]) {g.phases[e.to]=phase;g.parents[e.to]=x;todo.push(e.to);}
                else if(g.phases[e.to]!=phase && g.negativeCycle.empty()) {
                    for(int v=x;v!=-1;v=g.parents[v]) g.negativeCycle.push_back(v);
                    std::reverse(g.negativeCycle.begin(),g.negativeCycle.end());g.negativeCycle.push_back(e.to);
                    for(int v=g.parents[e.to];v!=-1;v=g.parents[v]) g.negativeCycle.push_back(v);
                }
            }
        }
    }
    return g;
}
ZigzagConfiguration configuration(const ZigzagLadder& m,const Graph& g,const std::vector<int>& ids,
                                     const std::vector<int>& parities,double beta=.9) {
    ZigzagConfiguration c;c.beta=beta;c.parities=parities;
    for(int v:ids) c.path.push_back(decode(v,m.cells()));
    for(size_t k=0;k+1<ids.size();++k) {
        auto it=std::find_if(g.adjacent[ids[k]].begin(),g.adjacent[ids[k]].end(),[&](auto e){return e.to==ids[k+1];});
        require(it!=g.adjacent[ids[k]].end(),"Invalid validation graph walk");
        c.events.push_back({beta*double(k+1)/double(ids.size()),it->term});
    }
    return c;
}
struct Recorder {
    std::ostream& out;std::uint64_t seed;
    int rows=0,positive=0,negative=0,zero=0;double maxError=0.,maxPhase=0.;
    Recorder(std::ostream& output,std::uint64_t seed):out(output),seed(seed) {
        out<<"case,sample,seed,rungs,spin_sites,complex_modes,dimer_basis,weight_real,weight_imag,reference_real,reference_imag,sign_real,sign_imag,expected,scaled_error,beta,events,path,parities,J0,J1,J2,J3,sign_free_parameters\n"<<std::setprecision(17);
    }
    void add(const std::string& name,int sample,const ZigzagLadder& m,const ZigzagConfiguration& c,
             bool dimer,DataType w,DataType ref,int expected) {
        require(std::isfinite(std::abs(w)) && std::isfinite(std::abs(ref)),"Nonfinite zigzag weight: "+name);
        double error=std::abs(w-ref)/std::max(1.,std::abs(ref));maxError=std::max(maxError,error);
        require(error<3e-9,"Zigzag Pfaffian/reference disagreement: "+name);
        DataType phase=std::abs(w)>0. ? w/std::abs(w):DataType(std::numeric_limits<double>::quiet_NaN(),0.);
        if(expected==0) {require(std::abs(w)<1e-11 && std::abs(ref)<1e-11,"Nonzero zero control: "+name);++zero;}
        else {
            require(expected*ref.real()>0. && std::abs(ref.imag())<1e-10*std::max(1.,std::abs(ref)),"Incorrect reference sign: "+name);
            require(std::abs(w)>1e-14,"Unresolved small zigzag signed weight: "+name);
            require(std::abs(phase-double(expected))<1e-8,"Incorrect measured zigzag sign: "+name);
            if(expected>0) {++positive;maxPhase=std::max(maxPhase,std::abs(phase-1.));} else ++negative;
        }
        ++rows;out<<name<<','<<sample<<','<<seed<<','<<m.rungs<<','<<m.spinSites()<<','<<m.modes()<<','<<dimer<<','
            <<w.real()<<','<<w.imag()<<','<<ref.real()<<','<<ref.imag()<<','<<phase.real()<<','<<phase.imag()<<','
            <<expected<<','<<error<<','<<c.beta<<",\"";
        for(auto e:c.events) out<<e.time<<':'<<e.term<<';';out<<"\",\"";
        for(const auto& s:c.path) {for(int v:s) out<<v;out<<';';}out<<"\",\"";
        for(int p:c.parities) out<<p<<';';out<<"\","<<m.J0<<','<<m.J1<<','<<m.J2<<','<<m.J3<<','<<m.signFreeParameters()<<'\n';
    }
    void check(const std::string& name,int sample,const ZigzagLadder& m,const ZigzagConfiguration& c,
               bool dimer=true,int expected=1,bool fock=false) {
        add(name,sample,m,c,dimer,m.configurationPfaffian(c,dimer),
            fock ? m.configurationFockReference(c,dimer):m.configurationLocalReference(c,dimer),expected);
    }
};
void mustThrow(const std::function<void()>& fn,const std::string& what) {
    bool threw=false;try {fn();} catch(const std::exception&) {threw=true;}require(threw,"Missing guard: "+what);
}
void verifyKT(const ZigzagLadder& m,const Graph& graph) {
    require(m.signFreeParameters(),"KT test parameters outside theorem cone");
    auto terms=m.terms();ZigzagLadder pair;pair.rungs=2;
    for(int x=0;x<int(graph.adjacent.size());++x) for(auto edge:graph.adjacent[x]) {
        auto before=decode(x,m.cells()),after=decode(edge.to,m.cells());const auto& t=terms[edge.term];
        int full=m.ktPhase(before)*m.ktPhase(after);
        require(full*edge.coefficient>=-1e-12,"KT transformed Hamiltonian is not stoquastic");
        if(t.cells.size()==2) {
            int local=pair.ktPhase({before[t.cells[0]],before[t.cells[1]]})*pair.ktPhase({after[t.cells[0]],after[t.cells[1]]});
            require(full==local,"Global nonlocal KT transform disagrees with local bond transform");
        }
    }
    require(graph.negativeCycle.empty(),"Finite Bell graph has a negative cycle");
}
void fixedAndControls(Recorder& r) {
    ZigzagLadder m;auto graph=makeGraph(m);verifyKT(m,graph);
    for(int bits=0;bits<8;++bits) {
        std::vector<int> parities(3);for(int i=0;i<3;++i) parities[i]=bits&(1<<i) ? -1:1;
        auto c=configuration(m,graph,{0,5,10,0},parities);
        r.check("fixed_Bell_three_vertex_cycle",bits,m,c);
        if(bits==0 || bits==7) r.check("fixed_full_9mode_Fock",bits,m,c,true,1,true);
        if(bits==0) require(std::abs(m.configurationPfaffian(c)-.04879947523541948)<1e-12,"Independent Python fixed result mismatch");
    }
    auto outside=m;outside.J0=.8;outside.J1=.7;outside.J3=.2;
    require(!outside.signFreeParameters(),"Outside-cone parameters claimed sign free");
    auto og=makeGraph(outside);require(!og.negativeCycle.empty(),"Outside cone should have a negative cycle");
    auto negative=configuration(outside,og,{0,20,27,15,0},{1,-1,1});
    r.check("negative_outside_KT_cone",0,outside,negative,true,-1,true);
    require(std::abs(outside.configurationPfaffian(negative)+.0006974145661328646)<1e-12,"Independent negative Python result mismatch");
    auto plain=makeGraph(m,false);require(!plain.negativeCycle.empty(),"Original Sz basis should have a sign problem");
    auto sz=configuration(m,plain,plain.negativeCycle,{1,-1,1});r.check("negative_original_Sz_basis",0,m,sz,false,-1,true);
    auto c=configuration(m,graph,{0,5,10,0},{1,1,1});
    auto ferroRung=m;ferroRung.J2=-1.2;require(ferroRung.signFreeParameters(),"KT proof should allow negative rung exchange");
    r.check("ferromagnetic_rung_exchange",0,ferroRung,c,true,1,true);
    auto boundary=m;boundary.J3=boundary.J0+boundary.J1;require(boundary.signFreeParameters(),"Cone boundary rejected");
    verifyKT(boundary,makeGraph(boundary));r.check("KT_cone_boundary",0,boundary,c);
    auto vanish=configuration(m,graph,{0,15,0},{1,1,1});r.check("zero_boundary_d_vertex",0,boundary,vanish,true,0);
    // Full Fock space is eight faithful copies, not a projected spin subspace.
    MatType hf=m.exactHamiltonian(),hs=m.exactSpinHamiltonian();
    require((hf-hf.adjoint()).norm()<1e-11,"Quartic zigzag Hamiltonian is not Hermitian");
    Eigen::SelfAdjointEigenSolver<MatType> ef(hf),es(hs);require(ef.info()==Eigen::Success && es.info()==Eigen::Success,"Spectrum oracle failed");
    double spectrumError=0.;for(int i=0;i<512;++i) spectrumError=std::max(spectrumError,std::abs(ef.eigenvalues()(i)-es.eigenvalues()(i/8)));
    require(spectrumError<2e-11,"Full zigzag Fock spectrum not eight spin copies");
    std::cerr<<"Zigzag spectrum_error="<<std::setprecision(15)<<spectrumError<<" spin_ground_energy="<<es.eigenvalues()(0)<<'\n';
    // Never use a fixed absolute vertex cutoff to discard small physical couplings.
    auto tiny=m;tiny.J0*=1e-16;tiny.J1*=1e-16;tiny.J3*=1e-16;
    require(tiny.validConfiguration(c),"Tiny physical coupling was cut off");
    auto tw=tiny.configurationPfaffian(c),tr=tiny.configurationLocalReference(c);
    require(tr.real()>0. && std::abs(tw/tr-1.)<1e-9,"Tiny-coupling signed weight mismatch");
    std::cerr<<"Zigzag tiny_exchange_weight="<<tw<<" reference="<<tr<<'\n';
    auto zero=[&](const std::string& name,const ZigzagConfiguration& v){r.check(name,0,m,v,true,0);};
    auto bad=c;bad.path.pop_back();zero("zero_path_length",bad);
    bad=c;bad.parities[0]=0;zero("zero_invalid_parity",bad);
    bad=c;bad.parities.pop_back();zero("zero_missing_parity",bad);
    bad=c;bad.path[1][0]=4;zero("zero_state_label",bad);
    bad=c;bad.path[1].pop_back();zero("zero_state_size",bad);
    bad=c;bad.path.back()[0]=1;zero("zero_nonclosed",bad);
    bad=c;bad.path[1]=bad.path[0];zero("zero_diagonal_vertex",bad);
    bad=c;bad.path[1][0]=1;zero("zero_transition_outside_support",bad);
    bad=c;bad.events[0].term=5;zero("zero_invalid_term",bad);
    bad=c;bad.events[0].time=-.1;zero("zero_negative_time",bad);
    bad=c;bad.events[0].time=.8;zero("zero_unsorted_time",bad);
    bad=c;bad.events.back().time=1.;zero("zero_after_beta",bad);
    bad=c;bad.events[0].time=std::numeric_limits<double>::quiet_NaN();zero("zero_nonfinite_time",bad);
    bad=c;bad.beta=0.;zero("zero_nonpositive_beta",bad);
    auto absent=m;absent.J0=absent.J1=absent.J3=0.;r.check("zero_absent_interrung_vertex",0,absent,c,true,0);
    ZigzagConfiguration empty{.9,{}, {1,-1,1},{{0,1,2}}};r.check("empty_diagonal_word",0,m,empty,true,1,true);
    auto equal=c;for(auto& e:equal.events) e.time=.3;r.check("equal_time_ordered_vertices",0,m,equal);
    auto invalid=m;invalid.rungs=1;mustThrow([&]{invalid.validate();},"too few rungs");
    invalid=m;invalid.J0=-1.;mustThrow([&]{invalid.validate();},"negative interrung exchange");
    invalid=m;invalid.J2=std::numeric_limits<double>::quiet_NaN();mustThrow([&]{invalid.validate();},"nonfinite rung exchange");
}
void randomPaths(Recorder& r,int samples,std::mt19937_64& rng) {
    std::uniform_real_distribution<double> coupling(.15,.45),time(.04,.86);
    for(int rungs:{2,3,4}) {
        ZigzagLadder m;m.rungs=rungs;auto graph=makeGraph(m);verifyKT(m,graph);
        std::cerr<<"Zigzag finite_KT_check rungs="<<rungs<<" states="<<graph.adjacent.size()
            <<" positive_directed_edges="<<graph.positiveEdges<<" negative_directed_edges="<<graph.negativeEdges<<'\n';
        for(int sample=0;sample<samples;++sample) {
            m.J0=coupling(rng);m.J1=coupling(rng);m.J3=m.J0+m.J1+coupling(rng);
            m.J2=(rng()%2 ? 1.:-1.)*(.5+coupling(rng));
            int start;do {start=rng()%graph.adjacent.size();} while(graph.adjacent[start].empty());
            std::vector<int> parents(graph.adjacent.size(),-2);parents[start]=-1;std::queue<int> todo;todo.push(start);
            while(!todo.empty()) {int x=todo.front();todo.pop();for(auto e:graph.adjacent[x]) if(parents[e.to]==-2) {parents[e.to]=x;todo.push(e.to);}}
            std::vector<int> walk{start};int steps=3+rng()%4;
            for(int k=0;k<steps;++k) {const auto& options=graph.adjacent[walk.back()];walk.push_back(options[rng()%options.size()].to);}
            for(int x=parents[walk.back()];x!=-1;x=parents[x]) walk.push_back(x);
            require(walk.front()==walk.back(),"Random walk did not close");
            std::vector<int> parity(m.cells());for(int& p:parity) p=rng()%2 ? 1:-1;
            auto c=configuration(m,graph,walk,parity);std::vector<double> times(c.events.size());
            for(double& t:times) t=time(rng);std::sort(times.begin(),times.end());
            for(size_t k=0;k<times.size();++k) c.events[k].time=times[k];
            r.check("random_closed_Bell_path",sample,m,c);
            if(rungs<=3 && sample==0) r.check("random_full_Fock_path",sample,m,c,true,1,true);
        }
    }
}
}
int main(int argc,char** argv) {
    try {
        int samples=8;std::uint64_t seed=9707019;std::string csv;
        for(int i=1;i<argc;++i) {
            std::string arg=argv[i];if(arg=="--help") {
                std::cout<<"sign_free_zigzag --samples N --seed N --csv PATH\n"
                    <<"Signed CT-INT weight verifier for the open J0-J1-J2-J3 zigzag ladder.\n"
                    <<"Bell-rung Gaussian states preserve all original vertex signs. J3>=J0+J1.\n"
                    <<"The all-fermion Hamiltonian is local Majorana quartic, with no hard projection.\n"
                    <<"This executable is not a production Monte Carlo updater.\n";return 0;
            }
            require(i+1<argc,"Missing argument value");std::string value=argv[++i];
            if(arg=="--samples") samples=std::stoi(value);else if(arg=="--seed") seed=std::stoull(value);
            else if(arg=="--csv") csv=value;else throw std::runtime_error("Unknown argument: "+arg);
        }
        require(samples>0,"samples must be positive");mkl_set_num_threads(1);mkl_set_dynamic(0);
        std::ofstream file;if(!csv.empty()) {file.open(csv);require(file.good(),"Cannot open CSV");}
        Recorder r(csv.empty() ? std::cout:file,seed);std::mt19937_64 rng(seed);fixedAndControls(r);randomPaths(r,samples,rng);
        r.out.flush();require(r.out.good(),"Failed to write CSV");
        std::cerr<<std::setprecision(15)<<"PASS ZigzagLadder rows="<<r.rows<<" positive="<<r.positive<<" zero="<<r.zero
            <<" negative_controls="<<r.negative<<" max_scaled_error="<<r.maxError<<" max_positive_phase_error="<<r.maxPhase
            <<" samples_per_rungs="<<samples<<" seed="<<seed<<'\n';
    } catch(const std::exception& e) {std::cerr<<"FAIL ZigzagLadder: "<<e.what()<<'\n';return 1;}
    return 0;
}
