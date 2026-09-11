#include "catalog/frustrated_impurity.h"
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
Graph makeGraph(const FrustratedImpurity& m,bool dimer=true) {
    require(m.length<=4,"Validation graph is restricted to length<=4");int dim=1<<(2*m.cells());
    Graph g;g.adjacent.resize(dim);auto terms=m.terms(dimer);
    for(int v=0;v<dim;++v) {
        auto state=decode(v,m.cells());
        for(int e=0;e<int(terms.size());++e) {
            const auto& t=terms[e];int before=0;for(int i:t.cells) before=4*before+state[i];
            for(int after=0;after<t.matrix.rows();++after) if(after!=before && std::abs(t.matrix(after,before))>1e-12) {
                require(std::abs(t.matrix(after,before).imag())<1e-12,"Complex local impurity vertex");
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
ImpurityConfiguration configuration(const FrustratedImpurity& m,const Graph& g,const std::vector<int>& ids,
                                     const std::vector<int>& parities,double beta=1.2) {
    ImpurityConfiguration c;c.beta=beta;c.parities=parities;
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
        out<<"case,sample,seed,length,spin_sites,complex_modes,dimer_basis,weight_real,weight_imag,reference_real,reference_imag,sign_real,sign_imag,expected,scaled_error,beta,events,path,parities,chain,impurity,arm_imbalance,sign_free_parameters\n"<<std::setprecision(17);
    }
    void add(const std::string& name,int sample,const FrustratedImpurity& m,const ImpurityConfiguration& c,
             bool dimer,DataType w,DataType ref,int expected) {
        require(std::isfinite(std::abs(w)) && std::isfinite(std::abs(ref)),"Nonfinite impurity weight: "+name);
        double error=std::abs(w-ref)/std::max(1.,std::abs(ref));maxError=std::max(maxError,error);
        require(error<3e-9,"Impurity Pfaffian/reference disagreement: "+name);
        DataType phase=std::abs(w)>0. ? w/std::abs(w):DataType(std::numeric_limits<double>::quiet_NaN(),0.);
        if(expected==0) {require(std::abs(w)<1e-11 && std::abs(ref)<1e-11,"Nonzero zero control: "+name);++zero;}
        else {
            require(expected*ref.real()>0. && std::abs(ref.imag())<1e-10*std::max(1.,std::abs(ref)),"Incorrect reference sign: "+name);
            require(std::abs(w)>1e-14,"Unresolved small impurity signed weight: "+name);
            require(std::abs(phase-double(expected))<1e-8,"Incorrect measured impurity sign: "+name);
            if(expected>0) {++positive;maxPhase=std::max(maxPhase,std::abs(phase-1.));} else ++negative;
        }
        ++rows;out<<name<<','<<sample<<','<<seed<<','<<m.length<<','<<m.spinSites()<<','<<m.modes()<<','<<dimer<<','
            <<w.real()<<','<<w.imag()<<','<<ref.real()<<','<<ref.imag()<<','<<phase.real()<<','<<phase.imag()<<','
            <<expected<<','<<error<<','<<c.beta<<",\"";
        for(auto e:c.events) out<<e.time<<':'<<e.term<<';';out<<"\",\"";
        for(const auto& s:c.path) {for(int v:s) out<<v;out<<';';}out<<"\",\"";
        for(int p:c.parities) out<<p<<';';out<<"\","<<m.chain<<','<<m.impurity<<','<<m.armImbalance<<','<<m.signFreeParameters()<<'\n';
    }
    void check(const std::string& name,int sample,const FrustratedImpurity& m,const ImpurityConfiguration& c,
               bool dimer=true,int expected=1,bool fock=false) {
        add(name,sample,m,c,dimer,m.configurationPfaffian(c,dimer),
            fock ? m.configurationFockReference(c,dimer):m.configurationLocalReference(c,dimer),expected);
    }
};
void mustThrow(const std::function<void()>& fn,const std::string& what) {
    bool threw=false;try {fn();} catch(const std::exception&) {threw=true;}require(threw,"Missing guard: "+what);
}
void fixedAndControls(Recorder& r) {
    FrustratedImpurity m;auto graph=makeGraph(m);
    require(graph.negativeCycle.empty() && graph.negativeEdges>0 && graph.positiveEdges>0,"Dimer graph screening failed");
    for(int bits=0;bits<8;++bits) {
        std::vector<int> parities(3);for(int i=0;i<3;++i) parities[i]=bits&(1<<i) ? -1:1;
        auto c=configuration(m,graph,{0,20,38,20,0},parities);
        r.check("fixed_all_parity_sectors",bits,m,c);
        if(bits==0 || bits==7) r.check("fixed_full_9mode_Fock",bits,m,c,true,1,true);
    }
    // A non-retracing THREE-event positive loop exercises the turning signs.
    bool found=false;ImpurityConfiguration triangle;
    for(int x=0;x<int(graph.adjacent.size()) && !found;++x) for(auto a:graph.adjacent[x]) {
        for(auto b:graph.adjacent[a.to]) {
            auto back=std::find_if(graph.adjacent[b.to].begin(),graph.adjacent[b.to].end(),[&](auto e){return e.to==x;});
            if(back!=graph.adjacent[b.to].end()) {triangle=configuration(m,graph,{x,a.to,b.to,x},{1,-1,1});found=true;break;}
        }
        if(found) break;
    }
    require(found,"No non-retracing triangular configuration loop");r.check("positive_three_vertex_cycle",0,m,triangle,true,1,true);
    auto plain=makeGraph(m,false);require(!plain.negativeCycle.empty(),"Sz basis should have a sign problem");
    auto negative=configuration(m,plain,{2,16,32,2},{1,1,1});
    r.check("negative_original_Sz_triangle",0,m,negative,false,-1,true);
    auto unequal=m;unequal.armImbalance=.27;require(!unequal.signFreeParameters(),"Unequal arms claimed sign free");
    auto ug=makeGraph(unequal);require(!ug.negativeCycle.empty(),"Unequal arms did not create a frustrated cycle");
    auto un=configuration(unequal,ug,ug.negativeCycle,{1,-1,1});r.check("negative_unequal_arms",0,unequal,un,true,-1,true);
    // The nine-mode quartic Majorana H contains exactly eight spin copies.
    MatType hf=m.exactHamiltonian(),hs=m.exactSpinHamiltonian();
    require((hf-hf.adjoint()).norm()<1e-11,"Quartic Fock Hamiltonian is not Hermitian");
    Eigen::SelfAdjointEigenSolver<MatType> ef(hf),es(hs);require(ef.info()==Eigen::Success && es.info()==Eigen::Success,"Spectrum oracle failed");
    double spectrumError=0.;for(int i=0;i<512;++i) spectrumError=std::max(spectrumError,std::abs(ef.eigenvalues()(i)-es.eigenvalues()(i/8)));
    require(spectrumError<2e-11,"Full Fock space is not eight faithful spin copies");
    std::cerr<<"Impurity spectrum_error="<<std::setprecision(15)<<spectrumError<<" spin_ground_energy="<<es.eigenvalues()(0)<<'\n';
    auto c=configuration(m,graph,{0,20,38,20,0},{1,1,1});
    auto tiny=m;tiny.chain=1e-16;
    require(tiny.validConfiguration(c),"Tiny nonzero exchange was cut off as a forbidden vertex");
    auto tinyWeight=tiny.configurationPfaffian(c),tinyReference=tiny.configurationLocalReference(c);
    require(tinyReference.real()>0. && std::abs(tinyWeight/tinyReference-1.)<1e-9,"Tiny exchange signed weight mismatch");
    std::cerr<<"Impurity tiny_exchange_weight="<<tinyWeight<<" reference="<<tinyReference<<'\n';
    auto zero=[&](const std::string& name,const ImpurityConfiguration& v){r.check(name,0,m,v,true,0);};
    auto bad=c;bad.path.pop_back();zero("zero_path_length",bad);
    bad=c;bad.parities[0]=0;zero("zero_invalid_parity",bad);
    bad=c;bad.parities.pop_back();zero("zero_missing_parity",bad);
    bad=c;bad.path[1][0]=4;zero("zero_state_label",bad);
    bad=c;bad.path[1].pop_back();zero("zero_state_size",bad);
    bad=c;bad.path.back()[0]=1;zero("zero_nonclosed",bad);
    bad=c;bad.path[1]=bad.path[0];zero("zero_diagonal_vertex",bad);
    bad=c;bad.path[1][2]=1;zero("zero_transition_outside_support",bad);
    bad=c;bad.events[0].term=3;zero("zero_invalid_term",bad);
    bad=c;bad.events[0].time=-.1;zero("zero_negative_time",bad);
    bad=c;bad.events[0].time=1.;zero("zero_unsorted_time",bad);
    bad=c;bad.events.back().time=1.3;zero("zero_after_beta",bad);
    bad=c;bad.events[0].time=std::numeric_limits<double>::quiet_NaN();zero("zero_nonfinite_time",bad);
    bad=c;bad.beta=0.;zero("zero_nonpositive_beta",bad);
    auto absent=m;absent.chain=0.;r.check("zero_absent_chain_vertex",0,absent,c,true,0);
    ImpurityConfiguration empty{1.2,{}, {1,-1,1},{{0,1,2}}};r.check("empty_diagonal_word",0,m,empty,true,1,true);
    auto equal=c;for(auto& e:equal.events) e.time=.3;r.check("equal_time_ordered_vertices",0,m,equal);
    auto invalid=m;invalid.length=3;mustThrow([&]{invalid.validate();},"odd chain length");
    invalid=m;invalid.chain=-1.;mustThrow([&]{invalid.validate();},"negative exchange");
    invalid=m;invalid.armImbalance=-2.;mustThrow([&]{invalid.validate();},"negative first-arm exchange");
}
void randomPaths(Recorder& r,int samples,std::mt19937_64& rng) {
    std::uniform_real_distribution<double> coupling(.6,1.2),time(.04,1.16);
    for(int length:{2,4}) {
        FrustratedImpurity m;m.length=length;
        auto graph=makeGraph(m);require(graph.negativeCycle.empty(),"Finite dimer graph has a negative cycle");
        // This finite graph check is validation only. The proof for arbitrary
        // chain length is the source's paired-loop parity theorem.
        std::cerr<<"Impurity finite_graph length="<<length<<" states="<<graph.adjacent.size()
            <<" positive_directed_edges="<<graph.positiveEdges<<" negative_directed_edges="<<graph.negativeEdges<<'\n';
        for(int sample=0;sample<samples;++sample) {
            m.chain=coupling(rng);m.impurity=coupling(rng);
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
            r.check("random_closed_dimer_path",sample,m,c);
            if(length==2 && sample==0) r.check("random_full_Fock_path",sample,m,c,true,1,true);
        }
    }
}
}
int main(int argc,char** argv) {
    try {
        int samples=8;std::uint64_t seed=160805144;std::string csv;
        for(int i=1;i<argc;++i) {
            std::string arg=argv[i];if(arg=="--help") {
                std::cout<<"sign_free_frustrated_impurity --samples N --seed N --csv PATH\n"
                    <<"Signed CT-INT weight verifier for three equal-arm spin chains with a frustrated triangular impurity.\n"
                    <<"The all-fermion Hamiltonian is local Majorana quartic; no projection is required.\n"
                    <<"Includes Gaussian dimer projectors, full Fock spectra and actual negative controls.\n"
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
        std::cerr<<std::setprecision(15)<<"PASS FrustratedImpurity rows="<<r.rows<<" positive="<<r.positive<<" zero="<<r.zero
            <<" negative_controls="<<r.negative<<" max_scaled_error="<<r.maxError<<" max_positive_phase_error="<<r.maxPhase
            <<" samples_per_length="<<samples<<" seed="<<seed<<'\n';
    } catch(const std::exception& e) {std::cerr<<"FAIL FrustratedImpurity: "<<e.what()<<'\n';return 1;}
    return 0;
}
