#include "catalog/bilayer_models.h"
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
Graph makeGraph(const AnisotropicBilayer& m,bool dimer=true) {
    require(m.sites<=4,"Validation graph is restricted to sites<=4");int dim=1<<(2*m.cells());
    Graph g;g.adjacent.resize(dim);auto terms=m.terms(dimer);
    for(int v=0;v<dim;++v) {
        auto state=decode(v,m.cells());
        for(int e=0;e<int(terms.size());++e) {
            const auto& t=terms[e];int before=0;for(int i:t.cells) before=4*before+state[i];
            for(int after=0;after<t.matrix.rows();++after) if(after!=before && std::abs(t.matrix(after,before))>1e-12) {
                require(std::abs(t.matrix(after,before).imag())<1e-12,"Complex local bilayer vertex");
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
BilayerConfiguration configuration(const AnisotropicBilayer& m,const Graph& g,const std::vector<int>& ids,
                                     const std::vector<int>& parities,double beta=1.2) {
    BilayerConfiguration c;c.beta=beta;c.parities=parities;
    for(int v:ids) c.path.push_back(decode(v,m.cells()));
    for(size_t k=0;k+1<ids.size();++k) {
        auto it=std::find_if(g.adjacent[ids[k]].begin(),g.adjacent[ids[k]].end(),[&](auto e){return e.to==ids[k+1];});
        require(it!=g.adjacent[ids[k]].end(),"Invalid validation graph walk");
        c.events.push_back({beta*double(k+1)/double(ids.size()),it->term});
    }
    return c;
}

struct Recorder {
    std::ostream& out;int rows=0,pos=0,neg=0,zero=0;double error=0.,phaseError=0.;
    Recorder(std::ostream& o):out(o){out<<"case,sample,sites,Jz,Jperp,Kz,Kperp,Dz,Dperp,real_pf,imag_pf,real_reference,imag_reference,relative_error,expected,path,parities,events\n"<<std::setprecision(17);}
    void check(std::string name,int sample,const AnisotropicBilayer& m,const BilayerConfiguration& c,int expected=1,bool fock=false,bool dimer=true) {
        DataType w=m.configurationPfaffian(c,dimer),r=fock?m.configurationFockReference(c,dimer):m.configurationLocalReference(c,dimer);
        double e=std::abs(w-r)/std::max(1.,std::abs(r));error=std::max(error,e);
        require(std::isfinite(e)&&e<3e-9,"Signed bilayer Pf/reference mismatch: "+name);
        if(expected==0){require(w==DataType(0.)&&r==DataType(0.),"Nonzero forbidden configuration");++zero;}
        else {
            require(expected*r.real()>0&&std::abs(r.imag())<1e-10*std::abs(r),"Incorrect bilayer reference sign: "+name);
            require(std::abs(w)>0.,"Unresolved bilayer weight");double ph=std::abs(w/std::abs(w)-double(expected));
            require(ph<1e-8,"Incorrect measured bilayer sign: "+name);phaseError=std::max(phaseError,ph);
            if(expected>0)++pos;else ++neg;
        }
        ++rows;out<<name<<','<<sample<<','<<m.sites<<','<<m.Jz<<','<<m.Jperp<<','<<m.Kz<<','<<m.Kperp<<','<<m.Dz<<','<<m.Dperp<<','
            <<w.real()<<','<<w.imag()<<','<<r.real()<<','<<r.imag()<<','<<e<<','<<expected<<",\"";
        for(auto s:c.path){for(int v:s)out<<v;out<<';';}out<<"\",\"";for(int p:c.parities)out<<p<<';';out<<"\",\"";
        for(auto e:c.events)out<<e.time<<':'<<e.term<<';';out<<"\"\n";
    }
};
std::vector<int> randomLoop(const Graph& graph,std::mt19937_64& rng) {
    int start;do{start=rng()%graph.adjacent.size();}while(graph.adjacent[start].empty());
    std::vector<int> parents(graph.adjacent.size(),-2);parents[start]=-1;std::queue<int> todo;todo.push(start);
    while(!todo.empty()){int x=todo.front();todo.pop();for(auto e:graph.adjacent[x])if(parents[e.to]==-2){parents[e.to]=x;todo.push(e.to);}}
    std::vector<int> walk{start};for(int k=0,n=3+rng()%3;k<n;++k){auto& options=graph.adjacent[walk.back()];walk.push_back(options[rng()%options.size()].to);}
    for(int x=parents[walk.back()];x!=-1;x=parents[x])walk.push_back(x);return walk;
}
void algebra(const AnisotropicBilayer& m) {
    auto two=m;two.sites=2;two.edges={{0,1}};
    MatType h=two.exactHamiltonian(),spin=two.exactSpinHamiltonian();
    Eigen::SelfAdjointEigenSolver<MatType> ef(h),es(spin);
    double e=0.;for(int k=0;k<64;++k)e=std::max(e,std::abs(ef.eigenvalues()(k)-es.eigenvalues()(k/4)));
    require(e<1e-11,"Full Fock spectrum is not four spin copies");
    MatType sz=two.exactSpinHamiltonian(false);Eigen::SelfAdjointEigenSolver<MatType> ez(sz);
    require((es.eigenvalues()-ez.eigenvalues()).norm()<1e-11,"Dimer basis changed physical spectrum");
    // The diagonal and every off-diagonal CT vertex reconstruct physical H.
    MatType d=spin.diagonal().asDiagonal();MatType k=d-spin;
    double target=0.;for(auto x:es.eigenvalues())target+=std::exp(-.9*x);
    double previous=1e99;
    for(int slices:{32,64,128}) {
        MatType step=(MatType::Identity(16,16)+(.9/slices)*k)*matrixExp(-(.9/slices)*d),prod=MatType::Identity(16,16);
        for(int i=0;i<slices;++i)prod=step*prod;
        double err=std::abs(prod.trace().real()-target);require(err<previous,"CT Trotter completeness failed to converge");previous=err;
    }
    // Unlike the conserved-sector subfamily, every tested face has a truly
    // dynamic singlet label on a cell: [H,P_s(cell 0)] is nonzero.
    MatType singlet=MatType::Zero(16,16);for(int i=0;i<4;++i)singlet(i,i)=1.;
    require((spin*singlet-singlet*spin).norm()>.1,"Only a static dimer sector was tested");
    std::cerr<<"bilayer face="<<(m.Kz==0?0:m.Kperp==0?1:2)<<" full_Fock_copy_error="<<e
        <<" dynamic_singlet_commutator="<<(spin*singlet-singlet*spin).norm()<<'\n';
}
void fixed(Recorder& r) {
    for(int face=0;face<3;++face) {
        auto m=squareBilayer(face);require(m.signFreeParameters(),"Rejected valid bilayer face");algebra(m);
        auto graph=makeGraph(m);require(graph.negativeCycle.empty(),"Negative loop on a sign-free face");
        int from=0;while(graph.adjacent[from].empty())++from;int to=graph.adjacent[from][0].to;
        for(int bits=0;bits<16;++bits) {
            std::vector<int> parity(4);for(int i=0;i<4;++i)parity[i]=bits>>i&1?-1:1;
            auto c=configuration(m,graph,{from,to,from},parity,.9);
            r.check("all_full_Fock_parities_face_"+std::to_string(face),bits,m,c);
            if(bits==0||bits==15)r.check("full_12mode_occupation_face_"+std::to_string(face),bits,m,c,1,true);
        }
        auto two=m;two.sites=2;two.edges={{0,1}};auto tg=makeGraph(two);int index=0;
        for(int x=0;x<int(tg.adjacent.size());++x)for(auto e:tg.adjacent[x]) {
            auto c=configuration(two,tg,{x,e.to,x},{1,-1},.9);
            r.check("all_two_event_vertices_face_"+std::to_string(face),index++,two,c);
        }
        // Every physical parity sector, including all empty-time paths.
        for(int bits=0;bits<4;++bits)for(int state=0;state<16;++state) {
            auto c=configuration(two,tg,{state},{bits&1?-1:1,bits&2?-1:1},.9);
            r.check("all_diagonal_states_face_"+std::to_string(face),16*bits+state,two,c);
        }
    }
    auto bad=squareBilayer(3);require(!bad.signFreeParameters(),"All-nonzero couplings incorrectly certified");
    auto bg=makeGraph(bad);require(!bg.negativeCycle.empty(),"Square negative cycle missing");
    auto c=configuration(bad,bg,{1,4,135,39,45,1},{1,-1,1,-1},.9);
    r.check("negative_all_nonzero_square",0,bad,c,-1);r.check("negative_full_12mode_occupation",0,bad,c,-1,true);
    auto two=bad;two.sites=2;two.edges={{0,1}};
    require(makeGraph(two).negativeCycle.empty(),"Two-cell caveat changed unexpectedly");
    auto m=squareBilayer();auto plain=makeGraph(m,false);
    require(!plain.negativeCycle.empty(),"Original Sz basis negative loop missing");
    auto pn=configuration(m,plain,plain.negativeCycle,{1,1,1,1},.9);r.check("negative_original_Sz_basis",0,m,pn,-1,true,false);
    auto graph=makeGraph(m);int from=0;while(graph.adjacent[from].empty())++from;
    auto good=configuration(m,graph,{from,graph.adjacent[from][0].to,from},{1,-1,1,-1},.9);
    auto zero=good;zero.path.back()[0]=(zero.path.back()[0]+1)%4;r.check("zero_open_path",0,m,zero,0);
    zero=good;zero.parities[0]=0;r.check("zero_invalid_parity",0,m,zero,0);
    zero=good;zero.path[1]=zero.path[0];r.check("zero_diagonal_vertex",0,m,zero,0);
    zero=good;zero.events[0].time=-.1;r.check("zero_invalid_time",0,m,zero,0);
    auto disconnected=m;disconnected.Jperp=disconnected.Kperp=disconnected.Kz=0.;r.check("zero_disabled_vertex",0,disconnected,good,0);
    auto nonbip=m;nonbip.sites=3;nonbip.edges={{0,1},{1,2},{2,0}};require(!nonbip.signFreeParameters(),"Nonbipartite graph certified");
    auto full=two;full.Kz=full.Kperp=0.;MatType hs=full.exactSpinHamiltonian(),p=MatType::Zero(16,16);for(int i=0;i<4;++i)p(i,i)=1.;
    require((hs*p-p*hs).norm()<1e-12,"Fully frustrated conservation check failed");
    // A genuine zero signed Pf contraction, independent of configuration guards:
    // orthogonal normalized dimer/parity states have Tr(P_s P_t)=0.
    auto ps=two.stateProjector({0,0},{1,-1}),pt=two.stateProjector({1,0},{1,-1});
    auto actualZero=traceGaussianFactors({ps,pt});
    require(std::abs(actualZero)<1e-12,"Orthogonal Gaussian projectors have a nonzero trace");
    ++r.rows;++r.zero;r.error=std::max(r.error,std::abs(actualZero));
    r.out<<"zero_orthogonal_projector_Pf,0,"<<two.sites<<','<<two.Jz<<','<<two.Jperp<<','<<two.Kz<<','<<two.Kperp<<','<<two.Dz<<','<<two.Dperp
        <<','<<actualZero.real()<<','<<actualZero.imag()<<",0,0,"<<std::abs(actualZero)<<",0,\"00;10\",\"1;-1\",\"\"\n";
}
void randomTests(Recorder& r,int samples,std::mt19937_64& rng) {
    std::uniform_real_distribution<double> amplitude(.6,1.4),time(.02,.88);
    for(int face=0;face<3;++face)for(int signs=0;signs<4;++signs) {
        auto m=squareBilayer(face);
        double* first=face==2?&m.Kz:&m.Jperp;double* second=face==1?&m.Kz:&m.Kperp;
        if(signs&1)*first *= -1;if(signs&2)*second *= -1;
        m.Dz*=signs&1?-1.:1.;m.Dperp*=signs&2?-1.:1.;
        auto graph=makeGraph(m);require(graph.negativeCycle.empty(),"Sign-independent parameter face failed");
        for(int sample=0;sample<samples;++sample) {
            auto varied=m;varied.Jz*=amplitude(rng);varied.Jperp*=amplitude(rng);varied.Kz*=amplitude(rng);varied.Kperp*=amplitude(rng);
            std::vector<int> parity(4);for(int& p:parity)p=rng()%2?1:-1;
            auto c=configuration(varied,graph,randomLoop(graph,rng),parity,.9);
            std::vector<double> times(c.events.size());for(double& t:times)t=time(rng);std::sort(times.begin(),times.end());
            for(size_t i=0;i<times.size();++i)c.events[i].time=times[i];
            r.check("random_face_"+std::to_string(face)+"_signs_"+std::to_string(signs),sample,varied,c);
            if(sample==0&&signs==0)r.check("random_full_12mode_face_"+std::to_string(face),sample,varied,c,1,true);
        }
    }
    // A connected non-square bipartite graph, with nonuniform coordination.
    for(int face=0;face<3;++face) {
        auto m=squareBilayer(face);m.edges={{0,1},{0,2},{0,3}};auto graph=makeGraph(m);
        require(graph.negativeCycle.empty()&&m.signFreeParameters(),"Bipartite tree check failed");
        for(int sample=0;sample<samples;++sample) {
            auto c=configuration(m,graph,randomLoop(graph,rng),{1,-1,1,-1},.9);r.check("star_graph_face_"+std::to_string(face),sample,m,c);
        }
    }
}
}
int main(int argc,char** argv) {try {
    int samples=16;std::uint64_t seed=731;std::string csv;
    for(int i=1;i<argc;++i){std::string arg=argv[i];require(i+1<argc,"Missing argument value");std::string value=argv[++i];
        if(arg=="--samples")samples=std::stoi(value);else if(arg=="--seed")seed=std::stoull(value);else if(arg=="--csv")csv=value;else throw std::runtime_error("Unknown argument");}
    require(samples>0,"samples must be positive");mkl_set_num_threads(1);mkl_set_dynamic(0);
    std::ofstream file;if(!csv.empty()){file.open(csv);require(file.good(),"Cannot open CSV");}
    Recorder r(csv.empty()?std::cout:file);std::mt19937_64 rng(seed);fixed(r);randomTests(r,samples,rng);
    r.out.flush();require(r.out.good(),"Failed CSV write");
    std::cerr<<std::setprecision(17)<<"PASS Bilayer rows="<<r.rows<<" positive="<<r.pos<<" zero="<<r.zero<<" negative_controls="<<r.neg
       <<" max_scaled_error="<<r.error<<" max_phase_error="<<r.phaseError<<" samples_per_face_sign="<<samples<<" seed="<<seed<<'\n';
    return 0;
}catch(const std::exception& e){std::cerr<<"FAIL Bilayer: "<<e.what()<<'\n';return 1;}}
