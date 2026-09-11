#include "catalog/albert_models.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>

using namespace catalog;
namespace {
MatType eye(int n){return MatType::Identity(n,n);}
void close(DataType actual,DataType ref,const std::string& label,double tol=2e-7) {
    require(std::isfinite(std::abs(actual)) && std::isfinite(std::abs(ref)) &&
            std::abs(actual-ref)<=tol*std::max(1.,std::abs(ref)),label);
}
struct Report {
    std::ostream& out;long rows=0,positive=0,negative=0;double error=0.,phaseError=0.;
    Report(std::ostream& stream):out(stream) {
        out<<std::setprecision(17)<<"model,kind,index,modes,weight_real,weight_imag,reference_real,reference_imag,relative_error,sign_real,sign_imag,attempts,accepted\n";
    }
    void add(const std::string& name,const std::string& kind,int index,int modes,
             DataType z,DataType ref,int expected=1,long attempts=0,long accepted=0) {
        close(z,ref,name+" "+kind+" oracle mismatch");
        require(std::abs(z)>1e-12,"Unresolved Albert zero is not sign=1");
        DataType phase=z/std::abs(z);close(phase,double(expected),name+" wrong phase");
        double err=std::abs(z-ref)/std::max(1.,std::abs(ref));error=std::max(error,err);
        if(expected>0){++positive;phaseError=std::max(phaseError,std::abs(phase-1.));}else ++negative;
        ++rows;out<<name<<','<<kind<<','<<index<<','<<modes<<','<<z.real()<<','<<z.imag()<<','
            <<ref.real()<<','<<ref.imag()<<','<<err<<','<<phase.real()<<','<<phase.imag()<<','<<attempts<<','<<accepted<<'\n';
    }
};
Albert randomAlbert(rdGenerator& rng,double scale=.2) {
    Albert k{};for(auto& t:albertBasis())k=albertAdd(k,albertScale(t,scale*rng.rdNormal()));return k;
}
void controls(Report& report,rdGenerator& rng) {
    auto& basis=albertBasis();auto& directions=albertTracelessBasis();
    require(basis.size()==27 && directions.size()==26,"Full Albert dimension");
    Albert unit{};unit[0][0]=unit[4][0]=unit[8][0]=1.;
    require((albertAction(unit)-eye(27)).norm()<1e-12,"Albert unit action");
    for(int a=0;a<27;++a)for(int b=0;b<27;++b)
        require(std::abs(albertInner(basis[a],basis[b])-(a==b))<1e-12,"Albert orbital Gram matrix");
    for(int a=0;a<26;++a) {
        auto h=albertAction(directions[a]);
        require(std::abs(albertTrace(directions[a]))<1e-12 && std::abs(h.trace())<1e-12,"Albert traceless direction");
        require((h-h.adjoint()).norm()<1e-12,"Albert Hermitian orbital Hamiltonian");
        for(int b=0;b<26;++b)require(std::abs(albertInner(directions[a],directions[b])-(a==b))<1e-12,"Full 26 independent directions");
    }
    auto a=albertAction(directions[0]),b=albertAction(directions[2]);
    require((a*b-b*a).norm()>.1,"Albert family must be noncommuting");
    // Odd real eigenspaces rule out the common Wei contraction conditions for +/-D.
    Eigen::SelfAdjointEigenSolver<MatType> es(numberMajorana(eye(27)));int multiplicity=0;
    for(auto v:es.eigenvalues())if(std::abs(v-1.)<1e-12)++multiplicity;
    require(multiplicity==27,"Albert odd Majorana eigenspace");
    double invariantError=0.,actionError=0.;
    for(int test=0;test<24;++test) {
        auto x=randomAlbert(rng),y=randomAlbert(rng),k=randomAlbert(rng);
        k=albertAdd(k,albertScale(unit,-albertTrace(k)/3.));
        auto xx=albertJordan(x,x);
        auto difference=albertAdd(albertJordan(x,albertJordan(xx,y)),albertScale(albertJordan(xx,albertJordan(x,y)),-1.));
        require(std::sqrt(albertInner(difference,difference))<1e-10,"Albert Jordan identity");
        // The quadratic-representation oracle does not use the 27x27 matrix exponential.
        MatType oracle=albertExponentialAction(k),dense=matrixExp(albertAction(k));
        actionError=std::max(actionError,(oracle-dense).norm()/dense.norm());
        require(actionError<1e-11,"Independent Jordan exponential action");
        auto mapped=albertApply(oracle,x);double before=albertDeterminant(x),after=albertDeterminant(mapped);
        invariantError=std::max(invariantError,std::abs(before-after)/std::max(1.,std::abs(before)));
        require(invariantError<1e-10,"Albert cubic determinant invariant");
    }
    // Full real 3x3 subgroup: Sym^2(R^3) + 7 Lambda^2(R^3).
    // The two-dimensional commutator has x=tr(g); the full 27-orbital trace is
    // 2^9 x^2 (x+2)^8. This is an independent character, not |Pfaffian|.
    Albert x{},z{};x[1][0]=x[3][0]=1.;z[0][0]=1.;z[4][0]=-1.;
    for(int test=0;test<2;++test) {
        double t=std::log(test?3.:1.5),u=std::sinh(t),trace=2.-4.*u*u*u*u;
        auto A=numberMajorana(albertAction(albertScale(x,2.*t)));
        auto B=numberMajorana(albertAction(albertScale(z,2.*t)));
        std::vector<MatType> pulses={A,B,-A,-B};
        double ref=512.*trace*trace*std::pow(trace+2.,8);
        report.add("albert27","real_subgroup_character",test,27,traceProduct(pulses),ref);
    }
    // Removing the common-field restriction: one physical orbital filter inside
    // a two-orbital invariant block gives an exact negative configuration.
    double t=std::log(3.);MatType filter=MatType::Zero(27,27);filter(11,11)=1.;
    MatType A=numberMajorana(albertAction(albertScale(x,2.*t))),B=2.*t*numberMajorana(filter);
    std::vector<MatType> bad={A,B,-A,-B};
    MatType hop=MatType::Zero(2,2),leaf=hop;hop(0,1)=hop(1,0)=1.;leaf(0,0)=1.;
    std::vector<MatType> small={t*numberMajorana(hop),2.*t*numberMajorana(leaf),-t*numberMajorana(hop),-2.*t*numberMajorana(leaf)};
    auto negativeRef=std::ldexp(1.,25)*exactTrace(small);
    close(negativeRef,-std::ldexp(700.,25)/81.,"Exact two-mode Fock negative certificate");
    report.add("albert_independent_orbital","negative",0,27,traceProduct(bad),negativeRef,-1);
    std::cerr<<"PASS Albert controls action_error="<<actionError<<" cubic_error="<<invariantError<<'\n';
}
void verify(Report& report,rdGenerator& rng,int cells,int samples,int sweeps) {
    auto model=makeAlbert(cells);auto& h=model.h;const int slices=2;const double beta=.45;
    require(int(h.channels.size())==27*cells+2*(cells-1),"Albert full HS directions and local bonds");
    for(int sample=0;sample<samples;++sample) {
        std::vector<double> fields(h.fieldCount(slices));for(auto& f:fields)f=(sample%2?.7:1.5)*rng.rdNormal();
        auto pulses=h.pulses(slices,beta,fields);auto z=traceProduct(pulses);
        close(z,determinantTrace(pulses),"Albert independent number determinant");
        report.add(h.id,"albert_fields",sample,h.modes,z,model.reference(slices,beta,fields));
    }
    Walker walker(h,slices,beta,&rng);PfQMC q(&walker,5);
    for(int sweep=0;sweep<=sweeps;++sweep) {
        auto fields=walker.fields();auto pulses=h.pulses(slices,beta,fields);auto z=traceProduct(pulses);
        close(q.sign,z/std::abs(z),"Albert actual PFQMC tracked sign");
        close(q.getSignRaw(),q.sign,"Albert raw signed Pfaffian");
        MatType B=eye(q.nDim);for(auto& pulse:pulses)B=(matrixExp(pulse)*B).eval();
        MatType greens=2.*(eye(q.nDim)+B).inverse();
        require((greens-q.g).norm()<2e-7*std::max(1.,greens.norm()),"Albert actual sweep Greens");
        report.add(h.id,"albert_sweep",sweep,h.modes,z,model.reference(slices,beta,fields),1,walker.attempts(),walker.accepted());
        if(sweep<sweeps){if(sweep%2==0)q.rightSweep();else q.leftSweep();}
    }
    require(walker.attempts()==sweeps*h.fieldCount(slices) && walker.accepted()>0,"Albert walker never moved");
    std::cerr<<"PASS "<<h.id<<" channels="<<h.channels.size()<<'\n';
}
}
int main(int argc,char** argv) {
    try {
        mkl_set_num_threads(1);mkl_set_dynamic(0);int samples=2,sweeps=2,seed=9120401;std::string model="all",csv;
        for(int i=1;i<argc;++i) {
            std::string key=argv[i];require(i+1<argc,"Missing Albert option value");std::string value=argv[++i];
            if(key=="--samples")samples=std::stoi(value);else if(key=="--sweeps")sweeps=std::stoi(value);
            else if(key=="--seed")seed=std::stoi(value);else if(key=="--model")model=value;else if(key=="--csv")csv=value;
            else throw std::runtime_error("Unknown Albert option "+key);
        }
        require(samples>0 && sweeps>=2,"Need samples>0 and both sweep directions");
        require(model=="all" || model=="controls" || model=="albert27" || model=="albert27_lattice","Unknown Albert model");
        std::ofstream file;if(!csv.empty()){file.open(csv);require(file.good(),"Cannot open Albert CSV");}
        Report report(csv.empty()?std::cout:file);rdGenerator rng(seed);
        if(model=="all" || model=="controls")controls(report,rng);
        if(model=="all" || model=="albert27")verify(report,rng,1,samples,sweeps);
        if(model=="all" || model=="albert27_lattice")verify(report,rng,2,samples,sweeps);
        std::cerr<<"PASS rows="<<report.rows<<" positive="<<report.positive<<" negative="<<report.negative
            <<" max_relative_error="<<report.error<<" max_phase_error="<<report.phaseError<<" seed="<<seed<<'\n';
    }catch(const std::exception& error){std::cerr<<"FAIL: "<<error.what()<<'\n';return 1;}return 0;
}
