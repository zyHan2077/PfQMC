#include "catalog/poisson_models.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <random>
#include <sstream>
using namespace catalog;
namespace {
struct Records {
    std::ostream& out;int rows=0,positive=0,negative=0,zeros=0;double error=0,phase=0;
    explicit Records(std::ostream& o):out(o){out<<"cells,kind,sample,word,real,imag,reference_real,reference_imag,expected_sign,relative_error,vertex_parameters\n"<<std::setprecision(17);}
    void add(int cells,const std::string& kind,int sample,const std::vector<int>& word,DataType z,DataType ref,int expected=1,const std::string& parameters="") {
        double e=std::abs(z-ref)/(expected==0?1.:std::max(1e-290,std::abs(ref)));error=std::max(error,e);
        require(e<2e-7,"Poisson signed Pfaffian/polynomial Fock mismatch");
        if(expected==0){require(std::abs(z)<1e-8,"Poisson expected zero");++zeros;}
        else {require(std::abs(z)>0 && std::isfinite(std::abs(z)),"Poisson invalid nonzero weight");double p=std::abs(z/std::abs(z)-double(expected));require(p<2e-7,"Poisson wrong sign");if(expected==1){++positive;phase=std::max(phase,p);}else ++negative;}
        out<<cells<<','<<kind<<','<<sample<<',';for(int i=0;i<int(word.size());++i){if(i)out<<':';out<<word[i];}
        out<<','<<z.real()<<','<<z.imag()<<','<<ref.real()<<','<<ref.imag()<<','<<expected<<','<<e<<','<<parameters<<'\n';++rows;
    }
};
void physicalChecks(Records& rec) {
    auto model=poissonModel(2);auto h=model.exactHamiltonian();MatType gaussianH=MatType::Zero(64,64);
    for(int v=0;v<int(model.vertices.size());++v) {
        auto g=matrixExp(fockMatrix(model.generator(v)));
        require((g-model.polynomialVertex(v)).norm()<1e-11,"Poisson finite polynomial/Gaussian vertex mismatch");
        gaussianH-=model.vertices[v].kappa*g;
    }
    require((h-gaussianH).norm()<1e-11,"Poisson physical H mismatch");
    // The four-occupation finite difference is an eight-fermion coefficient.
    // A Hamiltonian quadratic plus quadratic-operator squares cannot have it.
    const int bits[]={0,2,3,5};DataType coefficient=0.;
    for(unsigned mask=0;mask<16;++mask){unsigned state=0;for(int j=0;j<4;++j)if(mask&(1u<<j))state|=1u<<bits[j];coefficient+=(__builtin_popcount(mask)%2?-1.:1.)*h(state,state);}
    auto v=model.vertices[0];double expected=-4*v.kappa*(std::cosh(v.a)-1)*(std::cosh(v.b)-1);
    require(std::abs(coefficient-expected)<1e-12 && std::abs(coefficient)>1e-4,"Missing genuine eight-fermion interaction");
    auto line=poissonModel(3);auto a=line.generator(1),b=line.generator(5);
    require((a*b-b*a).norm()>1e-3,"Adjacent-edge vertices accidentally commute");
    // A common Wei contraction structure for both signs would force every
    // real eigenspace to be even-dimensional. The density direction has an
    // eigenvalue +.19 of multiplicity three, invariant under basis changes.
    auto plus=model.generator(4),minus=model.generator(5);
    require((plus+minus).norm()<1e-12,"Missing opposite density vertices");
    Eigen::SelfAdjointEigenSolver<MatType> spectrum(plus);int multiplicity=0;
    for(int i=0;i<plus.rows();++i)if(std::abs(spectrum.eigenvalues()(i)-.19)<1e-10)++multiplicity;
    require(multiplicity==3,"Wei exclusion witness lost odd real eigenspace");
    auto forward=model.polynomialVertex(1),opposite=model.polynomialVertex(3);
    require(forward(1,2).real()*opposite(1,2).real()<0.,"Missing opposite-sign hopping witness against a common occupation-phase gauge");
    // Complete CT expansion through order three, including every ordered word.
    const double beta=.53;MatType power=MatType::Identity(64,64);double factor=1.;
    for(int k=0;k<=3;++k) {
        if(k){power=(-h)*power;factor*=beta/k;}
        DataType sum=0.;int count=1;for(int i=0;i<k;++i)count*=model.vertices.size();
        for(int code=0;code<count;++code) {
            std::vector<int> word(k);int c=code;for(auto& j:word){j=c%model.vertices.size();c/=model.vertices.size();}
            auto z=model.signedWeight(word,beta),ref=model.fockWeight(word,beta);
            rec.add(2,"complete_CT_order",code,word,z,ref);sum+=z;
        }
        require(std::abs(sum-factor*power.trace())<1e-9*(1+std::abs(sum)),"CT ordering or factorial normalization wrong");
    }
    // Matrix series summed to order 18 checks the same CT scalar normalization
    // against the exact finite-system partition function beyond the word audit.
    MatType term=MatType::Identity(64,64),sum=term;
    for(int k=1;k<=18;++k){term=(-beta/k)*h*term;sum+=term;}
    require((sum-matrixExp(-beta*h)).norm()<1e-10,"CT series does not reproduce physical H");
}
void controls(Records& rec) {
    auto model=poissonModel(2);double b=2*std::asinh(std::pow(.5,.25));
    model.vertices={{0,1,'x','z',b,.17,.3},{0,1,'z','z',b,.23,.4},{0,1,'x','z',-b,-.17,.3},{0,1,'z','z',-b,-.23,.4}};
    std::vector<int> word={0,1,2,3};
    rec.add(2,"zero_fermion_weight",0,word,model.signedWeight(word,.7),model.fockWeight(word,.7),0);
    model=poissonModel(2);model.vertices[0].kappa=0.;
    rec.add(2,"zero_vertex_scalar",0,{0},model.signedWeight({0},.7),0.,0);
    model.vertices[0].kappa=-.2;bool rejected=false;try{model.validate();}catch(const std::exception&){rejected=true;}
    require(rejected,"Negative Poisson coefficient was accepted as sign-free");
    // Deliberately invalid negative coefficient: the *full* CT weight changes
    // sign even though the fermionic Gaussian trace alone stays positive.
    auto good=poissonModel(2);double negativeScalar=-.2*.7;
    auto pulse=good.generator(0);rec.add(2,"negative_coefficient_control",0,{0},negativeScalar*traceProduct({pulse}),negativeScalar*matrixExp(fockMatrix(pulse)).trace(),-1);
    // The Pfaffian must also see a genuine Gaussian negative Spin lift, without
    // concealing it behind a positive determinant square root.
    auto density=numberMajorana(MatType::Identity(1,1));std::vector<MatType> lift(4,DataType(0,std::acos(-1.)/2)*density);
    double scalar=std::pow(.6*.3,4)/24.;rec.add(0,"negative_spin_lift",0,{0,0,0,0},scalar*traceProduct(lift),-2.*scalar,-1);
}
}
int main(int argc,char** argv) {
    try {
        mkl_set_num_threads(1);mkl_set_dynamic(0);int samples=32;unsigned seed=731;std::string csv;
        for(int i=1;i<argc;++i){std::string a=argv[i];if(a=="--help"){std::cout<<"sign_free_poisson --samples N --seed N --csv PATH\nLocal exponential-vertex Hamiltonian and signed CT weight validator; no production update algorithm.\n";return 0;}require(i+1<argc,"Missing CLI value");std::string v=argv[++i];if(a=="--samples")samples=std::stoi(v);else if(a=="--seed")seed=std::stoul(v);else if(a=="--csv")csv=v;else throw std::runtime_error("Unknown argument "+a);}
        require(samples>0,"samples must be positive");std::ofstream file;if(!csv.empty()){file.open(csv);require(file.good(),"Cannot write CSV");}
        Records rec(csv.empty()?std::cout:file);physicalChecks(rec);controls(rec);std::mt19937_64 rng(seed);
        for(int cells:{2,3}) {
            auto model=poissonModel(cells);
            for(int sample=0;sample<samples;++sample) {
                // Each coefficient remains positive; amplitudes change sign and
                // magnitude to exercise the full real SL(2) generator family.
                std::uniform_real_distribution<double> amplitude(-.8,.8),coupling(.12,.4);
                for(auto& v:model.vertices){v.a=amplitude(rng);v.b=amplitude(rng);v.kappa=coupling(rng);}
                int k=1+int(rng()%12);std::vector<int> word(k);for(auto& v:word)v=rng()%model.vertices.size();
                auto z=model.signedWeight(word,.53),ref=model.fockWeight(word,.53);
                std::ostringstream parameters;parameters<<std::setprecision(17)<<"beta=.53";
                for(const auto& v:model.vertices)parameters<<';'<<v.r<<':'<<v.s<<':'<<v.axisR<<':'<<v.axisS<<':'<<v.a<<':'<<v.b<<':'<<v.kappa;
                rec.add(cells,"random_CT_word",sample,word,z,ref,1,parameters.str());
            }
        }
        std::cerr<<std::setprecision(12)<<"PASS poisson rows="<<rec.rows<<" positive="<<rec.positive<<" zero="<<rec.zeros<<" negative_controls="<<rec.negative<<" max_relative_error="<<rec.error<<" max_positive_phase_error="<<rec.phase<<" samples="<<samples<<" seed="<<seed<<"\n";
    }catch(const std::exception& e){std::cerr<<"FAIL poisson: "<<e.what()<<"\n";return 1;}return 0;
}
