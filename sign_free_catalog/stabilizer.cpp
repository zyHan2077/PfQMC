#include "catalog/stabilizer_models.h"
#include <mkl.h>
#include <Eigen/Eigenvalues>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
using namespace catalog;
int main(int argc,char** argv) { try {
    int samples=32;unsigned seed=731;std::string csv;
    for(int i=1;i<argc;++i) {std::string a=argv[i];
        if(a=="--samples"&&i+1<argc) samples=std::stoi(argv[++i]);
        else if(a=="--seed"&&i+1<argc) seed=std::stoul(argv[++i]);
        else if(a=="--csv"&&i+1<argc) csv=argv[++i];else throw std::invalid_argument("Unknown CLI argument");}
    if(samples<0) throw std::invalid_argument("Negative sample count");mkl_set_num_threads(1);
    std::ofstream file;if(!csv.empty()){file.open(csv);if(!file)throw std::runtime_error("Cannot open CSV");
        file<<"case,real_pf,imag_pf,real_reference,imag_reference,relative_error,expected\n"<<std::setprecision(17);}
    int rows=0,pos=0,zero=0,neg=0,complex=0;double maxerr=0.,maxphase=0.;
    auto record=[&](std::string name,DataType p,DataType ref,int expected) {
        double err=std::abs(p-ref)/std::max(1.,std::abs(ref));maxerr=std::max(maxerr,err);
        if(expected!=2)maxphase=std::max(maxphase,std::abs(p.imag())/std::max(1.,std::abs(p)));
        if(err>2e-9||!std::isfinite(err))throw std::runtime_error(name+" signed Pf mismatch");
        if((expected==1&&ref.real()<=0)||(expected<0&&ref.real()>=0)||(expected==0&&ref!=DataType(0.))||(expected==2&&ref.imag()==0.))
            throw std::runtime_error(name+" invalid exact sign classification");
        ++rows;if(expected==2)++complex;else if(expected>0)++pos;else if(expected<0)++neg;else ++zero;
        if(file)file<<name<<','<<p.real()<<','<<p.imag()<<','<<ref.real()<<','<<ref.imag()<<','<<err<<','<<expected<<'\n';
    };
    auto test=[&](const StabilizerModel& m,const std::vector<int>& labels,std::string name,bool occupation=false) {
        auto ref=m.binaryOracle(labels,.7);auto pf=m.pfaffian(labels,.7);
        record(name,pf,ref,ref.real()>0?1:ref.real()<0?-1:0);
        if(occupation)record(name+"_occupation",pf,m.occupationOracle(labels,.7),ref.real()>0?1:ref.real()<0?-1:0);
    };
    auto m=toricStabilizers();
    for(int subset=0;subset<256;++subset) {
        std::vector<int> labels;for(int j=0;j<8;++j)if(subset>>j&1)labels.push_back(j);
        test(m,labels,"toric_subset_"+std::to_string(subset),subset==0||subset==15||subset==240||subset==17);
    }
    std::mt19937 rng(seed);
    for(int k=0;k<samples;++k) {
        std::vector<int> labels;for(int j=0,n=1+rng()%12;j<n;++j)labels.push_back(rng()%8);
        test(m,labels,"toric_random_"+std::to_string(k),k<4);
        if(k%2==0) {auto closed=labels;closed.insert(closed.end(),labels.begin(),labels.end());test(m,closed,"toric_closed_"+std::to_string(k));}
    }
    auto big=toricStabilizers(3,2);
    for(int k=0;k<samples;++k) {
        std::vector<int> labels;for(int j=0,n=1+rng()%8;j<n;++j)labels.push_back(rng()%12);
        test(big,labels,"toric_3x2_"+std::to_string(k));
    }
    auto bad=m;bad.vertices[0].sign=-1;
    if(bad.positiveGroup())throw std::runtime_error("Missed -I stabilizer relation");
    test(bad,{0,1,2,3},"negative_relation",true);
    auto disabled=m;disabled.vertices[0].coupling=0.;test(disabled,{0,0},"zero_coupling",true);
    for(auto entry:std::vector<std::pair<std::string,std::vector<StabilizerPauli>>>{
          {"global_negative_XZXZ",{{0,0},{0,2},{0,0},{0,2}}},
          {"global_positive_ZZ",{{0,2},{0,2}}},
          {"global_zero_X",{{0,0}}},
          {"global_phase_XYZ",{{0,0},{0,1},{0,2}}}}) {
        auto ref=stabilizerOccupationTrace(8,entry.second);
        // XYZ is a deliberately complex noncommuting word: retain its phase.
        if(entry.first=="global_phase_XYZ") {
            auto p=stabilizerPauliPfaffian(8,entry.second,true);
            if(std::abs(p-ref)>1e-8*std::abs(ref)||std::abs(ref.imag())!=65536.)throw std::runtime_error("Lost complex Pauli scalar");
            record(entry.first,p,ref,2);
        } else record(entry.first,stabilizerPauliPfaffian(8,entry.second,true),ref,ref.real()>0?1:ref.real()<0?-1:0);
    }
    StabilizerModel noncomm;noncomm.qubits=1;noncomm.vertices={{{{0,0}},1.,1},{{{0,2}},1.,1}};
    if(noncomm.positiveGroup())throw std::runtime_error("Missed anticommuting generators");
    test(noncomm,{0,1,0,1},"noncommuting_negative",true);
    // Exact complete physical Z: logical spectrum with the full 2^qubits
    // spectator degeneracy. No gauge/parity sector is discarded.
    Eigen::SelfAdjointEigenSolver<MatType> es(m.logicalHamiltonian());
    double z=0.;for(auto e:es.eigenvalues())z+=std::exp(-.7*e);z*=std::ldexp(1.,m.qubits);
    double formula=std::ldexp(1.,m.modes()),stars=1.,plaquettes=1.;
    for(int j=0;j<8;++j){double x=.7*m.vertices[j].coupling;formula*=std::cosh(x);(j<4?stars:plaquettes)*=std::tanh(x);}
    formula*=(1+stars)*(1+plaquettes);
    if(std::abs(z-formula)>1e-10*z)throw std::runtime_error("Physical partition function mismatch");
    // Local complete four-state encoding, including Y and spectator count.
    for(int a=0;a<3;++a)for(int b=0;b<3;++b) {
        auto ref=stabilizerOccupationTrace(1,{{0,a},{0,b}});
        record("local_pauli_"+std::to_string(a)+std::to_string(b),stabilizerPauliPfaffian(1,{{0,a},{0,b}}),ref,a==b?1:0);
    }
    // Bell XX,ZZ,-YY generators have a +I relation; their complete four-mode
    // physical Hamiltonian must contain four identical copies of logical H.
    StabilizerModel bell;bell.qubits=2;bell.vertices={{{{0,0},{1,0}},.3,1},{{{0,2},{1,2}},.2,1},{{{0,1},{1,1}},.1,-1}};
    if(!bell.positiveGroup())throw std::runtime_error("Missed signed +I relation");
    Eigen::SelfAdjointEigenSolver<MatType> physical(bell.physicalHamiltonian()),logical(bell.logicalHamiltonian());
    for(int i=0;i<16;++i)if(std::abs(physical.eigenvalues()(i)-logical.eigenvalues()(i/4))>1e-12)
        throw std::runtime_error("Complete physical spectrum multiplicity mismatch");
    for(int subset=0;subset<8;++subset){std::vector<int> w;for(int j=0;j<3;++j)if(subset>>j&1)w.push_back(j);test(bell,w,"bell_signed_relation_"+std::to_string(subset),true);}
    std::cout<<std::setprecision(17)<<"seed="<<seed<<" rows="<<rows<<" positive="<<pos<<" zero="<<zero<<" negative="<<neg
        <<" complex_controls="<<complex<<" max_relative_error="<<maxerr<<" max_phase="<<maxphase<<" full_partition="<<z<<" PASS\n";
    return 0;
} catch(const std::exception& e) {std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;} }
