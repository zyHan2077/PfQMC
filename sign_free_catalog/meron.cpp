#include "catalog/meron_models.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
using namespace catalog;
namespace {
struct Records {
    std::ostream& out;int rows=0,positive=0,negative=0,zeros=0;double error=0,phase=0;
    explicit Records(std::ostream& o):out(o){out<<"sites,kind,sample,beta,events,real,imag,reference_real,reference_imag,expected_sign,relative_error\n"<<std::setprecision(17);}
    void add(int sites,const std::string& kind,int sample,double beta,const std::vector<MeronEvent>& events,DataType z,DataType ref,int expected) {
        double e=std::abs(z-ref)/(expected==0?1.:std::max(1e-280,std::abs(ref)));error=std::max(error,e);
        std::ostringstream diagnostic;diagnostic<<std::setprecision(17)<<"Meron Pfaffian/Fock mismatch in "<<kind<<" sample="<<sample<<" z="<<z<<" ref="<<ref<<" relative="<<e<<" events=";
        for(const auto& event:events)diagnostic<<event.time<<':'<<int(event.kind)<<':'<<event.index<<';';
        require(e<3e-6,diagnostic.str());
        if(expected==0){require(std::abs(z)<1e-8,"Meron expected zero");++zeros;}
        else {require(std::abs(z)>0&&std::isfinite(std::abs(z)),"Meron invalid nonzero weight");double d=std::abs(z/std::abs(z)-double(expected));require(d<3e-6,"Meron wrong signed phase in "+kind);if(expected>0){++positive;phase=std::max(phase,d);}else ++negative;}
        out<<sites<<','<<kind<<','<<sample<<','<<beta<<',';for(int k=0;k<int(events.size());++k){if(k)out<<';';out<<events[k].time<<':'<<int(events[k].kind)<<':'<<events[k].index;}
        out<<','<<z.real()<<','<<z.imag()<<','<<ref.real()<<','<<ref.imag()<<','<<expected<<','<<e<<'\n';++rows;
    }
};
int expected(const MeronModel& m,const std::vector<MeronEvent>& events) {
    int monomers=0;for(const auto& e:events)if(e.kind==MeronKind::Monomer)++monomers;
    return (monomers%2 || m.clusters(events).merons)?0:1;
}
// Some symmetric small graphs have extra exact fermion zeros. Certify these
// with integer Fock matrices, not a floating-point threshold. If K^3=cK then
// every exp(tK) lies in span{I,K,K^2}. Vanishing of ALL word coefficients proves
// zero for arbitrary time intervals. This is only a bounded validation oracle.
bool integerZeroCertificate(const MeronModel& m,const std::vector<MeronEvent>& events) {
    if(m.sites>4 || events.size()>8)return false;
    using IntMatrix=Eigen::Matrix<long long,Eigen::Dynamic,Eigen::Dynamic>;
    int n=m.sites,dim=1<<n;IntMatrix k=IntMatrix::Zero(dim,dim);
    for(unsigned state=0;state<unsigned(dim);++state)for(auto [i,j]:m.edges)for(int direction=0;direction<2;++direction) {
        int a=direction?j:i,b=direction?i:j;
        if((state&(1u<<b))&&!(state&(1u<<a))){unsigned u=state^(1u<<b),v=u|(1u<<a);k(v,state)+=(__builtin_popcount(state&((1u<<b)-1))+__builtin_popcount(u&((1u<<a)-1)))%2?-1:1;}
    }
    IntMatrix k2=k*k,k3=k2*k;long long c=0;bool found=false;
    for(int i=0;i<dim&&!found;++i)for(int j=0;j<dim;++j)if(k(i,j)){if(k3(i,j)%k(i,j))return false;c=k3(i,j)/k(i,j);found=true;break;}
    if(!found || (k3-c*k).cwiseAbs().maxCoeff()!=0)return false;
    std::vector<std::vector<int>> diagonals;
    for(const auto& e:events) {
        if(e.kind==MeronKind::SingletBond)continue;std::vector<int> d(dim);
        for(int state=0;state<dim;++state){if(e.kind==MeronKind::Monomer)d[state]=2*((state>>e.index)&1)-1;else {auto [i,j]=m.edges[e.index];d[state]=(2*((state>>i)&1)-1)*(2*((state>>j)&1)-1);}}
        diagonals.push_back(d);
    }
    if(diagonals.empty())return false;
    std::vector<IntMatrix> powers={IntMatrix::Identity(dim,dim),k,k2};
    std::function<bool(int,const IntMatrix&)> allZero=[&](int depth,const IntMatrix& product) {
        if(depth==int(diagonals.size()))return product.trace()==0;
        for(const auto& power:powers){IntMatrix next=power*product;for(int i=0;i<dim;++i)next.row(i)*=diagonals[depth][i];if(!allZero(depth+1,next))return false;}return true;
    };
    return allZero(0,IntMatrix::Identity(dim,dim));
}
MatType annihilation(int n,int i) {
    MatType a=MatType::Zero(1<<n,1<<n);for(unsigned s=0;s<(1u<<n);++s)if(s&(1u<<i))a(s^(1u<<i),s)=__builtin_popcount(s&((1u<<i)-1))%2?-1.:1.;return a;
}
void physicalChecks() {
    MeronModel m;auto h=m.exactHamiltonian();require((h-h.adjoint()).norm()<1e-12,"Meron physical H not Hermitian");
    // Nonzero eight-fermion coefficient distinguishes this H from the old
    // quartic PairHopping/mobile-XXZ Hamiltonians, including after Gaussian bases.
    DataType coefficient=0.;for(unsigned bits=0;bits<16;++bits)coefficient+=(__builtin_popcount(bits)%2?-1.:1.)*h(bits<<2,bits<<2);
    require(std::abs(coefficient+m.exchange)<1e-11,"Missing -J E_i E_j/4 eight-fermion interaction");
    int n=6,dim=64;MatType id=MatType::Identity(dim,dim);std::vector<MatType> c,num,E,D,plus;
    for(int i=0;i<n;++i){c.push_back(annihilation(n,i));num.push_back(c.back().adjoint()*c.back());}
    MatType sy=MatType::Zero(dim,dim),marshall=sy;
    for(int i=0;i<2;++i){int a=2+2*i,b=a+1;E.push_back(id-num[a]-num[b]+2.*num[a]*num[b]);D.push_back(num[a]+num[b]-id);plus.push_back(c[a].adjoint()*c[b].adjoint());sy+=DataType(0,-.5)*(plus.back()-plus.back().adjoint());if(i==1)marshall+=.5*D.back();}
    double pi=std::acos(-1.);auto u=matrixExp(DataType(0,pi/2)*sy),v=matrixExp(DataType(0,pi)*marshall);
    MatType rotated=-m.hopping*(c[0].adjoint()*c[1]+c[1].adjoint()*c[0])+m.repulsion*(num[0]-.5*id)*(num[1]-.5*id);
    rotated+=.25*m.exchange*(D[0]*D[1]-E[0]*E[1])-.5*m.exchange*(plus[0]*plus[1].adjoint()+plus[0].adjoint()*plus[1]);
    for(int i=0;i<2;++i)rotated-=m.coupling*(i==0?1.:-1.)*(num[i]-.5*id)*num[2+2*i]*num[3+2*i];
    require((v*u*h*u.adjoint()*v.adjoint()-rotated).norm()<1e-10,"Meron rotation/Marshall Hamiltonian identity failed");
    // Rank-one pair singlet is not a Gaussian pure state: its one-body density
    // matrix is I/2, not an idempotent Slater projector at fixed particle number.
    MatType state=MatType::Zero(16,1);state(3,0)=1/std::sqrt(2.);state(12,0)=-1/std::sqrt(2.);MatType correlation=MatType::Zero(4,4);
    for(int i=0;i<4;++i)for(int j=0;j<4;++j)correlation(i,j)=(state.adjoint()*annihilation(4,i).adjoint()*annihilation(4,j)*state)(0,0);
    require((correlation-.5*MatType::Identity(4,4)).norm()<1e-12 && (correlation*correlation-correlation).norm()>.4,"Singlet nongGaussian rank witness failed");
}
void controls(Records& rec) {
    MeronModel m;m.repulsion=0.;m.exchange=m.coupling=1.;double beta=1.2;
    std::vector<MeronEvent> zero={{.1,MeronKind::Monomer,0},{.4,MeronKind::Monomer,1},{.9,MeronKind::SingletBond,0}};
    std::vector<std::uint64_t> negativePath={3,0,12,3};
    double opposite=.5*(std::cosh(.48)-1),same=.5*(std::cosh(.48)+1);
    rec.add(2,"negative_resolved_path",0,beta,zero,m.resolvedPfaffian(beta,zero,negativePath),-opposite/8.,-1);
    int negatives=0;auto sum=m.exhaustivePfaffian(beta,zero,&negatives);require(negatives>0,"Meron test never visited a negative resolved path");
    rec.add(2,"meron_signed_resummation_zero",0,beta,zero,sum,0.,0);
    rec.add(2,"meron_cluster_zero",0,beta,zero,m.groupedWeight(beta,zero),m.exactWordTrace(beta,zero),0);
    auto positive=zero;positive[1].index=0;
    rec.add(2,"positive_grouped_word",0,beta,positive,m.groupedWeight(beta,positive),same/2.,1);
    rec.add(2,"positive_full_Pfaffian_sum",0,beta,positive,m.exhaustivePfaffian(beta,positive),same/2.,1);
    // Splitting P=(1/2+Sx) and keeping only Sx changes the grouped CT word.
    // Calculate its spin factor independently; the result is a real negative
    // full word, not merely one negative resolved path of a positive group.
    MatType sx=MatType::Zero(4,4),sy=sx;for(int s=0;s<4;++s){sx(s^1,s)=.5;sy(s^2,s)=.5;}
    MatType singlet=MatType::Zero(4,4);singlet(1,1)=singlet(2,2)=.5;singlet(1,2)=singlet(2,1)=-.5;
    DataType spin=(singlet*sy*sx).trace();require(std::abs(spin+.25)<1e-12,"Split monomer control construction");
    rec.add(2,"split_projector_negative_word",0,beta,zero,spin*m.fermionPfaffian(beta,zero),-opposite/4.,-1);
    // Odd blockers are part of the full Fock space, not silently projected out.
    std::vector<MeronEvent> one={{.3,MeronKind::Monomer,0}};
    require(m.clusterTrace({})==16. && m.clusterTrace(one)==4.,"Full-Fock untouched/blocker multiplicity wrong");
    auto blocked=std::vector<std::uint64_t>{1,1};rec.add(2,"blocked_spin_path_zero",0,beta,one,m.resolvedPfaffian(beta,one,blocked),0.,0);
    auto bad=m;bad.edges={{0,1},{1,2},{2,0}};bad.sites=3;bool rejected=false;try{bad.validate();}catch(const std::exception&){rejected=true;}require(rejected,"Nonbipartite meron graph accepted");
    auto square=MeronModel::square();std::vector<MeronEvent> extraZero={{.2,MeronKind::DensityBond,1},{.3,MeronKind::DensityBond,0}};
    require(integerZeroCertificate(square,extraZero),"Uniform square exact-zero certificate failed");
    rec.add(4,"integer_certified_fermion_zero",0,1.1,extraZero,square.groupedWeight(1.1,extraZero),square.referenceWeight(1.1,extraZero),0);
}
void completeSmallWords(Records& rec) {
    MeronModel m;double beta=.73;
    for(int order=0;order<=2;++order) {
        int count=1<<(2*order);
        for(int code=0;code<count;++code) {
            std::vector<MeronEvent> events;int bits=code;
            for(int i=0;i<order;++i){int v=bits%4;bits/=4;events.push_back({beta*(i+1)/(order+1),v==0?MeronKind::DensityBond:v==1?MeronKind::SingletBond:MeronKind::Monomer,v<2?0:v-2});}
            int sign=expected(m,events);auto oracle=m.exactWordTrace(beta,events);
            if(sign==1&&std::abs(oracle)<1e-11){require(integerZeroCertificate(m,events),"Unresolved near-zero small word");sign=0;}
            rec.add(2,"all_small_words",code,beta,events,m.groupedWeight(beta,events),oracle,sign);
            rec.add(2,"all_small_resolved_sums",code,beta,events,m.exhaustivePfaffian(beta,events),oracle,sign);
        }
    }
}
}
int main(int argc,char** argv) {
    try {
        mkl_set_num_threads(1);mkl_set_dynamic(0);int samples=16;unsigned seed=731;std::string csv;
        for(int i=1;i<argc;++i){std::string a=argv[i];if(a=="--help"){std::cout<<"sign_free_meron --samples N --seed N --csv PATH\nFull-Fock local singlet/conditional-density Hamiltonian; grouped signed-Pfaffian CT validator, no production update algorithm.\n";return 0;}require(i+1<argc,"Missing CLI value");std::string v=argv[++i];if(a=="--samples")samples=std::stoi(v);else if(a=="--seed")seed=std::stoul(v);else if(a=="--csv")csv=v;else throw std::runtime_error("Unknown argument "+a);}
        require(samples>0,"samples must be positive");std::ofstream file;if(!csv.empty()){file.open(csv);require(file.good(),"Cannot open CSV");}
        Records rec(csv.empty()?std::cout:file);physicalChecks();controls(rec);completeSmallWords(rec);std::mt19937_64 rng(seed);
        for(auto model:{MeronModel{},MeronModel::square()})for(int sample=0;sample<samples;++sample) {
            double beta=1.1;int order=1+int(rng()%7);std::vector<double> times(order);std::uniform_real_distribution<double> uniform(.03,beta-.03);for(auto& t:times)t=uniform(rng);std::sort(times.begin(),times.end());
            std::vector<MeronEvent> events;for(int k=0;k<order;++k){int kind=rng()%3,index=kind==2?rng()%model.sites:rng()%model.edges.size();events.push_back({times[k],MeronKind(kind),index});}
            int sign=expected(model,events);auto z=model.groupedWeight(beta,events),ref=model.referenceWeight(beta,events);
            if(sign==1&&std::abs(ref)<1e-11){require(integerZeroCertificate(model,events),"Unresolved near-zero random word");sign=0;}
            rec.add(model.sites,"random_cluster_word",sample,beta,events,z,ref,sign);
            if(model.sites==2)rec.add(2,"random_full_Fock_word",sample,beta,events,z,model.exactWordTrace(beta,events),sign);
        }
        std::cerr<<std::setprecision(12)<<"PASS meron rows="<<rec.rows<<" positive="<<rec.positive<<" zero="<<rec.zeros<<" negative_controls="<<rec.negative<<" max_relative_error="<<rec.error<<" max_positive_phase_error="<<rec.phase<<" samples="<<samples<<" seed="<<seed<<"\n";
    }catch(const std::exception& e){std::cerr<<"FAIL meron: "<<e.what()<<"\n";return 1;}return 0;
}
