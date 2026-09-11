#include "catalog/projected_models.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
using namespace catalog;
namespace {
struct Records {
    std::ostream& out;int rows=0,positive=0,negative=0,zeros=0;double error=0.,phase=0.;
    explicit Records(std::ostream& o):out(o) {out<<"model,kind,sample,sector,real,imag,reference_real,reference_imag,expected_sign,relative_error\n"<<std::setprecision(17);}
    void add(const std::string& model,const std::string& kind,int sample,const std::string& sector,DataType z,DataType ref,int sign=1) {
        double e=std::abs(z-ref)/(1+std::abs(ref));error=std::max(error,e);
        require(e<3e-7,model+" projection/Fock mismatch");
        if(sign==0) {require(std::abs(z)<1e-7,model+" expected zero");++zeros;}
        else {require(std::abs(z)>1e-11,model+" unresolved zero");double d=std::abs(z/std::abs(z)-double(sign));require(d<3e-7,model+" wrong signed phase");if(sign==1){++positive;phase=std::max(phase,d);}else ++negative;}
        out<<model<<','<<kind<<','<<sample<<','<<sector<<','<<z.real()<<','<<z.imag()<<','<<ref.real()<<','<<ref.imag()<<','<<sign<<','<<e<<'\n';++rows;
    }
};
std::string label(const std::vector<int>& sector) {std::string out;for(int k:sector){if(!out.empty())out+=':';out+=std::to_string(k);}return out;}
// Rebuild each bilinear directly by fermion bit operations on a selected global
// occupation basis. This does not use fockMatrix or the production projection.
MatType physicalOracle(const ProjectedModel& model,const std::vector<int>& sector) {
    int n=model.base.modes,local=model.cellModes;std::vector<unsigned> states;
    for(unsigned s=0;s<(1u<<n);++s) {
        bool good=true;for(int r=0;r<model.cells;++r){int k=__builtin_popcount((s>>(r*local))&((1u<<local)-1));good&=(model.kind==ProjectedModel::Charge?k:k%2)==sector[r];}
        if(good)states.push_back(s);
    }
    // Match the tensor order: cell zero is the slow tensor index.
    std::sort(states.begin(),states.end(),[&](unsigned a,unsigned b){for(int r=0;r<model.cells;++r){unsigned x=(a>>(r*local))&((1u<<local)-1),y=(b>>(r*local))&((1u<<local)-1);if(x!=y)return x<y;}return false;});
    require(states.size()<=512,"Physical oracle too large");
    std::vector<int> index(1<<n,-1);for(int i=0;i<int(states.size());++i)index[states[i]]=i;
    auto bilinear=[&](const MatType& a) {
        auto h=numberMatrix(a);MatType q=MatType::Zero(states.size(),states.size());
        for(int col=0;col<int(states.size());++col) {
            unsigned s=states[col];q(col,col)-=.5*h.trace();
            for(int j=0;j<n;++j)if(s&(1u<<j)) {
                unsigned u=s^(1u<<j);int sj=__builtin_popcount(s&((1u<<j)-1));
                for(int i=0;i<n;++i)if(!(u&(1u<<i))) {
                    unsigned v=u|(1u<<i);int row=index[v];
                    if(row>=0)q(row,col)+=((sj+__builtin_popcount(u&((1u<<i)-1)))%2?-1.:1.)*h(i,j);
                }
            }
        }return q;
    };
    MatType h=-bilinear(model.base.q0);
    for(const auto& c:model.base.channels){auto o=bilinear(c.a);h-=.5*c.coupling*o*o;}
    for(int i=0;i<int(states.size());++i) {
        double energy=model.base.scalarEnergy;int previous=0;
        for(int r=0;r<model.cells;++r) {
            int k=__builtin_popcount((states[i]>>(r*local))&((1u<<local)-1));
            if(model.kind==ProjectedModel::Charge){energy+=model.onsite*(k-model.target)*(k-model.target)+model.bias*k;if(r)energy+=model.neighbor*k*previous;previous=k;}
            else {int p=k%2?-1:1;energy+=model.bias*p;if(r)energy+=model.neighbor*p*previous;previous=p;}
        }h(i,i)+=energy;
    }return h;
}
void verifyPhysical(ProjectedModel model) {
    std::vector<int> sector(model.cells,model.kind==ProjectedModel::Charge?1:0);
    for(int variant=0;variant<2;++variant) {
        sector.back()=model.kind==ProjectedModel::Charge?variant+1:variant;
        auto a=model.exactSectorHamiltonian(sector),b=physicalOracle(model,sector);
        require((a-b).norm()<1e-10*(1+b.norm()),"Physical charging/parity Hamiltonian mismatch");
        require((a-a.adjoint()).norm()<1e-11,"Physical projected Hamiltonian not Hermitian");
    }
    // Integrate the actual real-HS matrix factors with Gaussian quadrature in
    // a small sector, then check convergence to the independently rebuilt H.
    // The charging/parity scalar is exact; finite-step noncommuting factors
    // still carry the ordinary Trotter error.
    sector.assign(model.cells,model.kind==ProjectedModel::Charge?1:0);
    auto h=physicalOracle(model,sector);int dim=h.rows();double beta=.18;
    auto quadratic=[&](const MatType& a) {
        auto m=model;m.base.q0=a;m.base.channels.clear();m.base.scalarEnergy=0.;m.onsite=m.neighbor=m.bias=0.;
        return MatType(-m.exactSectorHamiltonian(sector));
    };
    auto q0=quadratic(model.base.q0);
    std::vector<MatType> directions;double commutator=0.;
    for(const auto& c:model.base.channels) {
        auto o=quadratic(c.a);commutator=std::max(commutator,(q0*o-o*q0).norm());directions.push_back(o);
    }
    require(commutator>1e-5,"Projection test accidentally used only commuting fields");
    MatType jac=MatType::Zero(12,12);for(int k=1;k<12;++k)jac(k-1,k)=jac(k,k-1)=std::sqrt(double(k));
    Eigen::SelfAdjointEigenSolver<MatType> gh(jac);
    auto exact=matrixExp(-beta*h);double previous=1e100;
    for(int slices:{2,4,8}) {
        double dt=beta/slices;auto kinetic=matrixExp(.5*dt*q0);MatType step=kinetic;
        for(int c=0;c<int(directions.size());++c) {
            MatType integral=MatType::Zero(dim,dim);
            for(int k=0;k<12;++k) integral+=std::norm(gh.eigenvectors()(0,k))*matrixExp(std::sqrt(dt*model.base.channels[c].coupling)*gh.eigenvalues()(k)*directions[c]);
            step=integral*step;
        }
        step=std::exp(-dt*(model.sectorEnergy(sector)+model.base.scalarEnergy))*kinetic*step;
        MatType product=MatType::Identity(dim,dim);for(int k=0;k<slices;++k)product=step*product;
        double error=(product-exact).norm()/(1+exact.norm());
        require(error<.8*previous,"Projected HS does not converge to physical Hamiltonian");previous=error;
    }
    require(previous<.002,"Projected HS physical-H error too large at eight slices");
}
void verify(Records& rec,ProjectedModel model,std::mt19937_64& rng,int samples) {
    verifyPhysical(model);std::normal_distribution<double> normal;
    const double beta=.61;const int slices=3;int range=model.kind==ProjectedModel::Charge?model.cellModes+1:2;
    std::string id=model.base.id+"_"+std::to_string(model.cells)+(model.neighbor<0?"_negative_neighbor":"_positive_neighbor");
    for(int s=0;s<samples;++s) {
        std::vector<double> fields(model.base.fieldCount(slices));for(auto& x:fields)x=normal(rng)*1.2;
        auto pulses=model.base.pulses(slices,beta,fields);
        std::vector<std::vector<DataType>> pf,oracle;
        for(int r=0;r<model.cells;++r){pf.push_back(model.cellWeights(pulses,r));oracle.push_back(model.cellWeights(pulses,r,true));}
        int combinations=1;for(int r=0;r<model.cells;++r)combinations*=range;
        for(int code=0;code<combinations;++code) {
            std::vector<int> sector(model.cells);int c=code;for(int r=0;r<model.cells;++r){sector[r]=c%range;c/=range;}
            double scalar=std::exp(-beta*(model.sectorEnergy(sector)+model.base.scalarEnergy));DataType z=scalar,ref=scalar;
            for(int r=0;r<model.cells;++r){z*=pf[r][sector[r]];ref*=oracle[r][sector[r]];}
            rec.add(id,"all_sectors",s,label(sector),z,ref);
        }
        // Positive charging is not a negative-square HS: it is applied once as
        // the exact conserved-sector scalar. Test an actual hard restriction.
        model.hard.assign(model.cells,1);std::vector<int> chosen(model.cells,1),excluded=chosen;excluded[0]=0;
        DataType z=model.weight(pulses,chosen,beta),ref=model.weight(pulses,chosen,beta,true);
        rec.add(id,"hard_sector",s,label(chosen),z,ref);
        rec.add(id,"hard_excluded",s,label(excluded),model.weight(pulses,excluded,beta),0.,0);
        model.hard.clear();
    }
}
void controls(Records& rec) {
    auto x=numberMajorana(spinX(2)),z=numberMajorana(spinZ(2));double a=2*std::log(2.);
    std::vector<MatType> p={a*x,a*z,-a*x,-a*z};
    auto pf=canonicalPfaffians(p),fock=canonicalFock(p);
    rec.add("spin1_fixed_number_failure","negative",0,"1",pf[1],-1887./4096.,-1);
    for(int k=0;k<4;++k)rec.add("spin1_fixed_number_oracle","control",0,std::to_string(k),pf[k],fock[k],k==1||k==2?-1:1);
    rec.add("spin1_parity_survives","even",0,"0",pf[0]+pf[2],2209./4096.);
    rec.add("spin1_parity_survives","odd",0,"1",pf[1]+pf[3],2209./4096.);
    double b=2*std::asinh(std::pow(.5,.25));p={b*x,b*z,-b*x,-b*z};
    pf=canonicalPfaffians(p);fock=canonicalFock(p);
    for(int parity=0;parity<2;++parity)rec.add("spin1_parity_zero","zero",0,std::to_string(parity),pf[parity]+pf[parity+2],fock[parity]+fock[parity+2],0);
    // An extra real chemical potential destroys one projected parity weight.
    p.push_back(std::log(2.)*numberMajorana(MatType::Identity(3,3)));
    pf=canonicalPfaffians(p);fock=canonicalFock(p);
    for(int parity=0;parity<2;++parity)rec.add("spin1_mu_breaks_parity","control",0,std::to_string(parity),pf[parity]+pf[parity+2],fock[parity]+fock[parity+2],parity?1:-1);
    // The actual Spin lift must survive the additional Fourier phases.
    auto density=numberMajorana(MatType::Identity(1,1));p.assign(4,DataType(0,std::acos(-1.)/2)*density);
    pf=canonicalPfaffians(p);for(int k=0;k<2;++k)rec.add("projected_spin_lift","negative",0,std::to_string(k),pf[k],-1.,-1);
}
}
int main(int argc,char** argv) {
    try {
        mkl_set_num_threads(1);mkl_set_dynamic(0);int samples=8;unsigned seed=731;std::string csv;
        for(int i=1;i<argc;++i){std::string a=argv[i];if(a=="--help"){std::cout<<"sign_free_projected --samples N --seed N --csv PATH\nSigned canonical/parity weight validator; not a production Monte Carlo sampler.\n";return 0;}require(i+1<argc,"Missing CLI value");std::string v=argv[++i];if(a=="--samples")samples=std::stoi(v);else if(a=="--seed")seed=std::stoul(v);else if(a=="--csv")csv=v;else throw std::runtime_error("Unknown argument "+a);}
        require(samples>0,"samples must be positive");std::ofstream file;if(!csv.empty()){file.open(csv);require(file.good(),"Cannot open CSV");}
        Records rec(csv.empty()?std::cout:file);std::mt19937_64 rng(seed);controls(rec);
        for(int cells:{1,2}){verify(rec,chargingModel(cells),rng,samples);verify(rec,parityModel(cells),rng,samples);}
        // Both signs of local intercell charging/parity interactions are legal.
        auto charge=chargingModel(2);charge.neighbor=-.24;verify(rec,charge,rng,1);
        auto parity=parityModel(2);parity.neighbor=.27;verify(rec,parity,rng,1);
        std::cerr<<std::setprecision(12)<<"PASS projected rows="<<rec.rows<<" positive="<<rec.positive<<" zero="<<rec.zeros<<" negative_controls="<<rec.negative<<" max_relative_error="<<rec.error<<" max_positive_phase_error="<<rec.phase<<" samples="<<samples<<" seed="<<seed<<"\n";
    }catch(const std::exception& e){std::cerr<<"FAIL projected: "<<e.what()<<"\n";return 1;}return 0;
}
