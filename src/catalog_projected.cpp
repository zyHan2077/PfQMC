#include "catalog/projected_models.h"
#include <cmath>
namespace catalog {
namespace {
constexpr double pi=3.14159265358979323846;
std::vector<MatType> localPulses(const std::vector<MatType>& pulses,int cell,int n,int cells) {
    std::vector<MatType> out;
    for(const auto& a:pulses) {
        require(a.rows()==2*n*cells,"Projected pulse dimension");
        MatType block=MatType::Zero(a.rows(),a.cols());
        for(int r=0;r<cells;++r) block.block(2*n*r,2*n*r,2*n,2*n)=a.block(2*n*r,2*n*r,2*n,2*n);
        require((a-block).norm()<1e-11*(1+a.norm()),"Projection requires conserved cell particle numbers");
        out.push_back(a.block(2*n*cell,2*n*cell,2*n,2*n));
    }
    return out;
}
MatType kron(const MatType& a,const MatType& b) {
    MatType out(a.rows()*b.rows(),a.cols()*b.cols());
    for(int i=0;i<a.rows();++i) for(int j=0;j<a.cols();++j)
        out.block(i*b.rows(),j*b.cols(),b.rows(),b.cols())=a(i,j)*b;
    return out;
}
}
std::vector<DataType> canonicalPfaffians(const std::vector<MatType>& pulses) {
    require(!pulses.empty(),"Canonical projection needs a pulse");
    int n=pulses.front().rows()/2,m=n+1;if(m%2==0)++m;
    for(const auto& a:pulses) (void)numberMatrix(a);
    auto density=numberMajorana(MatType::Identity(n,n));
    std::vector<DataType> result(n+1,0.);
    for(int j=0;j<m;++j) {
        double theta=2*pi*j/m;
        auto p=pulses;
        // Splitting avoids an individually traceless exp(i*pi*(N-n/2)).
        p.push_back(DataType(0,theta/2)*density);
        p.push_back(DataType(0,theta/2)*density);
        DataType z=std::exp(DataType(0,theta*n/2))*traceProduct(p);
        for(int k=0;k<=n;++k) result[k]+=std::exp(DataType(0,-theta*k))*z/double(m);
    }
    return result;
}
std::vector<DataType> canonicalFock(const std::vector<MatType>& pulses) {
    require(!pulses.empty(),"Canonical Fock needs pulses");
    int n=pulses.front().rows()/2;require(n<=8,"Independent Fock oracle is limited to eight modes");
    MatType product=MatType::Identity(1<<n,1<<n);
    for(const auto& a:pulses) product=matrixExp(fockMatrix(a))*product;
    std::vector<DataType> result(n+1,0.);
    for(unsigned s=0;s<(1u<<n);++s)result[__builtin_popcount(s)]+=product(s,s);
    return result;
}
bool ProjectedModel::allowed(const std::vector<int>& sector) const {
    if(int(sector.size())!=cells)return false;
    for(int r=0;r<cells;++r) {
        if(sector[r]<0 || sector[r]>(kind==Charge?cellModes:1))return false;
        if(!hard.empty() && (int(hard.size())!=cells || sector[r]!=hard[r]))return false;
    }
    return true;
}
double ProjectedModel::sectorEnergy(const std::vector<int>& sector) const {
    require(int(sector.size())==cells,"Sector dimension");double e=0.;
    for(int r=0;r<cells;++r) {
        if(kind==Charge) e+=onsite*std::pow(sector[r]-target,2)+bias*sector[r];
        else e+=bias*(sector[r]?-1.:1.);
        if(r) e+=neighbor*(kind==Charge?double(sector[r]*sector[r-1]):((sector[r]+sector[r-1])%2?-1.:1.));
    }
    return e;
}
std::vector<DataType> ProjectedModel::cellWeights(const std::vector<MatType>& pulses,int cell,bool fock) const {
    require(cell>=0 && cell<cells,"Invalid cell");
    auto local=localPulses(pulses,cell,cellModes,cells);
    auto z=fock?canonicalFock(local):canonicalPfaffians(local);
    if(kind==Charge)return z;
    std::vector<DataType> p(2,0.);for(int k=0;k<=cellModes;++k)p[k%2]+=z[k];return p;
}
DataType ProjectedModel::weight(const std::vector<MatType>& pulses,const std::vector<int>& sector,double beta,bool fock) const {
    if(!allowed(sector))return 0.;
    DataType z=std::exp(-beta*(sectorEnergy(sector)+base.scalarEnergy));
    for(int r=0;r<cells;++r)z*=cellWeights(pulses,r,fock)[sector[r]];
    return z;
}
MatType ProjectedModel::exactSectorHamiltonian(const std::vector<int>& sector) const {
    require(allowed(sector),"Requested physical sector is excluded");
    std::vector<std::vector<int>> states(cells);
    int dim=1;
    for(int r=0;r<cells;++r) {
        for(int s=0;s<(1<<cellModes);++s) {
            int k=__builtin_popcount(unsigned(s));
            if((kind==Charge?k:k%2)==sector[r])states[r].push_back(s);
        }
        dim*=states[r].size();
    }
    require(dim>0 && dim<=512,"Exact tensor-sector oracle dimension exceeds 512");
    auto projected=[&](const MatType& a) {
        auto local=localPulses({a},0,cellModes,cells);(void)local;
        MatType total=MatType::Zero(dim,dim);
        for(int r=0;r<cells;++r) {
            auto f=fockMatrix(a.block(2*cellModes*r,2*cellModes*r,2*cellModes,2*cellModes));
            MatType b(states[r].size(),states[r].size());
            for(int i=0;i<b.rows();++i)for(int j=0;j<b.cols();++j)b(i,j)=f(states[r][i],states[r][j]);
            MatType term=MatType::Ones(1,1);
            for(int c=0;c<cells;++c)term=kron(term,c==r?b:MatType(MatType::Identity(states[c].size(),states[c].size())));
            total+=term;
        }
        return total;
    };
    MatType h=-projected(base.q0);
    for(const auto& c:base.channels) {auto o=projected(c.a);h-=.5*c.coupling*o*o;}
    h+=(sectorEnergy(sector)+base.scalarEnergy)*MatType::Identity(dim,dim);
    return h;
}
ProjectedModel chargingModel(int cells) {
    require(cells==1 || cells==2,"Charging reference supports one or two cells");
    ProjectedModel p;p.kind=ProjectedModel::Charge;p.cells=cells;p.cellModes=6;
    p.base=symmetricPowerModel(3,2,cells,1).h;p.base.id="projected_charging";return p;
}
ProjectedModel parityModel(int cells) {
    require(cells==1 || cells==2,"Parity reference supports one or two cells");
    ProjectedModel p;p.kind=ProjectedModel::Parity;p.cells=cells;p.cellModes=3;p.onsite=0.;p.neighbor=-.31;p.bias=.13;
    auto cell=spinCell(2,false);
    p.base.id="projected_parity";p.base.description="Local parity interactions and optional hard parity constraints";
    p.base.modes=3*cells;p.base.q0=MatType::Zero(6*cells,6*cells);
    for(int r=0;r<cells;++r) {
        p.base.q0.block(6*r,6*r,6,6)=cell.q0;
        for(const auto& c:cell.channels) {MatType a=MatType::Zero(6*cells,6*cells);a.block(6*r,6*r,6,6)=c.a;p.base.channels.push_back({c.name,a,c.coupling});}
    }
    if(cells==2) {
        MatType a=MatType::Zero(12,12);a.block(0,0,6,6)=numberMajorana(spinX(2));a.block(6,6,6,6)=numberMajorana(spinZ(2));
        p.base.channels.push_back({"shared_xz",a,.17});
    }
    p.base.validate();return p;
}
}
