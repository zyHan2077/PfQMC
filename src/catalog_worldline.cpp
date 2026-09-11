#include "catalog/catalog.h"
namespace catalog {
std::vector<std::pair<int,int>> PairHopping::edges() const {
    require(lx>0 && ly>0 && 2*lx*ly<=10,"Worldline reference supports up to 10 fermion modes");
    std::vector<std::pair<int,int>> e;
    for(int x=0;x<lx;++x) for(int y=0;y<ly;++y) {
        int r=x*ly+y; if(x+1<lx) e.emplace_back(r,r+ly); if(y+1<ly) e.emplace_back(r,r+1);
    }
    return e;
}
double PairHopping::diagonalEnergy(unsigned v) const {
    double e=-mu*__builtin_popcount(v);
    for(auto bond:edges()) {
        int ni=__builtin_popcount((v>>(2*bond.first))&3u), nj=__builtin_popcount((v>>(2*bond.second))&3u);
        e+=repulsion*(ni-1)*(nj-1);
    }
    return e;
}
static MatType pairRotation(int n,int site,double theta) {
    MatType c=pairMajorana(n,2*site,2*site+1,1.,false);
    // Q(A^dagger)=Q(A)^dagger, including anti-Hermitian generators.
    return theta*(c-c.adjoint()).eval();
}
static Gaussian projector(int n,unsigned v) {
    Gaussian p; p.trace=1.; p.g=MatType::Zero(2*n,2*n);
    for(int i=0;i<n;++i) { p.g(2*i,2*i+1)=((v&(1u<<i)) ? -1.i:1.i); p.g(2*i+1,2*i)=-p.g(2*i,2*i+1); }
    return p; // b is deliberately unused: a pure projector has no invertible B.
}
static MatType moveRotation(const PairHopping& model,unsigned before,unsigned after) {
    int n=2*model.lx*model.ly;
    for(auto e:model.edges()) for(int reverse=0;reverse<2;++reverse) {
        int source=reverse ? e.second:e.first, dest=reverse ? e.first:e.second;
        unsigned s=3u<<(2*source),d=3u<<(2*dest);
        if((before&s)==s && (before&d)==0 && after==(before^s^d))
            return pairRotation(n,dest,.5*std::acos(-1.))+pairRotation(n,source,-.5*std::acos(-1.));
    }
    throw std::runtime_error("Illegal pair hop (single occupations block motion)");
}
MatType PairHopping::exactHamiltonian() const {
    auto e=edges(); int n=2*lx*ly,dim=1<<n;
    require(hopping>=0.,"Pair hopping t must be nonnegative");
    MatType h=MatType::Zero(dim,dim);
    for(unsigned v=0;v<unsigned(dim);++v) {
        h(v,v)=diagonalEnergy(v);
        for(auto bond:e) for(int reverse=0;reverse<2;++reverse) {
            int i=reverse ? bond.second:bond.first,j=reverse ? bond.first:bond.second;
            unsigned si=3u<<(2*i),sj=3u<<(2*j);
            if((v&si)==si && (v&sj)==0) h(v^si^sj,v)=-hopping;
        }
    }
    return h;
}
DataType PairHopping::pathPfaffian(const std::vector<unsigned>& path,double dt) const {
    (void)edges(); require(path.size()>=2 && path.front()==path.back(),"Worldline must close");
    require(dt>0. && hopping>=0.,"Invalid worldline step or hopping");
    int n=2*lx*ly; std::vector<Gaussian> factors; double scalar=1.;
    for(size_t k=0;k+1<path.size();++k) {
        require(path[k]<(1u<<n),"Worldline occupation outside Fock space");
        factors.push_back(projector(n,path[k]));
        scalar*=std::exp(-dt*diagonalEnergy(path[k]));
        if(path[k]!=path[k+1]) {
            factors.push_back(gaussian(moveRotation(*this,path[k],path[k+1])));
            scalar*=dt*hopping;
        }
    }
    return scalar*traceGaussianFactors(factors);
}
DataType PairHopping::closedHopPfaffian(unsigned before,unsigned after) const {
    int n=2*lx*ly;
    // A two-vertex closed loop; this is not an off-diagonal matrix element.
    return traceGaussianFactors({projector(n,before),gaussian(moveRotation(*this,before,after)),
                                 projector(n,after),gaussian(moveRotation(*this,after,before))});
}
}
