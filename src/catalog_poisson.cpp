#include "catalog/poisson_models.h"
#include <cmath>
namespace catalog {
namespace {
MatType oneBody(char axis) {
    require(axis=='x' || axis=='z' || axis=='d',"Poisson vertex axis must be x, z, or d");
    if(axis=='d')return MatType::Identity(3,3);
    return axis=='x'?spinX(2):spinZ(2);
}
// Independent occupation-basis bilinear construction for the polynomial oracle.
MatType localSpin(char axis) {
    auto h=oneBody(axis);MatType q=-.5*h.trace()*MatType::Identity(8,8);
    for(unsigned state=0;state<8;++state)for(int j=0;j<3;++j)if(state&(1u<<j)) {
        unsigned after=state^(1u<<j);int sj=__builtin_popcount(state&((1u<<j)-1));
        for(int i=0;i<3;++i)if(!(after&(1u<<i))) {
            unsigned target=after|(1u<<i);
            q(target,state)+=((sj+__builtin_popcount(after&((1u<<i)-1)))%2?-1.:1.)*h(i,j);
        }
    }return q;
}
MatType polynomial(char axis,double a) {
    if(axis=='d') {
        MatType out=MatType::Zero(8,8);
        for(unsigned state=0;state<8;++state) {
            double value=1.;for(int i=0;i<3;++i)value*=std::cosh(a/2)+((state&(1u<<i))?1.:-1.)*std::sinh(a/2);
            out(state,state)=value;
        }return out;
    }
    auto s=localSpin(axis);
    require((s*s*s-s).norm()<1e-12,"Spin-one Fock cubic identity failed");
    return MatType::Identity(8,8)+std::sinh(a)*s+(std::cosh(a)-1)*s*s;
}
MatType kron(const MatType& a,const MatType& b) {
    MatType out(a.rows()*b.rows(),a.cols()*b.cols());
    for(int i=0;i<a.rows();++i)for(int j=0;j<a.cols();++j)
        out.block(i*b.rows(),j*b.cols(),b.rows(),b.cols())=a(i,j)*b;
    return out;
}
}
void PoissonModel::validate() const {
    require(cells>=2,"Poisson reference needs at least two cells");
    require(!vertices.empty(),"Poisson model has no vertices");
    for(const auto& v:vertices) {
        require(v.r>=0 && v.r<cells && v.s>=0 && v.s<cells && v.r!=v.s,"Invalid Poisson vertex cells");
        require(std::abs(v.r-v.s)==1,"Poisson reference vertices must be nearest neighbors");
        (void)oneBody(v.axisR);(void)oneBody(v.axisS);
        require(std::isfinite(v.a) && std::isfinite(v.b) && std::isfinite(v.kappa) && v.kappa>=0,"Poisson vertex requires finite real amplitudes and nonnegative kappa");
    }
}
MatType PoissonModel::generator(int vertex) const {
    require(vertex>=0 && vertex<int(vertices.size()),"Invalid Poisson vertex index");
    const auto& v=vertices[vertex];MatType h=MatType::Zero(3*cells,3*cells);
    h.block(3*v.r,3*v.r,3,3)=v.a*oneBody(v.axisR);
    h.block(3*v.s,3*v.s,3,3)=v.b*oneBody(v.axisS);
    return numberMajorana(h);
}
MatType PoissonModel::polynomialVertex(int vertex) const {
    require(cells<=2,"Full polynomial Fock oracle is limited to two cells");
    require(vertex>=0 && vertex<int(vertices.size()),"Invalid Poisson vertex index");
    const auto& v=vertices[vertex];MatType out=MatType::Ones(1,1);
    // Cell zero occupies the least significant three occupation bits.
    for(int r=cells-1;r>=0;--r)out=kron(out,r==v.r?polynomial(v.axisR,v.a):r==v.s?polynomial(v.axisS,v.b):MatType(MatType::Identity(8,8)));
    return out;
}
MatType PoissonModel::exactHamiltonian() const {
    validate();require(cells<=2,"Exact physical H limited to six modes");
    MatType h=MatType::Zero(1<<(3*cells),1<<(3*cells));
    for(int v=0;v<int(vertices.size());++v)h-=vertices[v].kappa*polynomialVertex(v);
    return h;
}
double PoissonModel::orderedTimeScalar(const std::vector<int>& word,double beta) const {
    require(beta>=0 && std::isfinite(beta),"Invalid beta");double scalar=1.;
    for(int i=0;i<int(word.size());++i) {
        require(word[i]>=0 && word[i]<int(vertices.size()),"Invalid vertex word");
        scalar*=beta*vertices[word[i]].kappa/(i+1);
    }
    return scalar; // Integrates 0<tau_1<...<tau_k<beta, including k=0.
}
DataType PoissonModel::signedWeight(const std::vector<int>& word,double beta) const {
    validate();std::vector<MatType> pulses;
    for(int v:word)pulses.push_back(generator(v));
    if(word.empty())pulses.push_back(MatType::Zero(6*cells,6*cells));
    return orderedTimeScalar(word,beta)*traceProduct(pulses);
}
DataType PoissonModel::fockWeight(const std::vector<int>& word,double beta) const {
    validate();std::vector<MatType> products(cells,MatType::Identity(8,8));
    for(int index:word) {
        require(index>=0 && index<int(vertices.size()),"Invalid vertex word");const auto& v=vertices[index];
        products[v.r]=polynomial(v.axisR,v.a)*products[v.r];
        products[v.s]=polynomial(v.axisS,v.b)*products[v.s];
    }
    DataType z=orderedTimeScalar(word,beta);for(const auto& p:products)z*=p.trace();return z;
}
PoissonModel poissonModel(int cells) {
    PoissonModel m;m.cells=cells;
    for(int r=0;r<cells-1;++r) {
        m.vertices.push_back({r,r+1,'z','z',.31,-.27,.23});
        m.vertices.push_back({r,r+1,'x','z',.38,.21,.19});
        m.vertices.push_back({r,r+1,'z','x',-.24,.33,.17});
        // Opposite x hopping on the same cell prevents one occupation-phase
        // choice from making every Gaussian vertex entrywise nonnegative.
        m.vertices.push_back({r,r+1,'x','z',-.29,-.25,.13});
    }
    for(int r=0;r<cells;++r) {
        int s=r==cells-1?r-1:r+1;
        m.vertices.push_back({r,s,'d','z',.19,0.,.11});
        m.vertices.push_back({r,s,'d','z',-.19,0.,.09});
    }
    m.validate();return m;
}
}
