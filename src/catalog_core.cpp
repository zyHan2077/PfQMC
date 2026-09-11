#include "catalog/catalog.h"
#include <unsupported/Eigen/MatrixFunctions>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace catalog {
void require(bool ok,const std::string& s) { if(!ok) throw std::runtime_error(s); }
MatType matrixExp(const MatType& a) {
    require(a.allFinite(),"Nonfinite exponential input");
    if ((a-a.adjoint()).norm()<1e-12*std::max(1.,a.norm())) {
        Eigen::SelfAdjointEigenSolver<MatType> es(a);
        require(es.info()==Eigen::Success,"Hermitian diagonalization failed");
        return es.eigenvectors()*es.eigenvalues().array().exp().matrix().asDiagonal()*es.eigenvectors().adjoint();
    }
    return a.exp();
}
MatType numberMajorana(const MatType& h) {
    const int n=h.rows(); MatType a=MatType::Zero(2*n,2*n);
    for(int i=0;i<n;++i) for(int j=0;j<n;++j) {
        a(2*i,2*j)=a(2*i+1,2*j+1)=.5*(h(i,j)-h(j,i));
        a(2*i,2*j+1)=.5i*(h(i,j)+h(j,i));
        a(2*i+1,2*j)=-.5i*(h(i,j)+h(j,i));
    }
    return a;
}
MatType numberMatrix(const MatType& a) {
    int n=a.rows()/2; MatType h(n,n);
    for(int i=0;i<n;++i) for(int j=0;j<n;++j) h(i,j)=a(2*i,2*j)-1.i*a(2*i,2*j+1);
    require((numberMajorana(h)-a).norm()<1e-9*std::max(1.,a.norm()),"Generator contains pairing");
    return h;
}
MatType pairMajorana(int n,int i,int j,DataType delta,bool annihilate) {
    require(i!=j && i>=0 && j>=0 && i<n && j<n,"Invalid pairing indices");
    MatType a=MatType::Zero(2*n,2*n);
    auto term=[&](int p,bool pd,int q,bool qd,DataType coefficient) {
        DataType u[2]={.5,pd ? -.5i:.5i}, v[2]={.5,qd ? -.5i:.5i};
        for(int k=0;k<2;++k) for(int l=0;l<2;++l) {
            a(2*p+k,2*q+l)+=2.*coefficient*u[k]*v[l];
            a(2*q+l,2*p+k)-=2.*coefficient*u[k]*v[l];
        }
    };
    term(i,true,j,true,delta);
    if(annihilate) term(j,false,i,false,std::conj(delta));
    return a;
}
static std::pair<unsigned,DataType> majoranaAction(int a,unsigned v) {
    unsigned mask=1u<<(a/2);
    double sign=__builtin_popcount(v&(mask-1))%2 ? -1.:1.;
    return {v^mask,sign*(a%2 ? ((v&mask) ? -1.i:1.i):DataType(1.))};
}
MatType fockMatrix(const MatType& a) {
    int n=a.rows()/2; require(n<=10,"Exact Fock oracle is limited to 10 modes");
    int dim=1<<n; MatType q=MatType::Zero(dim,dim);
    for(int i=0;i<2*n;++i) for(int j=i+1;j<2*n;++j) if(std::abs(a(i,j))>1e-16) {
        for(unsigned v=0;v<unsigned(dim);++v) {
            auto first=majoranaAction(j,v), second=majoranaAction(i,first.first);
            q(second.first,v)+=.5*a(i,j)*first.second*second.second;
        }
    }
    return q;
}
MatType occupationProjector(int modes,unsigned occupation) {
    require(modes<=10,"Fock projector size limit");
    MatType p=MatType::Zero(1<<modes,1<<modes); p(occupation,occupation)=1.; return p;
}
Gaussian gaussian(const MatType& a) {
    require(a.rows()==a.cols() && a.rows()%2==0 && (a+a.transpose()).norm()<1e-9*std::max(1.,a.norm()),"Invalid Majorana matrix");
    int d=a.rows(), n=d/2; MatType id=MatType::Identity(d,d);
    Gaussian out; out.b=matrixExp(a);
    MatType sinh=std::sqrt(.5)*(matrixExp(-a/4.)-matrixExp(a/4.));
    MatType eta=MatType::Zero(2*d,2*d);
    eta.topLeftCorner(d,d)=sinh; eta.bottomRightCorner(d,d)=sinh;
    eta.topRightCorner(d,d)=-id; eta.bottomLeftCorner(d,d)=id;
    out.trace=std::ldexp(1.,n)*(n%2 ? -1.:1.)*pfaf(d,eta);
    Eigen::FullPivLU<MatType> lu(id+out.b);
    require(lu.isInvertible(),"Individual Gaussian has zero trace: split this pulse first");
    out.g=2.*lu.inverse()-id;
    out.g=.5*(out.g-out.g.transpose()).eval();
    require(out.b.allFinite() && out.g.allFinite() && std::isfinite(std::abs(out.trace)),"Nonfinite Gaussian");
    return out;
}
DataType traceGaussianFactors(const std::vector<Gaussian>& supplied) {
    require(!supplied.empty(),"Empty product needs an explicit mode count");
    auto factors=supplied; int d=factors[0].g.rows(), n=d/2;
    MatType id=MatType::Identity(d,d);
    if(factors.size()%2) factors.push_back({id,MatType::Zero(d,d),std::ldexp(1.,n)});
    const int len=factors.size(); require(d*len<=8192,"Spacetime Pfaffian exceeds reference memory limit");
    MatType block=MatType::Zero(d*len,d*len);
    DataType prefactor=std::ldexp(1.,n)*((n*(len/2))%2 ? -1.:1.);
    for(int i=0;i<len;++i) {
        // Input lists pulses in chronological order, matching left multiplication.
        const auto& f=factors[len-1-i];
        block.block(i*d,i*d,d,d)=f.g;
        prefactor*=f.trace/std::ldexp(1.,n);
        for(int j=i+1;j<len;++j) {
            double sign=(j-i)%2 ? -1.:1.;
            block.block(i*d,j*d,d,d)=sign*id;
            block.block(j*d,i*d,d,d)=-sign*id;
        }
    }
    return prefactor*pfaf(block.rows()/2,block);
}
DataType traceProduct(const std::vector<MatType>& pulses,bool spacetime) {
    require(!pulses.empty(),"No pulses");
    std::vector<Gaussian> factors; factors.reserve(pulses.size());
    for(const auto& a:pulses) factors.push_back(gaussian(a));
    if(spacetime) return traceGaussianFactors(factors);
    int d=pulses[0].rows(), n=d/2;
    MatType id=MatType::Identity(d,d), g=MatType::Zero(d,d);
    UDT product(d); DataType w=std::ldexp(1.,n);
    for(size_t l=0;l<factors.size();++l) {
        MatType block=MatType::Zero(2*d,2*d);
        block.topLeftCorner(d,d)=factors[l].g; block.bottomRightCorner(d,d)=g;
        block.topRightCorner(d,d)=-id; block.bottomLeftCorner(d,d)=id;
        w*=factors[l].trace/std::ldexp(1.,n)*(n%2 ? -1.:1.)*pfaf(d,block);
        if(l+1<factors.size()) {
            product=factors[l].b*product;
            if(std::abs(w)<1e-10) return traceGaussianFactors(factors);
            product.onePlusInv(g); g-=id;
            if(!g.allFinite() || (g+g.transpose()).norm()>1e-7*std::max(1.,g.norm())) return traceGaussianFactors(factors);
            g=.5*(g-g.transpose()).eval();
        }
    }
    require(std::isfinite(std::abs(w)),"Nonfinite signed Pfaffian trace"); return w;
}
DataType exactTrace(const std::vector<MatType>& pulses) {
    int n=pulses[0].rows()/2; require(n<=10,"Exact Fock mode limit");
    MatType p=MatType::Identity(1<<n,1<<n);
    for(const auto& a:pulses) p=(matrixExp(fockMatrix(a))*p).eval();
    return p.trace();
}
DataType determinantTrace(const std::vector<MatType>& pulses) {
    int n=pulses[0].rows()/2; MatType b=MatType::Identity(n,n); DataType c=0.;
    for(const auto& a:pulses) {
        MatType h=numberMatrix(a); c-=.5*h.trace(); b=(matrixExp(h)*b).eval();
    }
    return std::exp(c)*(MatType::Identity(n,n)+b).determinant();
}
void Hamiltonian::validate() const {
    require(modes>0 && q0.rows()==2*modes && q0.cols()==2*modes,"Hamiltonian dimensions");
    auto check=[&](const MatType& a) {
        require(a.rows()==2*modes && a.cols()==2*modes,"Hamiltonian channel dimensions");
        require(a.allFinite() && (a+a.transpose()).norm()<1e-10*std::max(1.,a.norm()),"Non-skew generator");
        if(hermitian) require((a-a.adjoint()).norm()<1e-10*std::max(1.,a.norm()),"Non-Hermitian Hamiltonian channel");
        if(conserving) (void)numberMatrix(a);
    };
    check(q0);
    for(const auto& c:channels) { require(std::isfinite(c.coupling) && c.coupling>=0.,"HS coupling must be nonnegative"); check(c.a); }
}
MatType Hamiltonian::exactHamiltonian() const {
    require(modes<=10,"Exact Hamiltonian oracle is limited to 10 modes");
    MatType h=-fockMatrix(q0)+scalarEnergy*MatType::Identity(1<<modes,1<<modes);
    for(const auto& c:channels) { MatType o=fockMatrix(c.a); h-=.5*c.coupling*o*o; }
    for(const auto& gate:controlledPhases) {
        int i,j; double theta; std::tie(i,j,theta)=gate;
        cVecType phases=cVecType::Ones(1<<modes);
        for(int v=0;v<(1<<modes);++v) if((v&(1<<i)) && (v&(1<<j))) phases(v)=std::exp(1.i*theta);
        h=phases.asDiagonal()*h*phases.conjugate().asDiagonal();
    }
    return h;
}
int Hamiltonian::fieldCount(int slices) const { return slices*channels.size(); }
std::vector<MatType> Hamiltonian::pulses(int slices,double beta,const std::vector<double>& fields) const {
    require(slices>0 && std::isfinite(beta) && beta>0 && int(fields.size())==fieldCount(slices),"Invalid time grid/field count");
    double dt=beta/slices; std::vector<MatType> out; int idx=0;
    for(int l=0;l<slices;++l) {
        out.push_back(.5*dt*q0);
        for(const auto& c:channels) out.push_back(std::sqrt(dt*c.coupling)*fields[idx++]*c.a);
        out.push_back(.5*dt*q0);
    }
    return out;
}
FieldOperator::FieldOperator(const MatType& a,rdGenerator* r,bool variable):direction(a),field(variable),rng(r) {
    require((a-a.adjoint()).norm()<1e-10,"Dense HS walker requires Hermitian pulses");
    Eigen::SelfAdjointEigenSolver<MatType> es(a);
    require(es.info()==Eigen::Success,"HS eigensolver failed");
    vectors=es.eigenvectors(); eigenvalues=es.eigenvalues();
    setValue(field ? rng->rdNormal():1.);
}
void FieldOperator::setValue(double value) {
    phi=value;
    B=vectors*(phi*eigenvalues.array()).exp().matrix().asDiagonal()*vectors.adjoint();
    inverseB=vectors*(-phi*eigenvalues.array()).exp().matrix().asDiagonal()*vectors.adjoint();
    G=vectors*(-.5*phi*eigenvalues.array()).tanh().matrix().asDiagonal()*vectors.adjoint();
    G=.5*(G-G.transpose()).eval();
}
void FieldOperator::left_multiply(const MatType& a,MatType& out) { out=B*a; }
void FieldOperator::right_multiply(const MatType& a,MatType& out) { out=a*B; }
void FieldOperator::left_propagate(MatType& g,MatType& tmp) { tmp=B*g; g=tmp*inverseB; }
void FieldOperator::right_propagate(MatType& g,MatType& tmp) { tmp=inverseB*g; g=tmp*B; }
void FieldOperator::stabilizedLeftMultiply(UDT& f) { f=B*f; }
DataType FieldOperator::update(MatType& g) {
    if(!field) return 1.;
    ++attempts; double proposal=rng->rdNormal(), delta=proposal-phi;
    int d=g.rows(), n=d/2; MatType id=MatType::Identity(d,d);
    MatType gm=vectors*(-.5*delta*eigenvalues.array()).tanh().matrix().asDiagonal()*vectors.adjoint();
    gm=.5*(gm-gm.transpose()).eval();
    MatType block=MatType::Zero(2*d,2*d);
    block.topLeftCorner(d,d)=gm; block.bottomRightCorner(d,d)=g-id;
    block.topRightCorner(d,d)=-id; block.bottomLeftCorner(d,d)=id;
    double logEta=0.;
    for(int i=0;i<d;++i) logEta+=.5*std::log(std::cosh(.5*delta*eigenvalues(i)));
    DataType ratio=(n%2 ? -1.:1.)*std::exp(logEta)*pfaf(d,block);
    require(std::isfinite(std::abs(ratio)),"Nonfinite HS proposal ratio");
    if(rng->rdUniform01()>=std::min(1.,std::abs(ratio))) return 1.;
    MatType m=vectors*(delta*eigenvalues.array()).exp().matrix().asDiagonal()*vectors.adjoint();
    MatType system=id+(id-.5*g)*(m-id);
    g=system.partialPivLu().solve(g).eval();
    setValue(proposal); ++accepted;
    return ratio/std::abs(ratio);
}
Walker::Walker(const Hamiltonian& h,int slices,double beta,rdGenerator* r) {
    h.validate(); require(h.hermitian && slices>0 && beta>0,"Invalid HS walker model/time grid");
    nDim=2*h.modes; double dt=beta/slices;
    for(int l=0;l<slices;++l) {
        op_array.push_back(new FieldOperator(.5*dt*h.q0,r,false));
        for(const auto& c:h.channels) op_array.push_back(new FieldOperator(std::sqrt(dt*c.coupling)*c.a,r,true));
        op_array.push_back(new FieldOperator(.5*dt*h.q0,r,false));
    }
}
std::vector<double> Walker::fields() const {
    std::vector<double> f;
    for(auto op:op_array) { auto p=static_cast<FieldOperator*>(op); if(p->isField()) f.push_back(p->value()); }
    return f;
}
long Walker::attempts() const { long a=0; for(auto p:op_array) a+=static_cast<FieldOperator*>(p)->attempts; return a; }
long Walker::accepted() const { long a=0; for(auto p:op_array) a+=static_cast<FieldOperator*>(p)->accepted; return a; }
}
