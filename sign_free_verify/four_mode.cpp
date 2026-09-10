// Independent, signed validation of notes §3.2. No positive square root is used.
#include "skewMatUtils.h"
#include <array>
#include <fstream>
#include <iomanip>
#include <stdexcept>

using Params = std::array<double, 4>;
const MatType I16 = MatType::Identity(16, 16);
const MatType I8 = MatType::Identity(8, 8);
MatType comm(const MatType& a, const MatType& b) { return a*b-b*a; }
void require(bool ok, const std::string& message) {
    if (!ok) throw std::runtime_error(message);
}
MatType hexp(const MatType& q, double scale=1.) {
    Eigen::SelfAdjointEigenSolver<MatType> es(q);
    require(es.info() == Eigen::Success, "Hermitian eigensolver failed");
    return es.eigenvectors() * (scale*es.eigenvalues().array()).exp().matrix().asDiagonal()
           * es.eigenvectors().adjoint();
}

struct Model {
    std::array<MatType, 4> c, q, a;
    std::array<MatType, 8> gamma;
    Model() {
        for (int j=0; j<4; ++j) {
            c[j] = MatType::Zero(16,16);
            for (int v=0; v<16; ++v) if (v & (1<<j))
                c[j](v^(1<<j),v) = (__builtin_popcount(unsigned(v & ((1<<j)-1)))%2) ? -1.:1.;
            gamma[2*j] = c[j]+c[j].adjoint();
            gamma[2*j+1] = -1.i*(c[j]-c[j].adjoint());
        }
        q[0] = c[0].adjoint()*c[0]+c[1].adjoint()*c[1]-I16;
        q[1] = c[1].adjoint()*c[1]+c[2].adjoint()*c[2]-I16;
        MatType b = c[1].adjoint()*c[2]-c[3].adjoint()*(c[0]-c[0].adjoint());
        q[2] = b+b.adjoint();
        b = -c[0].adjoint()*c[1]+std::exp(-1.i*std::acos(-1.)/3.)*c[3].adjoint()*(c[2]-c[2].adjoint());
        q[3] = b+b.adjoint();
        for (int k=0; k<4; ++k) {
            require((q[k]-q[k].adjoint()).norm()<1e-13, "Non-Hermitian generator");
            a[k] = majorana(q[k]);
        }
        for (int i=0; i<8; ++i) for (int j=0; j<8; ++j)
            require((gamma[i]*gamma[j]+gamma[j]*gamma[i]-(i==j ? 2.*I16 : 0.*I16)).norm()<1e-13, "CAR failed");
        std::array<MatType,2> e, f, h;
        for (int i=0;i<2;++i) {
            e[i]=(q[i+2]+comm(q[i],q[i+2]))/2.;
            f[i]=(q[i+2]-comm(q[i],q[i+2]))/2.;
            h[i]=2.*q[i]-q[1-i];
        }
        double residual=comm(h[0],h[1]).norm();
        for (int i=0;i<2;++i) for (int j=0;j<2;++j) {
            double cartan=i==j ? 2.:-1.;
            residual=std::max(residual,(comm(h[i],e[j])-cartan*e[j]).norm());
            residual=std::max(residual,(comm(h[i],f[j])+cartan*f[j]).norm());
            residual=std::max(residual,(comm(e[i],f[j])-(i==j ? h[i] : MatType::Zero(16,16))).norm());
            if (i!=j) {
                residual=std::max(residual,comm(e[i],comm(e[i],e[j])).norm());
                residual=std::max(residual,comm(f[i],comm(f[i],f[j])).norm());
            }
        }
        require(residual<1e-12,"Chevalley-Serre relations failed");
        std::cout << "algebra_residual=" << residual << '\n';
    }
    MatType majorana(const MatType& op) const {
        // Q = gamma^T A gamma / 4, hence A_ab = -2 Tr(gamma_a gamma_b Q)/16.
        MatType result=MatType::Zero(8,8), reconstructed=MatType::Zero(16,16);
        for(int i=0;i<8;++i) for(int j=i+1;j<8;++j) {
            result(i,j)=-(gamma[i]*gamma[j]*op).trace()/8.;
            result(j,i)=-result(i,j);
            reconstructed+=.5*result(i,j)*gamma[i]*gamma[j];
        }
        require((reconstructed-op).norm()<1e-11,"Majorana normalization/reconstruction failed");
        return result;
    }
    MatType combine(const Params& p, const std::array<MatType,4>& basis) const {
        MatType out=MatType::Zero(basis[0].rows(),basis[0].cols());
        for(int i=0;i<4;++i) out+=p[i]*basis[i];
        return out;
    }
    MatType auxiliary(const Params& p) const {
        double u=p[0],v=p[1];
        MatType x=MatType::Zero(3,3);
        x(0,0)=(2*u+v)/3.; x(1,1)=(v-u)/3.; x(2,2)=(-u-2*v)/3.;
        x(0,1)=x(1,0)=p[2]; x(1,2)=x(2,1)=p[3];
        return x;
    }
};

struct Product {
    MatType fock=I16, b=I8;
    DataType weight=16., phase=1.;
    void append(const MatType& q, const MatType& a) {
        MatType bi=hexp(a), h=-a;
        MatType g=2.*(I8+b).inverse()-I8;
        MatType gi=2.*(I8+bi).inverse()-I8;
        // Eq. (5): eta12 = (-1)^N eta1 eta2 Pf([G1,-I;I,G2]). N=4.
        // pfaffianForEta in this repository returns the full trace, 2^N eta.
        MatType block=MatType::Zero(16,16);
        block.topLeftCorner(8,8)=gi; block.bottomRightCorner(8,8)=g;
        block.topRightCorner(8,8)=-I8; block.bottomLeftCorner(8,8)=I8;
        weight *= pfaffianForEta(h)*pfaf(8,block)/16.;
        phase *= pfaffianForSignOfEta(h)*pfaffianForSignOfProduct(gi,g);
        b=(bi*b).eval(); fock=(hexp(q)*fock).eval();
    }
};

int main(int argc,char** argv) try {
    int samples=1000; unsigned seed=240810311;
    std::string output;
    for(int i=1;i<argc;++i) {
        std::string arg=argv[i];
        require(i+1<argc,"Expected --samples N, --seed N, or --csv PATH");
        if(arg=="--samples") samples=std::stoi(argv[++i]);
        else if(arg=="--seed") seed=std::stoul(argv[++i]);
        else if(arg=="--csv") output=argv[++i];
        else throw std::runtime_error("Unknown argument: "+arg);
    }
    require(samples>0,"samples must be positive");
    mkl_set_num_threads(1);
    std::cout << std::setprecision(16);
    Model model;
    // Verify the Gaussian identity and its centered-density constant independently.
    MatType jacobi=MatType::Zero(64,64);
    for(int k=1;k<64;++k) jacobi(k-1,k)=jacobi(k,k-1)=std::sqrt(double(k));
    Eigen::SelfAdjointEigenSolver<MatType> quadrature(jacobi);
    require(quadrature.info()==Eigen::Success,"Gaussian quadrature failed");
    double hsError=0.;
    for(int i=0;i<2;++i) for(double coupling : {0.,.7,8.}) {
        MatType average=MatType::Zero(16,16);
        for(int k=0;k<64;++k)
            average+=std::norm(quadrature.eigenvectors()(0,k))*hexp(model.q[i],std::sqrt(coupling)*quadrature.eigenvalues()(k));
        MatType ni=model.c[i].adjoint()*model.c[i]-.5*I16;
        MatType nj=model.c[i+1].adjoint()*model.c[i+1]-.5*I16;
        MatType exact=hexp(ni*nj,coupling);
        hsError=std::max(hsError,(std::exp(-coupling/4.)*average-exact).norm()/exact.norm());
    }
    require(hsError<1e-11,"Gaussian HS/constant normalization failed");
    // Pairing contributes to the actual trace: delete all number-changing entries.
    MatType paired=model.combine({.2,-.15,.4,.3},model.q), conserving=paired;
    for(int i=0;i<16;++i) for(int j=0;j<16;++j)
        if(__builtin_popcount(unsigned(i))!=__builtin_popcount(unsigned(j))) conserving(i,j)=0.;
    double pairingEffect=std::abs(hexp(paired).trace()-hexp(conserving).trace());
    require(pairingEffect>1e-3,"Pairing unexpectedly does not affect the trace");
    std::cout << "HS_relative_error=" << hsError << " pairing_trace_difference=" << pairingEffect << '\n';
    std::mt19937 rng(seed); std::normal_distribution<double> normal;
    std::uniform_real_distribution<double> uniform(0.,1.);
    std::ofstream csv;
    if(!output.empty()) {
        csv.open(output); require(bool(csv),"Cannot open CSV");
        csv << std::setprecision(17)
            << "ensemble,sample,seed,slices,Z_fock_re,Z_fock_im,Z_character,Z_pf_re,Z_pf_im,phase_re,phase_im,positive_HS_factor\n";
    }
    double maxCharacter=0,maxPf=0,maxPhase=0,maxImag=0;
    long count=0;
    // Identity, commuting diagonal and one-edge boundaries, plus four random ensembles.
    for(int ensemble=0;ensemble<7;++ensemble) {
        double minimum=std::numeric_limits<double>::infinity();
        int n=ensemble<3 ? 1:samples;
        for(int trial=0;trial<n;++trial) {
            int length=ensemble<3 ? 4 : std::array<int,4>{1,8,32,64}[trial%4];
            Product product; MatType g=MatType::Identity(3,3), invg=g;
            auto append=[&](const Params& p) {
                product.append(model.combine(p,model.q),model.combine(p,model.a));
                MatType x=model.auxiliary(p);
                g=(hexp(x)*g).eval();
                // Inverse product by reversed factors, never by ill-conditioned inversion.
                invg=(invg*hexp(x,-1.)).eval();
            };
            double factor=1.;
            for(int l=0;l<length;++l) {
                if(ensemble==0) append({0,0,0,0});
                else if(ensemble==1) append({.2*(l-2),-.3*(l-1),0,0});
                else if(ensemble==2) append({.1,-.1,.2,0});
                else if(ensemble==3 || ensemble==4) {
                    double scale=(ensemble==3 ? .35:.8)/std::sqrt(double(length));
                    append({scale*normal(rng),scale*normal(rng),scale*uniform(rng),scale*uniform(rng)});
                } else {
                    // Symmetric Trotter for H=-Q0-U dn1 dn2-V dn2 dn3.
                    // Standard-normal independent fields; scalar e^{-beta(U+V)/4}.
                    double beta=ensemble==5 ? 1.:4., dt=beta/length;
                    double U=ensemble==5 ? .7:2., V=ensemble==5 ? 1.1:1.5;
                    Params half={dt*.2/2.,dt*(-.15)/2.,dt*.4/2.,dt*.3/2.};
                    append(half);
                    append({std::sqrt(dt*U)*normal(rng),std::sqrt(dt*V)*normal(rng),0,0});
                    append(half);
                    factor=std::exp(-beta*(U+V)/4.);
                }
            }
            DataType z=product.fock.trace();
            DataType character=2.*(g.trace()*invg.trace()-1.);
            double scale=std::max(16.,std::abs(z));
            double ec=std::abs(z-character)/scale, ep=std::abs(z-product.weight)/scale;
            double ei=std::abs(z.imag())/scale, es=std::abs(product.phase-1.);
            require(std::isfinite(std::abs(z)) && std::isfinite(ec) && std::isfinite(ep) && std::isfinite(es),"Nonfinite weight/error");
            if(csv) csv << ensemble << ',' << trial << ',' << seed << ',' << length << ',' << z.real() << ',' << z.imag() << ',' << character.real() << ',' << product.weight.real() << ',' << product.weight.imag() << ',' << product.phase.real() << ',' << product.phase.imag() << ',' << factor << '\n';
            require(z.real()>=16.-1e-9 && ei<1e-9,"Fock positivity/lower bound failed");
            require(ec<1e-8 && ep<1e-7 && es<1e-8,"Cross-check failed: ensemble="+std::to_string(ensemble)+" sample="+std::to_string(trial)+" character="+std::to_string(ec)+" pf="+std::to_string(ep));
            maxCharacter=std::max(maxCharacter,ec); maxPf=std::max(maxPf,ep);
            maxPhase=std::max(maxPhase,es); maxImag=std::max(maxImag,ei);
            minimum=std::min(minimum,z.real()); ++count;
        }
        std::cout << "ensemble=" << ensemble << " samples=" << n << " min_Z=" << minimum << '\n';
    }
    // Negative control outside the allowed family: unrestricted real hopping.
    // Search a reproducible finite set, then demand both Fock and Pfaffian detect - sign.
    bool negative=false;
    std::mt19937 controlRng(731); std::normal_distribution<double> controlNormal;
    for(int trial=0;trial<2000 && !negative;++trial) {
        Product product;
        for(int l=0;l<6;++l) {
            MatType h=MatType::Zero(4,4);
            for(int i=0;i<4;++i) for(int j=i;j<4;++j) h(i,j)=h(j,i)=.8*controlNormal(controlRng);
            MatType q=-.5*h.trace()*I16;
            for(int i=0;i<4;++i) for(int j=0;j<4;++j) q+=h(i,j)*model.c[i].adjoint()*model.c[j];
            product.append(q,model.majorana(q));
        }
        DataType z=product.fock.trace();
        if(z.real() < -1.) {
            require(std::abs(z-product.weight)/std::abs(z)<1e-7 && std::abs(product.phase+1.)<1e-8,"Negative control sign lost");
            std::cout << "negative_control_trial=" << trial << " Fock=" << z << " Pfaffian=" << product.weight << " phase=" << product.phase << '\n';
            negative=true;
        }
    }
    require(negative,"Negative control was not found");
    if(csv) { csv.flush(); require(bool(csv),"CSV write failed"); }
    std::cout << "PASS configurations=" << count << " seed=" << seed
              << " max_character_relative_error=" << maxCharacter
              << " max_pf_relative_error=" << maxPf << " max_phase_error=" << maxPhase
              << " max_relative_imag=" << maxImag << '\n';
    return 0;
} catch(const std::exception& e) {
    std::cerr << "FAIL: " << e.what() << '\n'; return 1;
}
