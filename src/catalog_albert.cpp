#include "catalog/albert_models.h"
#include <algorithm>
#include <cmath>

namespace catalog {
namespace {
// Cayley-Dickson convention (a,b)(c,d)=(ac-conj(d)b, da+b conj(c)).
// Recursive sizes are 1,2,4,8. All products below have only two factors.
void cdProduct(const double* x,const double* y,double* out,int n) {
    if(n==1){out[0]=x[0]*y[0];return;}
    int h=n/2;double ac[8]{},db[8]{},da[8]{},bc[8]{},dc[8]{},cc[8]{};
    for(int i=0;i<h;++i){dc[i]=(i?-1.:1.)*y[h+i];cc[i]=(i?-1.:1.)*y[i];}
    cdProduct(x,y,ac,h);cdProduct(dc,x+h,db,h);
    cdProduct(y+h,x,da,h);cdProduct(x+h,cc,bc,h);
    for(int i=0;i<h;++i){out[i]=ac[i]-db[i];out[h+i]=da[i]+bc[i];}
}
Albert unitAlbert(){Albert x{};for(int i=0;i<3;++i)x[3*i+i][0]=1.;return x;}
Albert albertProduct(const Albert& x,const Albert& y) {
    Albert z{};
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)for(int k=0;k<3;++k) {
        auto p=octonionProduct(x[3*i+k],y[3*k+j]);
        for(int a=0;a<8;++a)z[3*i+j][a]+=p[a];
    }
    return z;
}
Albert jordanExponential(Albert x) {
    double norm=std::sqrt(albertInner(x,x));int scale=0;
    if(norm>.5)scale=int(std::ceil(std::log2(norm/.5)));
    x=albertScale(x,std::ldexp(1.,-scale));
    Albert term=unitAlbert(),sum=term;
    for(int k=1;k<=80;++k) {
        term=albertScale(albertJordan(x,term),1./k);sum=albertAdd(sum,term);
        if(std::sqrt(albertInner(term,term))<2e-17)break;
        require(k<80,"Albert Jordan exponential failed to converge");
    }
    for(int i=0;i<scale;++i)sum=albertJordan(sum,sum);
    return sum;
}
MatType embedAlbert(const MatType& x,int cell,int cells) {
    MatType z=MatType::Zero(27*cells,27*cells);z.block(27*cell,27*cell,27,27)=x;return z;
}
}
Octonion octonionProduct(const Octonion& x,const Octonion& y) {
    Octonion z{};cdProduct(x.data(),y.data(),z.data(),8);return z;
}
Octonion octonionConjugate(const Octonion& x) {
    auto z=x;for(int i=1;i<8;++i)z[i]=-z[i];return z;
}
Albert albertJordan(const Albert& x,const Albert& y) {
    return albertScale(albertAdd(albertProduct(x,y),albertProduct(y,x)),.5);
}
Albert albertScale(const Albert& x,double a) {
    auto y=x;for(auto& q:y)for(auto& v:q)v*=a;return y;
}
Albert albertAdd(const Albert& x,const Albert& y) {
    auto z=x;for(int i=0;i<9;++i)for(int j=0;j<8;++j)z[i][j]+=y[i][j];return z;
}
double albertInner(const Albert& x,const Albert& y) {
    double s=0.;for(int i=0;i<9;++i)for(int j=0;j<8;++j)s+=x[i][j]*y[i][j];return s;
}
double albertTrace(const Albert& x){return x[0][0]+x[4][0]+x[8][0];}
double albertDeterminant(const Albert& x) {
    auto norm=[](const Octonion& q){double v=0.;for(double a:q)v+=a*a;return v;};
    double a=x[0][0],b=x[4][0],c=x[8][0];
    return a*b*c+2.*octonionProduct(octonionProduct(x[1],x[5]),x[6])[0]
        -a*norm(x[5])-b*norm(x[2])-c*norm(x[1]);
}
const std::vector<Albert>& albertBasis() {
    static const auto basis=[] {
        std::vector<Albert> ts;
        for(int i=0;i<3;++i){Albert t{};t[3*i+i][0]=1.;ts.push_back(t);}
        for(int i=0;i<3;++i)for(int j=i+1;j<3;++j)for(int a=0;a<8;++a) {
            Albert t{};t[3*i+j][a]=1./std::sqrt(2.);
            t[3*j+i][a]=(a?-1.:1.)/std::sqrt(2.);ts.push_back(t);
        }
        return ts;
    }();return basis;
}
const std::vector<Albert>& albertTracelessBasis() {
    static const auto basis=[] {
        auto& ts=albertBasis();std::vector<Albert> out;
        out.push_back(albertScale(albertAdd(ts[0],albertScale(ts[1],-1.)),1./std::sqrt(2.)));
        out.push_back(albertScale(albertAdd(albertAdd(ts[0],ts[1]),albertScale(ts[2],-2.)),1./std::sqrt(6.)));
        for(int i=3;i<27;++i)out.push_back(ts[i]);return out;
    }();return basis;
}
MatType albertAction(const Albert& k) {
    auto& ts=albertBasis();MatType h(27,27);
    for(int b=0;b<27;++b){auto y=albertJordan(k,ts[b]);for(int a=0;a<27;++a)h(a,b)=albertInner(ts[a],y);}
    return h;
}
Albert albertApply(const MatType& m,const Albert& x) {
    require(m.rows()==27 && m.cols()==27,"Albert action has wrong size");
    auto& ts=albertBasis();cVecType v(27);for(int a=0;a<27;++a)v(a)=albertInner(ts[a],x);
    v=(m*v).eval();require(v.imag().norm()<1e-9*std::max(1.,v.norm()),"Albert coordinates must remain real");
    Albert y{};for(int a=0;a<27;++a)y=albertAdd(y,albertScale(ts[a],v(a).real()));return y;
}
MatType albertExponentialAction(const Albert& k) {
    auto y=jordanExponential(albertScale(k,.5)),yy=albertJordan(y,y);
    auto& ts=albertBasis();MatType e(27,27);
    for(int b=0;b<27;++b) {
        auto z=albertAdd(albertScale(albertJordan(y,albertJordan(y,ts[b])),2.),albertScale(albertJordan(yy,ts[b]),-1.));
        for(int a=0;a<27;++a)e(a,b)=albertInner(ts[a],z);
    }
    return e;
}
AlbertCandidate makeAlbert(int cells) {
    require(cells>=1 && cells<=8,"Albert cells must be between one and eight");
    AlbertCandidate out;out.cells=cells;auto& h=out.h;
    h.id=cells==1?"albert27":"albert27_lattice";
    h.description="Full exceptional Albert Jordan cell with finite-range shared negative squares";
    h.modes=27*cells;h.q0=MatType::Zero(2*h.modes,2*h.modes);
    auto& directions=albertTracelessBasis();MatType kinetic=MatType::Zero(h.modes,h.modes);
    for(int r=0;r<cells;++r) {
        out.kinetic.push_back(albertAdd(albertScale(directions[0],.13),albertScale(directions[2+r%24],.09)));
        out.mu.push_back(r%2?-.07:.11);
        kinetic+=embedAlbert(albertAction(out.kinetic.back())+out.mu.back()*MatType::Identity(27,27),r,cells);
    }
    h.q0=numberMajorana(kinetic);
    auto add=[&](std::string name,double coupling,std::vector<Albert> ks,std::vector<double> mus) {
        MatType a=MatType::Zero(h.modes,h.modes);
        for(int r=0;r<cells;++r)a+=embedAlbert(albertAction(ks[r])+mus[r]*MatType::Identity(27,27),r,cells);
        h.channels.push_back({name,numberMajorana(a),coupling});out.channelK.push_back(std::move(ks));out.channelMu.push_back(std::move(mus));
    };
    for(int r=0;r<cells;++r) {
        std::vector<Albert> ks(cells);std::vector<double> mus(cells,0.);mus[r]=1.;
        add("density_"+std::to_string(r),.07,ks,mus);
        for(int a=0;a<26;++a) {
            ks.assign(cells,Albert{});mus.assign(cells,0.);ks[r]=directions[a];
            add("jordan_"+std::to_string(r)+"_"+std::to_string(a),.025,ks,mus);
        }
    }
    for(int r=0;r+1<cells;++r)for(int a=0;a<2;++a) {
        std::vector<Albert> ks(cells);std::vector<double> mus(cells,0.);
        ks[r]=directions[a?5:0];ks[r+1]=albertScale(directions[a?9:2],a?-1.:1.);
        add("local_exchange_"+std::to_string(r)+"_"+std::to_string(a),.035,ks,mus);
    }
    h.validate();return out;
}
DataType AlbertCandidate::reference(int slices,double beta,const std::vector<double>& fields) const {
    require(slices>0 && beta>0 && int(fields.size())==h.fieldCount(slices),"Albert reference fields");
    double dt=beta/slices;std::vector<MatType> b(cells,MatType::Identity(27,27)),half;
    std::vector<double> logR(cells,0.);
    for(int r=0;r<cells;++r)half.push_back(albertExponentialAction(albertScale(kinetic[r],dt/2.)));
    int index=0;
    for(int l=0;l<slices;++l) {
        for(int r=0;r<cells;++r){b[r]=(half[r]*b[r]).eval();logR[r]+=dt*mu[r]/2.;}
        for(int c=0;c<int(h.channels.size());++c) {
            double field=fields[index++]*std::sqrt(dt*h.channels[c].coupling);
            for(int r=0;r<cells;++r) {
                auto& k=channelK[c][r];
                if(albertInner(k,k)>0.)b[r]=(albertExponentialAction(albertScale(k,field))*b[r]).eval();
                logR[r]+=field*channelMu[c][r];
            }
        }
        for(int r=0;r<cells;++r){b[r]=(half[r]*b[r]).eval();logR[r]+=dt*mu[r]/2.;}
    }
    DataType w=1.;
    for(int r=0;r<cells;++r)w*=std::exp(-13.5*logR[r])*(MatType::Identity(27,27)+std::exp(logR[r])*b[r]).determinant();
    return w;
}
}
