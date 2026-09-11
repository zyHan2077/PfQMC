#include "catalog/stabilizer_models.h"
#include <cmath>
#include <array>
#include <stdexcept>

namespace catalog {
namespace {
void check(bool q,const char* s) { if(!q) throw std::invalid_argument(s); }
int parity(uint64_t q) { return __builtin_popcountll(q)&1; }
struct Binary { uint64_t x=0,z=0; int phase=0; };
// Operator is i^phase X^x Z^z; multiplication is LEFT insertion.
Binary left(Binary a,Binary b) {
    return {a.x^b.x,a.z^b.z,(a.phase+b.phase+2*parity(a.z&b.x))%4};
}
Binary pauli(StabilizerPauli p) {
    uint64_t bit=uint64_t(1)<<p.qubit;
    if(p.axis==0) return {bit,0,0};
    if(p.axis==1) return {bit,bit,1};
    return {0,bit,0};
}
Binary vertex(const StabilizerVertex& v) {
    Binary out; out.phase=v.sign<0 ? 2:0;
    for(auto p:v.factors) out=left(pauli(p),out);
    return out;
}
DataType phase(int p) { static const DataType v[]={1.,1.i,-1.,-1.i}; return v[p%4]; }
std::pair<int,int> pair(StabilizerPauli p) {
    int b=4*p.qubit;
    return p.axis==0 ? std::make_pair(b+1,b+2) : p.axis==1 ? std::make_pair(b+2,b) : std::make_pair(b,b+1);
}
DataType globalPf(int qubits,const std::vector<StabilizerPauli>& word) {
    int d=4*qubits; double dim=std::ldexp(1.,2*qubits);
    if(word.empty()) return dim;
    std::vector<Gaussian> factors;
    for(auto p:word) {
        auto ij=pair(p); MatType a=MatType::Zero(d,d);
        a(ij.first,ij.second)=M_PI/2; a(ij.second,ij.first)=-M_PI/2;
        // The individual Pauli has trace zero. Each half-angle has nonzero
        // signed Gaussian trace; do NOT replace the final trace by sqrt(det).
        Gaussian half=gaussian(a); factors.push_back(half);factors.push_back(half);
    }
    return std::pow(DataType(0.,-1.),int(word.size()))*traceGaussianFactors(factors);
}
}
void StabilizerModel::validate() const {
    check(qubits>0 && qubits<=30,"Stabilizer qubit count must be 1..30");
    for(const auto& v:vertices) {
        check(std::isfinite(v.coupling)&&v.coupling>=0,"Stabilizer coupling must be nonnegative");
        check(v.sign==1||v.sign==-1,"Stabilizer sign must be +/-1"); uint64_t used=0;
        for(auto p:v.factors) {
            check(p.qubit>=0&&p.qubit<qubits&&p.axis>=0&&p.axis<3,"Invalid Pauli factor");
            uint64_t bit=uint64_t(1)<<p.qubit;check(!(used&bit),"Repeated site in a stabilizer");used|=bit;
        }
    }
}
bool StabilizerModel::positiveGroup() const {
    validate();std::vector<Binary> vs;
    for(const auto& v:vertices) if(v.coupling>0) vs.push_back(vertex(v));
    for(auto a:vs) for(auto b:vs) if(parity((a.x&b.z)^(a.z&b.x))) return false;
    std::array<Binary,60> basis{};std::array<bool,60> occupied{};
    for(auto a:vs) {
        for(int k=2*qubits-1;k>=0;--k) {
            uint64_t bits=a.x|(a.z<<qubits);
            if(!(bits&(uint64_t(1)<<k))) continue;
            if(occupied[k]) a=left(basis[k],a);
            else { basis[k]=a;occupied[k]=true;break; }
        }
        if(a.x==0&&a.z==0&&a.phase!=0) return false;
    }
    return true;
}
std::vector<StabilizerPauli> StabilizerModel::word(const std::vector<int>& labels) const {
    validate();std::vector<StabilizerPauli> out;
    for(int a:labels) { check(a>=0&&a<int(vertices.size()),"Invalid stabilizer label");
        out.insert(out.end(),vertices[a].factors.begin(),vertices[a].factors.end()); }
    return out;
}
double StabilizerModel::scalar(const std::vector<int>& labels,double beta) const {
    check(beta>=0&&std::isfinite(beta),"Invalid inverse temperature");double out=1.;int k=0;
    for(int a:labels) { check(a>=0&&a<int(vertices.size()),"Invalid stabilizer label");
        out *= beta*vertices[a].coupling*vertices[a].sign/(++k); }
    return out;
}
DataType stabilizerPauliPfaffian(int qubits,const std::vector<StabilizerPauli>& word,bool global) {
    check(qubits>0&&qubits<=30,"Invalid Pauli trace mode count");
    for(auto p:word) check(p.qubit>=0&&p.qubit<qubits&&p.axis>=0&&p.axis<3,"Invalid Pauli word");
    if(global) return globalPf(qubits,word);
    // Even local operators on disjoint cells commute. Factorize only by cell,
    // retaining every local time ordering and every -i scalar in its Pfaffian.
    DataType out=1.;
    for(int q=0;q<qubits;++q) {
        std::vector<StabilizerPauli> local;
        for(auto p:word) if(p.qubit==q) local.push_back({0,p.axis});
        out *= globalPf(1,local);
    }
    return out;
}
DataType StabilizerModel::pfaffian(const std::vector<int>& labels,double beta,bool global) const {
    auto w=word(labels);return scalar(labels,beta)*stabilizerPauliPfaffian(qubits,w,global);
}
DataType StabilizerModel::binaryOracle(const std::vector<int>& labels,double beta) const {
    Binary out;for(auto p:word(labels)) out=left(pauli(p),out);
    return (out.x||out.z) ? 0. : scalar(labels,beta)*std::ldexp(1.,modes())*phase(out.phase);
}
DataType stabilizerOccupationTrace(int qubits,const std::vector<StabilizerPauli>& word) {
    check(qubits>0&&qubits<=8,"Occupation oracle limited to 16 modes");DataType tr=0.;
    for(uint64_t start=0;start<(uint64_t(1)<<(2*qubits));++start) {
        uint64_t state=start;DataType amp=1.;
        for(auto p:word) {
            auto ij=pair(p);amp *= -1.i;
            for(int g:{ij.second,ij.first}) {
                int mode=g/2,occ=(state>>mode)&1;
                amp *= parity(state&((uint64_t(1)<<mode)-1)) ? -1.:1.;
                if(g%2) amp *= 1.i*double(1-2*occ);
                state ^= uint64_t(1)<<mode;
            }
        }
        if(state==start) tr+=amp;
    }
    return tr;
}
DataType StabilizerModel::occupationOracle(const std::vector<int>& labels,double beta) const {
    return scalar(labels,beta)*stabilizerOccupationTrace(qubits,word(labels));
}
MatType StabilizerModel::logicalHamiltonian() const {
    validate();check(qubits<=10,"Logical Hamiltonian limited to ten spins");
    int dim=1<<qubits;MatType h=MatType::Zero(dim,dim);
    for(const auto& v:vertices) {
        Binary a=vertex(v);
        for(int state=0;state<dim;++state)
            h(state^a.x,state) -= v.coupling*phase(a.phase)*(parity(a.z&state) ? -1.:1.);
    }
    return h;
}
MatType StabilizerModel::physicalHamiltonian() const {
    validate();check(qubits<=4,"Dense physical Hamiltonian limited to eight modes");
    int dim=1<<modes();MatType h=MatType::Zero(dim,dim);
    for(const auto& v:vertices) for(int start=0;start<dim;++start) {
        uint64_t state=start;DataType amp=-v.coupling*v.sign;
        for(auto p:v.factors) {
            auto ij=pair(p);amp *= -1.i;
            for(int g:{ij.second,ij.first}) {
                int mode=g/2,occ=(state>>mode)&1;
                amp *= parity(state&((uint64_t(1)<<mode)-1)) ? -1.:1.;
                if(g%2) amp *= 1.i*double(1-2*occ);
                state ^= uint64_t(1)<<mode;
            }
        }
        h(state,start)+=amp;
    }
    return h;
}
StabilizerModel toricStabilizers(int lx,int ly) {
    check(lx>=2&&ly>=2&&2*lx*ly<=30,"Toric lattice must be at least 2x2 and at most 30 edges");
    StabilizerModel m;m.qubits=2*lx*ly;
    auto edge=[=](int x,int y,int d){return 2*((x+lx)%lx+lx*((y+ly)%ly))+d;};
    for(int axis:{0,2}) for(int y=0;y<ly;++y) for(int x=0;x<lx;++x) {
        std::vector<int> sites=axis==0 ? std::vector<int>{edge(x,y,0),edge(x-1,y,0),edge(x,y,1),edge(x,y-1,1)}
             : std::vector<int>{edge(x,y,0),edge(x+1,y,1),edge(x,y+1,0),edge(x,y,1)};
        StabilizerVertex v;v.coupling=.2+.03*m.vertices.size();
        for(int q:sites) v.factors.push_back({q,axis});m.vertices.push_back(v);
    }
    m.validate();check(m.positiveGroup(),"Toric generator relation check failed");return m;
}
}
