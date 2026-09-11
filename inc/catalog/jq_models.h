#ifndef PFQMC_CATALOG_JQ_MODELS_H
#define PFQMC_CATALOG_JQ_MODELS_H
#include "catalog.h"
#include <cstdint>
namespace catalog {
struct JQTerm {std::vector<std::pair<int,int>> bonds;double coupling;};
struct JQEvent {double time;int term;};
struct JQPath {double beta=.8;std::vector<JQEvent> events;std::vector<std::uint64_t> states;};
// Two complex fermions a,b per physical spin: Sz=1/2-na,
// S+=(b^dag+b)a, S-=a^dag(b^dag+b). Equivalently the first three of
// four onsite Majoranas give S^a=-i epsilon_abc gamma_b gamma_c/4.
// Every parity sector is one spin-1/2; FULL Fock Z=2^sites Z_spin.
// H=sum_i field_i Sz_i-sum_terms coupling*product_(ij in term) Pij,
// Pij=1/4-Si.Sj. Each multi-bond term is a disjoint matching of <=3 bonds.
// Nonnegative couplings and a bipartite union of bonds ensure positivity.
// The rectangular builder is finite range, including adjacent Q2/Q3 bonds.
// No production CT updater is supplied: this API calculates signed paths.
struct JQModel {
    int sites=4;
    std::vector<JQTerm> terms{{{{0,1}},.7},{{{1,2}},.7},{{{2,3}},.7},{{{3,0}},.7},
                             {{{0,1},{2,3}},1.1},{{{1,2},{3,0}},1.1}};
    std::vector<double> fields{0.,0.,0.,0.};
    int modes() const {return 2*sites;}
    static JQModel rectangle(int lx,int ly,int qOrder=2,double j=.7,double q=1.1);
    void validate() const;
    std::vector<int> sublattices() const; // empty for a nonbipartite control
    bool signFreeParameters() const;
    double diagonalEnergy(std::uint64_t state) const;
    std::uint64_t spinState(std::uint64_t occupation) const;
    static std::uint64_t occupation(std::uint64_t spinDown,std::uint64_t parityOdd,int sites);
    std::vector<std::uint64_t> nextStates(std::uint64_t before,int term) const;
    bool valid(const JQPath& path) const;
    MatType vertexHalfRotation(std::uint64_t before,std::uint64_t after,int term) const;
    DataType weight(const JQPath& path) const;
    DataType reference(const JQPath& path) const;
    // Independent canonical-Fock first-three-Majorana spin construction.
    static MatType localSpin(int axis);
    MatType vertexMatrix(int term) const; // full physical Fock, sites<=4
    MatType spinVertexMatrix(int term) const; // independent Pauli oracle
    MatType exactHamiltonian() const;
    MatType exactSpinHamiltonian() const;
};
}
#endif
