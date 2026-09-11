#ifndef PFQMC_CATALOG_FREDKIN_MODELS_H
#define PFQMC_CATALOG_FREDKIN_MODELS_H
#include "jq_models.h"
namespace catalog {
struct FredkinTerm {int control,first,second;bool up;double coupling;};
// term=-1 denotes the sum of all off-diagonal terms for this transition.
struct FredkinEvent {double time;int term;};
struct FredkinPath {double beta=.8;std::vector<FredkinEvent> events;std::vector<std::uint64_t> states;};
// The same complete local two-complex-fermion spin representation as JQ:
// Sz=1/2-na, S+=(b^dag+b)a, S-=a^dag(b^dag+b), with two spin copies.
// H=sum_i[gL_i Pu_i P_(i+1,i+2)+gR_i P_(i,i+1) Pd_(i+2)]
//   +boundaryLeft Pd_0+boundaryRight Pu_last+sum_i field_i Sz_i.
// Pu/d=1/2+/-Sz, P_ij=1/4-Si.Sj. Highest fermionic order is six.
// The complete diagonal part is propagated exactly; CT vertices are strictly
// off diagonal. All nonnegative bulk couplings give nonnegative path weights.
// No projection to Dyck paths or local parity sectors is required in H.
struct FredkinModel {
    int sites=4;
    std::vector<FredkinTerm> terms{{0,1,2,true,1.},{1,2,3,true,1.},
                                  {2,0,1,false,1.},{3,1,2,false,1.}};
    std::vector<double> fields{0.,0.,0.,0.};
    double boundaryLeft=1.,boundaryRight=1.;
    int modes() const {return 2*sites;}
    static FredkinModel chain(int sites,double left=1.,double right=1.);
    void validate() const;
    bool signFreeParameters() const;
    std::uint64_t spinState(std::uint64_t state) const;
    std::vector<std::uint64_t> nextStates(std::uint64_t before,int term) const;
    double diagonalEnergy(std::uint64_t state) const;
    double transitionCoefficient(std::uint64_t before,std::uint64_t after,int term) const;
    bool valid(const FredkinPath& path) const;
    MatType vertexHalfRotation(std::uint64_t before,std::uint64_t after) const;
    DataType weight(const FredkinPath& path) const;
    DataType reference(const FredkinPath& path) const;
    // Direct products of independent onsite canonical-Fock or Pauli spins.
    MatType termMatrix(int term,bool fullFock=true) const;
    MatType exactHamiltonian(bool fullFock=true) const;
};
}
#endif
