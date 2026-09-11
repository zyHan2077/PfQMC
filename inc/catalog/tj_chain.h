#ifndef PFQMC_CATALOG_TJ_CHAIN_H
#define PFQMC_CATALOG_TJ_CHAIN_H
#include "quantum_link.h"
namespace catalog {
// Full Fock local t-J extension: doublons are static blockers, not discarded.
// H=-t sum tilde(c_i,s)^dag tilde(c_j,s)+h.c.
//   +J sum(S_i.S_j-n_i n_j/4)+U sum d_i+V sum n_i n_j-mu sum n_i.
// tilde(c_i,s)=(1-n_i,bar(s))c_i,s. Open nearest-neighbor chain is sign-free.
// Periodic rings are available solely to evaluate out-of-domain controls.
struct TJChain {
 int sites=4;bool periodic=false;
 double hopping=.8,exchange=.7,doublon=.23,density=.13,chemical=.17;
 int modes()const{return 2*sites;}
 void validate() const;
 bool signFreeParameters()const;
 std::vector<std::pair<int,int>> edges()const;
 std::vector<std::uint64_t> nextStates(std::uint64_t state)const;
 bool valid(const LinkPath& p)const;
 int squeezedPhase(std::uint64_t state)const;
 double diagonalEnergy(std::uint64_t state)const;
 double vertexScalar(std::uint64_t before,std::uint64_t after)const;
 MatType rotation(std::uint64_t before,std::uint64_t after)const;
 DataType weight(const LinkPath& p)const;
 DataType reference(const LinkPath& p)const;
 MatType exactHamiltonian()const;
};
}
#endif
