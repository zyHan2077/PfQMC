#ifndef CATALOG_QUANTUM_LINK_H
#define CATALOG_QUANTUM_LINK_H
#include "catalog.h"
#include <cstdint>
namespace catalog {
struct LinkPath { double beta=.8; std::vector<double> times; std::vector<std::uint64_t> states; };
// One-dimensional path/cycle only. Each bond carries two auxiliary physical
// fermions a,b: B^dag=a^dag b^dag, D=n_a+n_b-1. No Fock-space projection.
// H=-t sum(c_j^dag B_e^dag c_i+h.c.) + V sum dn_i dn_j
//   +g/2 sum D_e(n_j-n_i) + mu sum dn_i.
// This is a continuous-time signed path-weight API, not a production sampler.
struct QuantumLink {
 int sites=3;bool periodic=true;
 double hopping=.7,repulsion=.4,electricDensity=.3,chemical=.12;
 std::vector<std::pair<int,int>> edges() const;
 int modes() const {return sites+2*int(edges().size());}
 void validate() const;
 double diagonalEnergy(std::uint64_t state) const;
 std::vector<std::uint64_t> nextStates(std::uint64_t state) const;
 bool valid(const LinkPath& path) const;
 MatType hopRotation(std::uint64_t before,std::uint64_t after) const;
 DataType weight(const LinkPath& path) const;
 DataType reference(const LinkPath& path) const;
 MatType exactHamiltonian() const;
};
Gaussian linkOccupationProjector(int modes,std::uint64_t state);
}
#endif
