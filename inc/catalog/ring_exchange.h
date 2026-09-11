#ifndef CATALOG_RING_EXCHANGE_H
#define CATALOG_RING_EXCHANGE_H
#include "quantum_link.h"
#include <array>
namespace catalog {
struct PairMove {std::vector<int> cells;double coupling;};
// b_r^dag=f_ra^dag f_rb^dag. Nearest-neighbor pair hopping plus
// -K sum_square (b_0^dag b_1 b_2^dag b_3+h.c.), J,K>=0, and real diagonal D.
// Full Fock odd occupations are static blockers. CT path weights retain signs
// outside the positive J,K domain to provide negative controls.
struct RingExchange {
 int lx=2,ly=2;double hopping=.4,ring=.3,repulsion=.17,chemical=.09;
 int cells() const{return lx*ly;}int modes() const{return 2*cells();}
 std::vector<std::pair<int,int>> edges() const;
 std::vector<std::array<int,4>> squares() const;
 std::vector<PairMove> moves() const;
 void validate() const;
 bool signFreeParameters() const{return hopping>=0 && ring>=0;}
 double diagonalEnergy(std::uint64_t state) const;
 std::vector<std::uint64_t> nextStates(std::uint64_t state) const;
 bool valid(const LinkPath& path) const;
 double moveCoupling(std::uint64_t before,std::uint64_t after) const;
 MatType rotation(std::uint64_t before,std::uint64_t after) const;
 DataType weight(const LinkPath& path) const;
 DataType reference(const LinkPath& path) const;
 MatType exactHamiltonian() const;
};
}
#endif
