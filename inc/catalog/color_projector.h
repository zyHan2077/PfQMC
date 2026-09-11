#ifndef PFQMC_CATALOG_COLOR_PROJECTOR_H
#define PFQMC_CATALOG_COLOR_PROJECTOR_H
#include "catalog.h"
#include <cstdint>

namespace catalog {
enum class ColorVertex { Singlet, Permutation };
struct ColorEvent { double time; int edge; ColorVertex kind; };
struct ColorPath {
    double beta=.8;
    std::vector<ColorEvent> events;
    std::vector<std::uint64_t> states;
};
// Three ordinary fermions per cell. S^a=-i epsilon_abc f_b^dag f_c,
// A=N(3-N)/2. Full local Fock space is 0 + 1 + 1 + 0 under SO(3).
// P_ij=[(S_i.S_j)^2-A_i A_j]/3; Pi_ij=(S_i.S_j)^2+S_i.S_j-A_i A_j.
// H=D(N)-sum_edges[J P_ij+K Pi_ij], J,K>=0. Any finite graph is allowed.
// N=0,3 cells are retained spin-zero blockers; no Hilbert-space projection.
// This is a continuous-time signed configuration-weight API, not a sampler.
struct ColorProjector {
    int sites=3;
    std::vector<std::pair<int,int>> edges{{0,1},{1,2},{2,0}};
    double singlet=.9,permutation=.17,charging=.23,chemical=.11,density=.07,activeDensity=0.;
    int modes() const {return 3*sites;}
    void validate() const;
    bool signFreeParameters() const;
    static int color(unsigned localOccupation); // -1 for the two spin-zero states
    static int occupationPhase(unsigned localOccupation); // hole color y has phase -1
    static unsigned occupation(int number,int color);
    double diagonalEnergy(std::uint64_t state) const;
    double coefficient(ColorVertex kind) const;
    std::vector<std::uint64_t> nextStates(std::uint64_t before,int edge,ColorVertex kind) const;
    bool valid(const ColorPath& path) const;
    MatType vertexRotation(std::uint64_t before,std::uint64_t after) const;
    DataType weight(const ColorPath& path) const;
    DataType reference(const ColorPath& path) const;
    // Independent local canonical-Fock S matrices and polynomial vertex oracle.
    static MatType localSpin(int axis);
    static MatType localActive();
    static MatType localBondVertex(ColorVertex kind); // 64x64, cell i bits then cell j
    MatType vertexMatrix(int edge,ColorVertex kind) const; // full Fock, sites<=3
    MatType exactHamiltonian() const;
};
}
#endif
