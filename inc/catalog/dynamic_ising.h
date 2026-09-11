#ifndef PFQMC_CATALOG_DYNAMIC_ISING_H
#define PFQMC_CATALOG_DYNAMIC_ISING_H
#include "catalog.h"
#include <cstdint>

namespace catalog {
enum class IsingVertex { DensityBond, SpinFlip };
struct IsingEvent {
    double time;
    IsingVertex kind;
    int index; // edge index for DensityBond; site index for SpinFlip
};
struct IsingConfiguration {
    double beta=1.2;
    std::vector<IsingEvent> events; // chronological, including a specified order at equal times
    // a,b occupations before the first vertex and after EVERY vertex.
    // A density vertex leaves them unchanged; a spin flip toggles 00 <-> 11.
    std::vector<std::uint64_t> spinPath;
};

// Three physical fermion orbitals c_i,a_i,b_i per site, without a projection:
// H = -t sum(c_i^dag c_j+h.c.) + V sum dn^c_i dn^c_j
//     -J/4 sum (n^a_i+n^b_i-1)(n^a_j+n^b_j-1)
//     +sum h_i/2 dn^c_i(a_i^dag b_i^dag+b_i a_i).
// V>=0; t,J,h_i real. Edges must form an undirected bipartite simple graph.
// Even a,b parity is spin 1/2; both odd states are static spin-zero blockers.
// This is a CT-INT configuration-weight API, not a production update algorithm.
struct DynamicIsing {
    int sites=2;
    std::vector<std::pair<int,int>> edges{{0,1}};
    double hopping=.8, repulsion=.5, ising=.7;
    std::vector<double> transverse{.6,-.4};

    static DynamicIsing square(); // four sites numbered around a cycle
    void validate() const;
    std::vector<int> sublattices() const; // BFS, including disconnected components
    // Invalid times, occupations, events, blockers or nonclosed paths return false.
    bool validConfiguration(const IsingConfiguration& configuration) const;
    // Full signed 3*sites-mode spacetime Pfaffian, with EVERY CT-INT scalar.
    DataType configurationPfaffian(const IsingConfiguration& configuration) const;
    // Exact c-Fock trace times the resolved spin path; independent sign oracle.
    DataType configurationReference(const IsingConfiguration& configuration) const;
    // Summed spin paths, using the c-Fock trace (sites<=5).
    DataType factorizedWordTrace(double beta,const std::vector<IsingEvent>& events) const;
    // Direct 3*sites-mode physical Hamiltonian and rotated CT-INT word (sites<=3).
    // Mode ordering is all c modes, then adjacent a_i,b_i pairs.
    MatType exactHamiltonian(bool staggeredSpinBasis=false) const;
    DataType exactWordTrace(double beta,const std::vector<IsingEvent>& events) const;

    // The following spin-only Gaussian factors act on 2*sites fermion modes.
    // Projectors have trace 1 and no invertible B. They are supported by
    // traceGaussianFactors, but must not be used by a regular Gaussian walker.
    Gaussian spinProjector(std::uint64_t occupation) const;
    Gaussian spinFlipRotation(int site,bool create) const;
    DataType spinPathPfaffian(const IsingConfiguration& configuration) const;
private:
    bool validEvents(double beta,const std::vector<IsingEvent>& events) const;
    double spinEnergy(std::uint64_t occupation) const;
    double spinPathAmplitude(const IsingConfiguration& configuration) const;
    MatType hoppingMatrix() const;
    DataType fermionWordTrace(double beta,const std::vector<IsingEvent>& events) const;
};
}
#endif
