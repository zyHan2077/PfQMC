#ifndef CATALOG_MERON_MODELS_H
#define CATALOG_MERON_MODELS_H
#include "catalog.h"
#include <cstdint>
namespace catalog {
enum class MeronKind { DensityBond, SingletBond, Monomer };
struct MeronEvent {double time;MeronKind kind;int index;};
struct MeronClusters {int clusters=0,merons=0,untouched=0,spinVertices=0;};
// Each site has c,a,b. E=1-na-nb+2na nb, S+=(a^dag b^dag), Sz=(na+nb-1)/2.
// P_i=E_i/2+Sx_i, P_ij=E_i E_j/4-S_i.S_j are rank-one projectors.
// H=H_tV-J sum P_ij-h sum sigma_i dn^c_i P_i, with V,J,h>=0.
// The complete Fock space includes two inert odd-parity states at each site.
// Meron grouping and a rotated pair-worldline representation describe ONE
// Hamiltonian construction. This API is a validator, not an update algorithm.
struct MeronModel {
    int sites=2;
    std::vector<std::pair<int,int>> edges{{0,1}};
    double hopping=.8,repulsion=.3,exchange=.6,coupling=.5;
    static MeronModel square();
    std::vector<int> sublattices() const;
    void validate() const;
    bool validEvents(double beta,const std::vector<MeronEvent>& events) const;
    MeronClusters clusters(const std::vector<MeronEvent>& events) const;
    double clusterTrace(const std::vector<MeronEvent>& events) const;
    DataType fermionPfaffian(double beta,const std::vector<MeronEvent>& events) const;
    // Density with respect to d tau_1 ... d tau_k on the ordered simplex.
    // All CT-INT signs and V/J/h factors are included; times are not integrated.
    DataType groupedWeight(double beta,const std::vector<MeronEvent>& events) const;
    DataType referenceWeight(double beta,const std::vector<MeronEvent>& events) const;
    // Each mask stores 2*sites spin-fermion occupations before/after every event.
    bool validPath(double beta,const std::vector<MeronEvent>& events,const std::vector<std::uint64_t>& path) const;
    DataType resolvedPfaffian(double beta,const std::vector<MeronEvent>& events,const std::vector<std::uint64_t>& path) const;
    DataType resolvedReference(double beta,const std::vector<MeronEvent>& events,const std::vector<std::uint64_t>& path) const;
    // Sum actual signed, full 3*sites-mode Pfaffians over all resolved spin paths.
    // Exponential small-system check, deliberately limited to sites<=2/events<=5.
    DataType exhaustivePfaffian(double beta,const std::vector<MeronEvent>& events,int* negativePaths=nullptr) const;
    MatType exactHamiltonian() const; // sites<=2, direct physical fermion operators
    DataType exactWordTrace(double beta,const std::vector<MeronEvent>& events) const;
};
}
#endif
