#ifndef PFQMC_CATALOG_PXP_MODELS_H
#define PFQMC_CATALOG_PXP_MODELS_H
#include "jq_models.h"
namespace catalog {
struct PXPEvent {double time;int site;};
struct PXPPath {double beta=.8;std::vector<PXPEvent> events;std::vector<std::uint64_t> states;};
// COMPLETE two-complex-mode spin cell, with no hard blockade projection:
// Z=1-2na, n=(1+Z)/2=1-na, P_down=na,
// X=a^dag b+b^dag a+a^dag b^dag+ba=-i gamma_1 gamma_2.
// H=-sum_i Omega_i P_down(i-1) X_i P_down(i+1)
//   +V sum_i n_i n_(i+1)+Delta sum_i n_i.
// Missing boundary controls are I. Each parity copy contains both spin states;
// the full physical partition function is 2^sites times the full spin trace.
// phaseFlux is ONLY an explicit out-of-domain test: multiply the site-0
// kinetic term by (1-2n_2), producing a frustrated closed four-flip loop.
struct PXPModel {
    int sites=4;
    std::vector<double> omega{.8,.8,.8,.8};
    double repulsion=.17,detuning=-.23;
    bool phaseFlux=false;
    int modes() const {return 2*sites;}
    static PXPModel chain(int sites,double omega=.8);
    void validate() const;
    bool signFreeParameters() const;
    static std::uint64_t occupation(std::uint64_t excitation,std::uint64_t odd,int sites);
    std::uint64_t excitations(std::uint64_t occupation) const;
    std::vector<int> enabledSites(std::uint64_t occupation) const;
    std::uint64_t flipped(std::uint64_t occupation,int site) const;
    double diagonalEnergy(std::uint64_t occupation) const;
    bool valid(const PXPPath& path) const;
    MatType halfRotation(int site,bool flux=false) const;
    // Density with respect to ordered event times: NO beta^k/k! factor here.
    DataType weight(const PXPPath& path) const;
    DataType reference(const PXPPath& path) const;
    // Complete kinetic SSE word, with beta^k/k! and every Omega/scalar.
    // Does not include the diagonal propagation; use weight() for full H.
    DataType kineticWordWeight(const std::vector<int>& labels,double beta=.8) const;
    MatType termMatrix(int site,bool fullFock=true) const;
    MatType exactHamiltonian(bool fullFock=true) const;
};
}
#endif
