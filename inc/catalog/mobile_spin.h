#ifndef PFQMC_CATALOG_MOBILE_SPIN_H
#define PFQMC_CATALOG_MOBILE_SPIN_H
#include "dynamic_ising.h"

namespace catalog {
enum class MobileVertex { DensityBond, SpinFlip, PairHop };
struct MobileEvent {
    double time;
    MobileVertex kind;
    int index; // site for SpinFlip, edge otherwise
};
struct MobileConfiguration {
    double beta=1.2;
    std::vector<MobileEvent> events;
    // Adjacent a_i,b_i occupations, before and after EVERY event.
    std::vector<std::uint64_t> spinPath;
};

// Three physical fermion modes c_i,a_i,b_i per site, full Fock space:
// D_i=n_ai+n_bi-1, p_i^dag=a_i^dag b_i^dag, dnci=n_ci-1/2.
// H=-t sum(c_i^dag c_j+h.c.)+V sum dnci dncj+Jz/4 sum D_i D_j
//   +Jxy/2 sum(p_i^dag p_j+h.c.)+sum_i h_i/2 dnci(p_i^dag+p_i).
// The graph must be bipartite. V,Jxy>=0, t,Jz real. Sign-free worldlines
// require h_i to have one sign on each connected hopping component; when
// Jxy=0 arbitrary real h_i are allowed (the DynamicIsing limit J=-Jz).
// Mixed-sign h is accepted by the weight API for physical negative controls;
// signFreeParameters() distinguishes it from the theorem's domain.
// This API verifies CT-INT weights, not production Monte Carlo updates.
struct MobileSpin {
    int sites=2;
    std::vector<std::pair<int,int>> edges{{0,1}};
    double hopping=.8,repulsion=.5,longitudinal=.7,exchange=.9;
    std::vector<double> transverse{.6,.4};
    static MobileSpin square();
    void validate() const;
    bool signFreeParameters() const;
    std::vector<int> sublattices() const;
    // Invalid times, transitions, blockers and nonclosed paths have zero weight.
    bool validConfiguration(const MobileConfiguration& c) const;
    MobileConfiguration resolve(double beta,const std::vector<MobileEvent>& events,
                                std::uint64_t initial) const;
    DataType configurationPfaffian(const MobileConfiguration& c) const;
    DataType configurationReference(const MobileConfiguration& c) const;
    DataType spinPathPfaffian(const MobileConfiguration& c) const;
    DataType factorizedWordTrace(double beta,const std::vector<MobileEvent>& events) const;
    // Independent physical Fock matrices, ordering all c then adjacent a,b pairs.
    MatType exactHamiltonian(bool staggeredSpinBasis=false) const; // sites<=3
    DataType exactWordTrace(double beta,const std::vector<MobileEvent>& events) const;
    Gaussian spinProjector(std::uint64_t occupation) const;
    Gaussian spinFlipRotation(int site,bool create) const;
private:
    DynamicIsing projectorModel() const;
    bool validEvents(double beta,const std::vector<MobileEvent>& events) const;
    double spinEnergy(std::uint64_t v) const;
    double spinPathAmplitude(const MobileConfiguration& c) const;
    DataType fermionWordTrace(double beta,const std::vector<MobileEvent>& events) const;
};
}
#endif
