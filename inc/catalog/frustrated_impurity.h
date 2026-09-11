#ifndef PFQMC_CATALOG_FRUSTRATED_IMPURITY_H
#define PFQMC_CATALOG_FRUSTRATED_IMPURITY_H
#include "catalog.h"

namespace catalog {
struct ImpurityTerm {
    std::vector<int> cells;
    MatType matrix; // one or two four-state spin cells, first cell most significant
};
struct ImpurityEvent { double time;int term; };
struct ImpurityConfiguration {
    double beta=1.2;
    std::vector<ImpurityEvent> events;
    // Each three-complex-mode cell has a conserved fermion parity +1 or -1.
    std::vector<int> parities;
    // Four-state cell labels before and after EVERY off-diagonal vertex.
    std::vector<std::vector<int>> path;
};

// Three open spin chains a=1,2,3, i=0,...,length-1, with a triangular root:
// H=lambda sum_a,i S_ai.S_a,i+1
//   +alpha[S_10.S_20+(S_10+S_20).S_30].
// S^a=-i epsilon_abc gamma^b gamma^c/4, {gamma,gamma}=2 delta.
// Every physical spin has THREE Majoranas, so H is explicitly local quartic
// in fermions. gamma_2r=f_r+f_r^dag, gamma_2r+1=i(f_r^dag-f_r).
// Pair spins (1,i),(2,i) into a three-complex-mode dimer cell; pair consecutive
// spins of chain 3 into the remaining cells. Even length avoids a dummy spin.
// The FULL Fock space equals 2^(3*length/2) copies of the physical spin space:
// there are no unphysical spin-zero blockers or required hard projections.
// lambda,alpha>=0 and equal couplings of the first two arms are the theorem
// domain. armImbalance changes arm 1 alone, for explicit negative controls.
// This is a signed CT-INT weight verifier, not a production updater.
struct FrustratedImpurity {
    int length=2;
    double chain=.8,impurity=1.1,armImbalance=0.;
    void validate() const;
    bool signFreeParameters() const;
    int spinSites() const {return 3*length;}
    int cells() const {return 3*length/2;}
    int modes() const {return 3*cells();}
    // dimerBasis=false exposes the ordinary Sz-basis sign problem as a control.
    std::vector<ImpurityTerm> terms(bool dimerBasis=true) const;
    bool validConfiguration(const ImpurityConfiguration& c,bool dimerBasis=true) const;
    double diagonalEnergy(const std::vector<int>& state,bool dimerBasis=true) const;
    DataType configurationPfaffian(const ImpurityConfiguration& c,bool dimerBasis=true) const;
    DataType configurationLocalReference(const ImpurityConfiguration& c,bool dimerBasis=true) const;
    // Independent full nine-mode Fock oracles (length=2).
    DataType configurationFockReference(const ImpurityConfiguration& c,bool dimerBasis=true) const;
    MatType exactHamiltonian() const;
    MatType exactSpinHamiltonian() const;
    // Gaussian projector onto every cell label and its parity, total trace 1.
    Gaussian stateProjector(const std::vector<int>& state,const std::vector<int>& parities,
                            bool dimerBasis=true) const;
};
}
#endif
