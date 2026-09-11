#ifndef PFQMC_CATALOG_ZIGZAG_MODELS_H
#define PFQMC_CATALOG_ZIGZAG_MODELS_H
#include "catalog.h"

namespace catalog {
struct ZigzagTerm {std::vector<int> cells;MatType matrix;};
struct ZigzagEvent {double time;int term;};
struct ZigzagConfiguration {
    double beta=.9;
    std::vector<ZigzagEvent> events;
    std::vector<int> parities; // fixed +1/-1 fermion parity of every three-mode rung
    std::vector<std::vector<int>> path; // all rung labels before/after each vertex
};

// Open ladder of two spin-1/2 moments sigma_i,tau_i on each rung:
// H=sum_i[J0 sigma_i.sigma_(i+1)+J1 tau_i.tau_(i+1)
//         +J3 tau_i.sigma_(i+1)] + J2 sum_i sigma_i.tau_i.
// Use THREE Majoranas per spin, S^a=-i epsilon_abc gamma^b gamma^c/4,
// gamma_2r=f_r+f_r^dag, gamma_2r+1=i(f_r^dag-f_r). Thus H is explicitly
// local quartic in ordinary fermions; one rung contains three complex modes.
// The full Fock space is 2^rungs copies of the physical spin space, no blockers.
// J0,J1,J3>=0, J2 real. The sign-free cone is J3>=J0+J1. Each rung uses
// Bell/dimer-R states Phi+,Phi-,Psi+,Psi-. A nonlocal Kennedy-Tasaki phase
// makes the complete finite-temperature Hamiltonian stoquastic on an OPEN chain.
// We retain the original signed Bell-basis vertices in the Pfaffian path.
// Outside-cone parameters and the original Sz basis are supported as controls.
// This is a configuration-weight API, not a production Monte Carlo updater.
struct ZigzagLadder {
    int rungs=3;
    double J0=.24,J1=.31,J2=.8,J3=1.2;
    void validate() const;
    bool signFreeParameters() const;
    int spinSites() const {return 2*rungs;}
    int cells() const {return rungs;}
    int modes() const {return 3*rungs;}
    // First rungs-1 terms are inter-rung bonds, then rungs on-rung terms.
    std::vector<ZigzagTerm> terms(bool bellBasis=true) const;
    // O(rungs) evaluation of the EXACT nonlocal KT sign in the Bell basis.
    int ktPhase(const std::vector<int>& state) const;
    bool validConfiguration(const ZigzagConfiguration& c,bool bellBasis=true) const;
    double diagonalEnergy(const std::vector<int>& state,bool bellBasis=true) const;
    Gaussian stateProjector(const std::vector<int>& state,const std::vector<int>& parities,
                            bool bellBasis=true) const;
    DataType configurationPfaffian(const ZigzagConfiguration& c,bool bellBasis=true) const;
    DataType configurationLocalReference(const ZigzagConfiguration& c,bool bellBasis=true) const;
    // Independent original quartic Majorana-Fock oracles, rungs<=3.
    DataType configurationFockReference(const ZigzagConfiguration& c,bool bellBasis=true) const;
    MatType exactHamiltonian() const;
    MatType exactSpinHamiltonian() const; // rungs<=4
};
}
#endif
