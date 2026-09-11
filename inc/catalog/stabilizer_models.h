#ifndef PFQMC_CATALOG_STABILIZER_MODELS_H
#define PFQMC_CATALOG_STABILIZER_MODELS_H
#include "catalog.h"
#include <cstdint>

namespace catalog {
// Every qubit is represented on the COMPLETE Fock space of two local modes.
// X=-i gamma_1 gamma_2, Y=-i gamma_2 gamma_0, Z=-i gamma_0 gamma_1.
// The four-state cell is spin(1/2) times a two-dimensional spectator.
struct StabilizerPauli { int qubit=0, axis=0; }; // X=0,Y=1,Z=2
struct StabilizerVertex {
    std::vector<StabilizerPauli> factors;
    double coupling=1.;
    int sign=1;
};
struct StabilizerModel {
    int qubits=0;
    std::vector<StabilizerVertex> vertices;
    int modes() const { return 2*qubits; }
    void validate() const;
    bool positiveGroup() const;
    std::vector<StabilizerPauli> word(const std::vector<int>& labels) const;
    // Full ordered-time-integrated CT scalar beta^k/k! product(J_a sign_a).
    double scalar(const std::vector<int>& labels,double beta) const;
    DataType pfaffian(const std::vector<int>& labels,double beta=1.,bool global=false) const;
    DataType binaryOracle(const std::vector<int>& labels,double beta=1.) const;
    DataType occupationOracle(const std::vector<int>& labels,double beta=1.) const;
    MatType logicalHamiltonian() const;
    MatType physicalHamiltonian() const; // independent occupation action, <=8 modes
};
StabilizerModel toricStabilizers(int lx=2,int ly=2);
DataType stabilizerPauliPfaffian(int qubits,const std::vector<StabilizerPauli>& word,bool global=false);
DataType stabilizerOccupationTrace(int qubits,const std::vector<StabilizerPauli>& word);
}
#endif
