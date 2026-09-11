#ifndef CATALOG_PROJECTED_MODELS_H
#define CATALOG_PROJECTED_MODELS_H
#include "search_models.h"
namespace catalog {
// Coefficients of Tr(z^N prod exp Q), evaluated with signed Gaussian Pfaffians.
// Fourier summands are complex intermediate values, never Monte Carlo weights.
std::vector<DataType> canonicalPfaffians(const std::vector<MatType>& pulses);
std::vector<DataType> canonicalFock(const std::vector<MatType>& pulses);
struct ProjectedModel {
    enum Kind { Charge, Parity } kind;
    Hamiltonian base;
    int cells=1, cellModes=6, target=2;
    double onsite=0.7, neighbor=0.2, bias=0.;
    // Empty: all sectors. Otherwise each entry is a required k (charge) or 0/1 parity.
    std::vector<int> hard;
    double sectorEnergy(const std::vector<int>& sector) const;
    bool allowed(const std::vector<int>& sector) const;
    std::vector<DataType> cellWeights(const std::vector<MatType>& pulses,int cell,bool fock=false) const;
    DataType weight(const std::vector<MatType>& pulses,const std::vector<int>& sector,double beta,bool fock=false) const;
    // Independent occupation-basis physical H in the requested tensor sector.
    // No Gaussian approximation to the positive charging/parity interaction.
    MatType exactSectorHamiltonian(const std::vector<int>& sector) const;
};
ProjectedModel chargingModel(int cells=1);
ProjectedModel parityModel(int cells=2);
}
#endif
