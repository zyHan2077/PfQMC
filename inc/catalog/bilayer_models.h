#ifndef PFQMC_CATALOG_BILAYER_MODELS_H
#define PFQMC_CATALOG_BILAYER_MODELS_H
#include "catalog.h"
namespace catalog {
struct BilayerTerm { std::vector<int> cells;MatType matrix; };
struct BilayerEvent { double time;int term; };
struct BilayerConfiguration {
    double beta=.9;
    std::vector<BilayerEvent> events;
    std::vector<int> parities;
    std::vector<std::vector<int>> path;
};
// One dimer per graph site, three physical complex modes per dimer.
// S^a_{2r+l}=-i gamma_{6r+3l+b} gamma_{6r+3l+c}/2 (a,b,c cyclic).
// H=sum_edges[Jz LzLz+Jperp(LxLx+LyLy)+Kz NzNz+Kperp(NxNx+NyNy)]
//   +sum_r[Dz S1z S2z+Dperp(S1x S2x+S1y S2y)], L=S1+S2,N=S1-S2.
// In layer couplings: scriptJ_mu=J_mu+K_mu, scriptK_mu=J_mu-K_mu.
// Uniform real couplings on any bipartite graph are sign-free if Kz=0 OR
// Kperp=0 OR Jperp=0. The fully frustrated conserved-dimer special case
// requires BOTH Kz=0 and Kperp=0; general members have dynamic dimers.
// Complete Fock = 2^cells faithful copies of the 2^(2*cells) spin space.
// CT path density includes every diagonal propagation and off-diagonal
// coefficient; times are NOT integrated, so no beta^k/k! is inserted.
struct AnisotropicBilayer {
    int sites=2;
    std::vector<std::pair<int,int>> edges{{0,1}};
    double Jz=.8,Jperp=.7,Kz=0.,Kperp=.35,Dz=.4,Dperp=1.1;
    int cells() const {return sites;}
    int modes() const {return 3*sites;}
    void validate() const;
    bool bipartite() const;
    bool signFreeParameters() const;
    std::vector<BilayerTerm> terms(bool dimerBasis=true) const;
    bool validConfiguration(const BilayerConfiguration& c,bool dimerBasis=true) const;
    double diagonalEnergy(const std::vector<int>& state,bool dimerBasis=true) const;
    Gaussian stateProjector(const std::vector<int>& state,const std::vector<int>& parities,bool dimerBasis=true) const;
    DataType configurationPfaffian(const BilayerConfiguration& c,bool dimerBasis=true) const;
    DataType configurationLocalReference(const BilayerConfiguration& c,bool dimerBasis=true) const;
    DataType configurationFockReference(const BilayerConfiguration& c,bool dimerBasis=true) const;
    MatType exactHamiltonian() const; // full physical Fock, up to two dimers
    MatType exactSpinHamiltonian(bool dimerBasis=true) const; // up to four dimers
};
AnisotropicBilayer squareBilayer(int face=0);
}
#endif
