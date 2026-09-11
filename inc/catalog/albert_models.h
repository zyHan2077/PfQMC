#ifndef PFQMC_ALBERT_MODELS_H
#define PFQMC_ALBERT_MODELS_H
#include "catalog.h"
#include <array>

namespace catalog {
// Ordinary real coordinates; octonions never replace physical complex fermions.
using Octonion = std::array<double,8>;
using Albert = std::array<Octonion,9>;
Octonion octonionProduct(const Octonion&,const Octonion&);
Octonion octonionConjugate(const Octonion&);
Albert albertJordan(const Albert&,const Albert&);
Albert albertScale(const Albert&,double);
Albert albertAdd(const Albert&,const Albert&);
double albertInner(const Albert&,const Albert&);
double albertTrace(const Albert&);
double albertDeterminant(const Albert&);
const std::vector<Albert>& albertBasis();
const std::vector<Albert>& albertTracelessBasis();
MatType albertAction(const Albert&); // L_K(X)=K o X, in 27 real orbital coordinates
Albert albertApply(const MatType&,const Albert&);
// Independent Jordan quadratic-representation formula P(exp_J(K/2)).
MatType albertExponentialAction(const Albert&);

struct AlbertCandidate {
    Hamiltonian h;
    int cells;
    std::vector<Albert> kinetic;
    std::vector<double> mu;
    std::vector<std::vector<Albert>> channelK;
    std::vector<std::vector<double>> channelMu;
    DataType reference(int slices,double beta,const std::vector<double>& fields) const;
};
// Full 26 noncommuting Jordan directions plus a common density at EVERY cell.
// A fixed-size cell has 27 physical complex fermion orbitals.
AlbertCandidate makeAlbert(int cells=1);
}
#endif
