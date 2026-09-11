#ifndef CATALOG_POISSON_MODELS_H
#define CATALOG_POISSON_MODELS_H
#include "catalog.h"
namespace catalog {
struct PoissonVertex {
    int r=0,s=1;
    char axisR='z',axisS='z';
    double a=.31,b=-.27,kappa=.23;
};
// H = -sum_v kappa_v exp(a_v S^{axisR}_r + b_v S^{axisS}_s).
// Three orbitals per cell; S_x^3=S_x and S_z^3=S_z on the eight-state Fock
// space. Axis d denotes D=N-3/2, whose exponential is the product of three
// [cosh(a/2)+(2n_i-1)sinh(a/2)] factors. Both signs of onsite d vertices
// exclude a common Wei structure; onsite vertices are not a separate family.
// This is a CT expansion weight evaluator, not an update algorithm.
struct PoissonModel {
    int cells=2;
    std::vector<PoissonVertex> vertices;
    void validate() const;
    MatType generator(int vertex) const;
    MatType polynomialVertex(int vertex) const;
    MatType exactHamiltonian() const;
    // Unnormalized CT-series measure beta^k/k! prod kappa, including all
    // ordered-time integrals. A normalized Poisson proposal instead needs its
    // compensating exp(beta*sum kappa); no such proposal is used here.
    double orderedTimeScalar(const std::vector<int>& word,double beta) const;
    DataType signedWeight(const std::vector<int>& word,double beta) const;
    DataType fockWeight(const std::vector<int>& word,double beta) const;
};
PoissonModel poissonModel(int cells=2);
}
#endif
