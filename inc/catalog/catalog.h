#ifndef PFQMC_CATALOG_H
#define PFQMC_CATALOG_H
#include "pfqmc.h"
#include <functional>
#include <memory>
#include <string>
#include <vector>
#include <tuple>

namespace catalog {
void require(bool ok, const std::string& message);
MatType matrixExp(const MatType& a);
MatType numberMajorana(const MatType& h);
MatType numberMatrix(const MatType& a);
MatType pairMajorana(int modes, int i, int j, DataType delta, bool annihilate=true);
MatType fockMatrix(const MatType& a);
MatType occupationProjector(int modes, unsigned occupation);

struct Gaussian {
    MatType b, g;
    DataType trace;
};
Gaussian gaussian(const MatType& a);
// Actual signed trace: recursion for regular products, spacetime Pfaffian at zeros.
DataType traceProduct(const std::vector<MatType>& pulses, bool spacetime=false);
DataType traceGaussianFactors(const std::vector<Gaussian>& factors);
DataType exactTrace(const std::vector<MatType>& pulses);
DataType determinantTrace(const std::vector<MatType>& pulses);

struct Channel {
    std::string name;
    MatType a;
    double coupling;
};
struct Hamiltonian {
    std::string id, description;
    int modes;
    MatType q0; // H = -Q(q0) - sum lambda Q(a)^2/2 + scalarEnergy
    std::vector<Channel> channels;
    double scalarEnergy=0.;
    bool conserving=true;
    bool hermitian=true;
    // An optional fixed local non-Gaussian unitary: physical H = U H_base U^dagger.
    std::vector<std::tuple<int,int,double>> controlledPhases;
    void validate() const;
    MatType exactHamiltonian() const;
    std::vector<MatType> pulses(int slices, double beta, const std::vector<double>& fields) const;
    int fieldCount(int slices) const;
};

// A general, intentionally dense reference walker using PfQMC's stabilized sweeps.
// Gaussian independence proposals cancel the standard-normal HS prior exactly.
class FieldOperator : public Operator {
    MatType direction, vectors;
    Eigen::VectorXd eigenvalues;
    bool field;
    double phi=1.;
    rdGenerator* rng;
    void setValue(double value);
public:
    MatType B, inverseB, G;
    long attempts=0, accepted=0;
    FieldOperator(const MatType& direction, rdGenerator* rng, bool field);
    void left_multiply(const MatType& a,MatType& out) override;
    void right_multiply(const MatType& a,MatType& out) override;
    void left_propagate(MatType& g,MatType& tmp) override;
    void right_propagate(MatType& g,MatType& tmp) override;
    void stabilizedLeftMultiply(UDT& f) override;
    void getGreensMat(MatType& g) override { g=G; }
    DataType getSignOfWeight() override { return 1.; } // exp(Hermitian) has positive trace
    DataType update(MatType& g) override;
    double value() const { return phi; }
    bool isField() const { return field; }
};
class Walker : public Spinless_tV {
public:
    Walker(const Hamiltonian& model,int slices,double beta,rdGenerator* rng);
    std::vector<double> fields() const;
    long attempts() const;
    long accepted() const;
};

Hamiltonian openHopping(int modes=6);
Hamiltonian oddPairRing(int modes=5, int pairingEdges=1);
Hamiltonian fourMode();
Hamiltonian adjointOdd(int d=5);
Hamiltonian spinCell(int twiceJ=2, bool chemicalPotential=true);
Hamiltonian spinTensor(int jx=1,int jy=1);
Hamiltonian spinLattice(int j=1,int lx=2,int ly=2,bool mixed=false);
Hamiltonian fiveOrbital(int lx=2,int ly=2);
Hamiltonian pairingLattice(int lx=2,int ly=1);
Hamiltonian exteriorSL(int d=4,int k=2);
Hamiltonian exteriorSp(int n=4,int k=3);
Hamiltonian sl2Orthogonal(const std::vector<int>& multiplicities);
Hamiltonian transformedHopping(int modes=6);
Hamiltonian makeModel(const std::string& name);
std::vector<std::string> modelNames();
// Non-Hermitian trace-preserving extensions; not accepted by the Hermitian walker.
Hamiltonian makeGenerator(const std::string& name);
std::vector<std::string> generatorNames();
MatType wedgeAction(const MatType& x,int k);
MatType primitiveBasis(int n,int k);
MatType spinX(int twiceJ);
MatType spinZ(int twiceJ);

struct StaticModel {
    int lx=3, ly=3;
    double hopping=.3,U=3.,muF=.2,flux=.27,pairing=0.;
    MatType generator(const std::vector<int>& s) const;
    double logScalar(const std::vector<int>& s,double beta) const;
};
// Positive worldline representation of H=-t sum(b_i^dagger b_j+h.c.)+D(n).
struct PairHopping {
    int lx=2,ly=2;
    double hopping=.4,repulsion=.2,mu=.1;
    MatType exactHamiltonian() const;
    double diagonalEnergy(unsigned occupation) const;
    std::vector<std::pair<int,int>> edges() const;
    DataType closedHopPfaffian(unsigned before,unsigned after) const;
    // Signed Pfaffian check of a whole closed resolved path, including projectors.
    DataType pathPfaffian(const std::vector<unsigned>& path,double dt) const;
};
}
#endif
