#ifndef CATALOG_SEARCH_MODELS_H
#define CATALOG_SEARCH_MODELS_H
#include "catalog.h"
namespace catalog {
// Candidate families: proof metadata provides independent small-matrix oracles.
struct SearchModel {
    enum Family { Lorentz, ConjugateTensor, SymmetricPower, QuaternionHermitian, QuaternionAntiHermitian, LorentzSymmetric } family;
    Hamiltonian h;
    int auxiliaryDim=0,degree=0,cells=0,cellModes=0;
    std::vector<MatType> auxiliary0;
    std::vector<double> chemical0;
    std::vector<std::vector<MatType>> auxiliaryChannels;
    std::vector<std::vector<double>> chemicalChannels;
    DataType reference(int slices,double beta,const std::vector<double>& fields) const;
};
MatType conjugateTensorAction(const MatType& k);
std::vector<std::vector<int>> symmetricStates(int d,int degree);
MatType symmetricPowerAction(const MatType& x,int degree);
std::vector<MatType> quaternionBasis(int n,bool antiHermitian=false);
MatType quaternionAction(const MatType& k,bool antiHermitian=false);
std::vector<MatType> lorentzSymmetricBasis(int n);
MatType lorentzSymmetricAction(const MatType& k);
SearchModel lorentzModel(int boosts=3,int lx=1,int ly=1);
SearchModel conjugateTensorModel(int d=3,int lx=1,int ly=1);
SearchModel symmetricPowerModel(int d=3,int degree=2,int lx=1,int ly=1);
SearchModel quaternionModel(int n=3,bool antiHermitian=false,int lx=1,int ly=1);
SearchModel lorentzSymmetricModel(int n=5,int lx=1,int ly=1);
std::vector<std::string> searchModelNames();
SearchModel makeSearchModel(const std::string& name);
}
#endif
