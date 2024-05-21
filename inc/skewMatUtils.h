#ifndef SkewMatUtils_H
#define SkewMatUtils_H

#include "types.h"
#include<iostream>
#include "../inc/pfapack/c_interface/pfapack.h"

// x^\dagger . x
void complexNorm2(const DataType* x, MKL_INT len, DataType* res);

/*
 * Householder is significantly slower than the Parlett-Reid approach
 * provided by the PFAPACK.
 * This routine is only for test purpose.
*/ 
void SkewMatHouseholder_PureMKL(const int N, DataType* A, DataType* temp, DataType* kVec);

DataType matDet(uint L, DataType* mat, lapack_int* temp);

DataType pfaf(const int N, MatType& A);

inline MatType expm(MatType &H, double lambda)
{
    Eigen::SelfAdjointEigenSolver<MatType> es(H);
    MatType V = es.eigenvectors();
    MatType D = es.eigenvalues();
    int N = H.rows();
    MatType expK(N, N);
    D = D.array() * lambda;
    D = D.array().exp();
    expK.noalias() = V * D.asDiagonal() * V.adjoint();
    return expK;
}

inline MatType sinhHQuarterSqrt2(MatType& H) {
    Eigen::SelfAdjointEigenSolver<MatType> es(H);
    MatType V = es.eigenvectors();
    MatType D = es.eigenvalues();
    int N = H.rows();
    MatType sinhK(N, N);
    D = D.array() * 0.25;
    MatType Dinv = - D.array();
    D = (D.array().exp() - Dinv.array().exp()) / sqrt(2);
    sinhK.noalias() = V * D.asDiagonal() * V.adjoint();
    return sinhK;
}

// calculate eta directly using
// eta = (-2)^N Pf [...]
// should only be used for testing
void generateMatForEta(const MatType& H, MatType& A);

// note that this calculation happens in place
DataType signOfPfaf(MatType& A);
DataType pfaffianForEta(const MatType &H);
DataType pfaffianForSignOfEta(const MatType &H);
DataType pfaffianForSignOfProduct(const MatType &G1, const MatType &G2 /*, bool diagno=false*/);

DataType signOfHamiltonian(const MatType &H);
DataType normalizeToPlusMinus1(DataType x);

#endif