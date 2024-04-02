#ifndef SkewMatUtils_H
#define SkewMatUtils_H

#include "types.h"
#include<iostream>

// x^\dagger . x
void complexNorm2(const DataType* x, MKL_INT len, DataType* res);

void SkewMatHouseholder_PureMKL(const int N, DataType* A, DataType* temp, DataType* kVec);

DataType matDet(uint L, DataType* mat, lapack_int* temp);

DataType pfaf(const int N, MatType& A, cVecType& temp);

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

#endif