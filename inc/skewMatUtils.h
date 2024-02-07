#ifndef SkewMatUtils_H
#define SkewMatUtils_H

#include "types.h"
#include<iostream>

// x^\dagger . x
void complexNorm2(const dtype* x, MKL_INT len, dtype* res);

void SkewMatHouseholder_PureMKL(const int N, dtype* A, dtype* temp, dtype* pfaf);

dtype matDet(uint L, dtype* mat, lapack_int* temp);

dtype pfaf(const int N, cMat& A, cVec& temp);

inline cMat expm(cMat &H, double lambda)
{
    Eigen::SelfAdjointEigenSolver<cMat> es(H);
    cMat V = es.eigenvectors();
    cMat D = es.eigenvalues();
    int N = H.rows();
    cMat expK(N, N);
    D = D.array() * lambda;
    D = D.array().exp();
    expK.noalias() = V * D.asDiagonal() * V.transpose();
    return expK;
}

#endif