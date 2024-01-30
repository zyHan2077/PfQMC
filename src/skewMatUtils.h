#ifndef SkewMatUtils_H
#define SkewMatUtils_H

#include "types.h"
// using namespace std;

#include<iostream>

void printMat(uint N, dtype* A);

template <typename T>
void printVec(uint L, T* x) {
    std::cout << "===vec===\n";
    for (int i = 0; i < L; i++) {
        std::cout << x[i] << "\n";
    }
    std::cout << "\n===end===\n";
}

// x^\dagger . x
void complexNorm2(const dtype* x, MKL_INT len, dtype* res);

void SkewMatMulVec(char uplo, uint L, const dtype* A, const uint lda, const char transx, const dtype *x, dtype *y);

void SkewMatHouseholder_PureMKL(const int N, dtype* A, dtype* temp, dtype* pfaf);

dtype matDet(uint L, dtype* mat, lapack_int* temp);

#endif