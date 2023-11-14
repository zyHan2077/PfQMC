#include<complex>
#include<float.h>
#define MKL_Complex16 std::complex<double>
#include "mkl_types.h"
#include "mkl.h"
// using namespace std;

#include<iostream>

typedef std::complex<double> dtype;
const dtype zero{0.0, 0.0};
const dtype one{1.0, 0.0};
#define ind(i, j, N) (((i)*N) + j)
const double thresholdDBL = 1.0e-300;

void printMat(uint N, dtype* A) {
    std::cout << "===skew mat===\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            std::cout << A[ind(j, i, N)] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "===end===\n";
}

void printVec(uint L, dtype* x) {
    std::cout << "===vec===\n";
    for (int i = 0; i < L; i++) {
        std::cout << x[i] << "\n";
    }
    std::cout << "\n===end===\n";
}



// x^\dagger . x
void complexNorm2(const dtype* x, MKL_INT len, dtype* res) {
    MKL_INT inc = 1;
    zdotc(res, &len, x, &inc, x, &inc);
    // (*res) = sqrt(res->real());
}

void SkewMatMulVec(char uplo, uint L, const dtype* A, const uint lda, const char transx, const dtype *x, dtype *y);