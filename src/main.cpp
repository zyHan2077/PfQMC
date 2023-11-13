#include<iostream>
#include<complex>
#include "mkl.h"
// using namespace std;

typedef std::complex<double> dtype;
const dtype zero{0.0, 0.0};
#define ind(i, j, N) (i*N + j)

void printMat(uint N, dtype* A) {
    std::cout << "===skew mat===\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            std::cout << A[ind(i, j, N)] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "===end===\n";
}

void printVec(uint L, dtype* x) {
    std::cout << "===vec===\n";
    for (int i = 0; i < L; i++) {
        std::cout << x[i] << " ";
    }
    std::cout << "\n===end===\n";
}

// L*L skew matrix A, 
void SkewMatMulVec(char uplo, uint L, dtype* A, dtype *x, dtype *y) {
    // clear y
    for (uint i=0; i<L; i++) y[i] = zero;
    if (uplo == 'U' ) {
        for (uint i=0; i<L; i++) {
            dtype temp = zero;
            for (uint j=i+1; j<L; j++) {
                y[j] -= A[ind(i, j, L)] * x[i];
                temp += A[ind(i, j, L)] * x[j];
            }
            y[i] += temp;
        }
    } else {
        printf("working ...\n");
    }
}

int main() {
    int N = 2;
    dtype *A, *B, *x, *y;
    A = (dtype *) mkl_malloc((4*N*N)*sizeof(dtype), 64);
    B = (dtype *) mkl_malloc((4*N*N)*sizeof(dtype), 64);
    x = (dtype *) mkl_malloc((2*N)*sizeof(dtype), 64);
    y = (dtype *) mkl_malloc((2*N)*sizeof(dtype), 64);

    // initialize
    for(int i=0; i<4*N*N; i++) {
        A[i] = i;
    }
    
    for(int i=0; i<2*N; i++) {
        x[i] = i;
    }

    SkewMatMulVec('U', 2*N, A, x, y);

    printMat(2*N, A); printVec(2*N, x); printVec(2*N, y);

    return 0;
}