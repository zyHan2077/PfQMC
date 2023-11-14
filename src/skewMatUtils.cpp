// Tempt to write skew matrix utilities, as
// what MKL has done for symmetric matrix
// not used yet

#include"skewMatUtils.h"

//  L*L skew matrix A,
//  y = A.x
//  only the upper or the lower triangle of A is used,
// terms are reorganized such that each independent 
// element is used only once
void SkewMatMulVec(char uplo, uint L, const dtype* A, const uint lda, const char transx, const dtype *x, dtype *y) {
    
    dtype temp;
    // clear y
    for (uint i=0; i<L; i++) y[i] = zero;

    uint aind = 0;
    if (uplo == 'U' ) {

        if(transx == 'N') {
            for (uint i=0; i<L; i++) {
                temp = zero;
                for (uint j=i+1; j<L; j++) {
                    y[j] -= A[aind+j] * x[i];
                    temp += A[aind+j] * x[j];
                }
                aind += lda;
                y[i] += temp;
            }
        } else if (transx == 'C') {
            for (uint i=0; i<L; i++) {
                temp = zero;
                for (uint j=i+1; j<L; j++) {
                    y[j] -= A[aind+j] * std::conj(x[i]);
                    temp += A[aind+j] * std::conj(x[j]);
                }
                aind += lda;
                y[i] += temp;
            }
        } else {
            printf("error!\n");    
        }
    } else {
        printf("not implemented ...\n");
    }
}