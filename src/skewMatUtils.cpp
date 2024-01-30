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

// after the Householder transformation
// A is brought into tri-diagonalized form
// ... P_2 P_1 A P^T_1 P^T_2 ... =
//        0     k1  *  ... ...
//        -k1   0   k2  *  ...
//         *   -k2  0  k3  ...
//         ..............
// output Anew is a compact representation such that
//        Anew = 
//         0   k1   *   *  ...
//        u11   0  k2   *  ...
//        u12  u21  0   k3 ...
//        u13  u22 u31     ...
// P_i = 
//        I |                0
//        ------------------------------------ 
//        0 | I - 2 u \otimes u^\dagger / |u|^2
//
// pfaffian is evaluated as a side effect
//
// ===================================
//
// A should be a 2N * 2N matrix
// while size of temp be at least 4*N
void SkewMatHouseholder_PureMKL(const int N, dtype* A, dtype* temp, dtype* pfaf) {
    
    const MKL_INT L = 2*N;
    const MKL_INT inc=1;
    dtype* temp2 = &temp[L];
    const dtype two = 2;
    const char mklNoTrans = 'N';
    dtype norm2, arg0, k, alpha;
    dtype *x, *B, *subA;
    double normx0, sqrtnorm2, normu;
    dtype pfaffian = 1;

    for(uint i = 0; i<L-2; i++) {
        x = &(A[ind(i, i+1, L)]);
        MKL_INT len = L - i - 1;
        // zcopy(&len, &(A[ind(i, i+1, L)]), &inc, x, &inc);
        complexNorm2(x, len, &norm2);
        sqrtnorm2 = sqrt(norm2.real());
        // printVec(len, x);
        // std::cout << "sqrtnorm=" << sqrtnorm2 << "\n";

        // in the case norm2 = 0, pfaffian is just 0
        if (norm2.real() < thresholdDBL) {
            pfaffian = 0;
            std::cout << "pfaffian zero\n";
        } else {
            normx0 = std::norm(x[0]);
            arg0 = 1;
            if (normx0 > thresholdDBL) {
                arg0 = x[0]/sqrt(normx0);
            }

            k = arg0 * sqrtnorm2;
            x[0] += k;
            A[ind(i+1, i, L)] = -k;

            normu = norm2.real() - normx0 + std::norm(x[0]);
            // std::cout << "normu=" << normu << "\n";

            // subA is just the lower-right sub matrix of A
            subA = &(A[ind(i+1, i+1, L)]);
            // temp = u*
            vzConj(len, x, temp);
            // printVec(len, temp);
            // temp = subA . u*
            zgemv(&mklNoTrans, &len, &len, &one, subA, &L, temp, &inc, &zero, temp2, &inc);
            // printVec(len, temp2);
            // subA = subA - (2/|u|^2) . temp . u^T + (2/|u|^2) u . temp^T
            alpha = (-2.0)/normu;
            zgeru(&len, &len, &alpha, temp2, &inc, x, &inc, subA, &L);
            alpha = -alpha;
            zgeru(&len, &len, &alpha, x, &inc, temp2, &inc, subA, &L);
            
            // (2/|u|^2)^2 u . (u^\dagger . temp) . u^T \equiv 0
            // if(i%2 == 0) 
            pfaffian *=  double(i%2) + double(1-i%2) * k;

            // printMat(L, A);
        }
        // std::cout << "=========================loop end \n";
    }
    *pfaf = pfaffian * A[ind(L-1, L-2, L)];
}

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

// template <typename T>
// void printVec(uint L, T* x) {
//     std::cout << "===vec===\n";
//     for (int i = 0; i < L; i++) {
//         std::cout << x[i] << "\n";
//     }
//     std::cout << "\n===end===\n";
// }



// x^\dagger . x
void complexNorm2(const dtype* x, MKL_INT len, dtype* res) {
    MKL_INT inc = 1;
    zdotc(res, &len, x, &inc, x, &inc);
    // (*res) = sqrt(res->real());
}

dtype matDet(uint L, dtype* mat, lapack_int* temp) {
    int info = LAPACKE_zgetrf(LAPACK_COL_MAJOR, L, L, mat, L, temp);
    printMat(L, mat);
    if (info == 0) {
        dtype r = 1;
        for (int i=0; i<L; i++) {
            r *= mat[i*L + i];
        }
        return r;
    } else {
        std::cout << "\n===error===\n";
        return 0;
    }
}