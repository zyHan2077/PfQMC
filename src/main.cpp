#include"skewMatUtils.h"

// after the Householder transformation
// A is brought into tri-diagonalized form
// ... P_2 P_1 A P^T_1 P^T_2 ... =
//        0     k1  *  ... ...
//        -k1   0   k2  *  ...
//         *   -k2  0  k3  ...
//         ..............
// output Anew = 
//         0   u11 u12 ... ...
//        -k1   0  u21 u22 ...
//         *   -k2  0  u31 ...
//         ..............
// P_i = 
//        I |                0
//        ------------------------------------ 
//        0 | I - 2 u \otimes u^\dagger / |u|^2
//
// pfaffian is evaluated as a side effect
// void SkewMatHouseholder(char uplo, uint L, const dtype* A, dtype* Anew, dtype* pfaf) {
    
//     MKL_INT inc=1;
//     dtype norm2, arg0, *x, k;
//     double normx0, sqrtnorm2, normu;
//     dtype pfaffian = 1;

//     for(uint i = 0; i<L-2; i++) {
//         x = &(Anew[ind(i, i+1, L)]);
//         MKL_INT len = L - i - 1;
//         zcopy(&len, &(A[ind(i, i+1, L)]), &inc, x, &inc);
//         complexNorm2(x, len, &norm2);
//         sqrtnorm2 = sqrt(norm2.real());
//         // printVec(len, x);
//         // std::cout << "norm=" << norm2 << "\n";

//         // in the case norm2 = 0, pfaffian is just 0
//         if (norm2.real() < thresholdDBL) {
//             std::cout << "pfaffian zero\n";
//         } else {
//             normx0 = std::norm(x[0]);
//             arg0 = 1;
//             if (normx0 > thresholdDBL) {
//                 arg0 = x[0]/sqrt(normx0);
//             }

//             k = arg0 * sqrtnorm2;
//             x[0] += k;

//             normu = norm2.real() - normx0 + std::norm(x[0]);



//         // arg1 = x[0]/
//         }
//     }
// }

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
            // std::cout << "alpha=" << alpha << "\n";
            zgeru(&len, &len, &alpha, temp2, &inc, x, &inc, subA, &L);
            alpha = -alpha;
            zgeru(&len, &len, &alpha, x, &inc, temp2, &inc, subA, &L);
            
            // (2/|u|^2)^2 u . (u^\dagger . temp) . u^T \equiv 0
            if(i%2 == 0) pfaffian *= k;

            // printMat(L, A);
        }
        // std::cout << "=========================loop end \n";
    }
    *pfaf = pfaffian * A[ind(L-1, L-2, L)];
}

int main() {
    int N = 3;
    dtype *A, *temp1, *temp2;
    A = (dtype *) mkl_malloc((4*N*N)*sizeof(dtype), 64);
    temp1 = (dtype *) mkl_malloc((4*N)*sizeof(dtype), 64);
    temp2 = (dtype *) mkl_malloc((2*N)*sizeof(dtype), 64);
    // y = (dtype *) mkl_malloc((2*N)*sizeof(dtype), 64);

    // initialize
    int count = 1;
    for(int i=0; i<2*N; i++) {
        for (int j=i+1; j<2*N; j++) {
            A[i*2*N + j] = dtype{(double)count, -0.0};
            A[j*2*N + i] = dtype{-(double)count, 0.0};
            count++;
        }
        temp1[i] = count;
    }
    
    // printVec(2*N, temp1);
    // const MKL_INT inc=1;
    // // const dtype two = 2;
    // MKL_INT len = 2*N;
    // const char mklNoTrans = 'N';
    // zgemv(&mklNoTrans, &len, &len, &one, A, &len, temp1, &inc, &zero, temp2, &inc);
    // printVec(2*N, temp2);

    // std::cout << "smallest double: " << thresholdDBL << "\n";

    printMat(2*N, A);

    dtype pfaf;
    SkewMatHouseholder_PureMKL(N, A, temp1, &pfaf);
    
    printMat(2*N, A);

    std::cout << "pfaffian= " << pfaf << "\n";

    MKL_free(A);MKL_free(temp1);MKL_free(temp2);

    return 0;
}