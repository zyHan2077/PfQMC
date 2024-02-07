// no optimizations based on 
// skew-matrix properties has been 
// implemented. Operations are based
// entirely on MKL/BLAS dense matrix utilities

#include"skewMatUtils.h"

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
    const MKL_INT inc = 1;
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

// Wrapper function for pfaffian calculation
// A is 2N*2N matrix, temp is a 4N vector
dtype pfaf(const int N, cMat& A, cVec& temp) {
    dtype r;
    SkewMatHouseholder_PureMKL(N, A.data(), temp.data(), &r);
    return r;
}


// x^\dagger . x
void complexNorm2(const dtype* x, MKL_INT len, dtype* res) {
    MKL_INT inc = 1;
    zdotc(res, &len, x, &inc, x, &inc);
}

