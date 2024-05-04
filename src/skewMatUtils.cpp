// no optimizations based on 
// skew-matrix properties has been 
// implemented. Operations are based
// entirely on MKL/BLAS dense matrix utilities

#include"skewMatUtils.h"
#include <fstream>

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
void SkewMatHouseholder_PureMKL(const int N, DataType* A, DataType* temp, DataType* kVec) {
    
    const MKL_INT L = 2*N;
    const MKL_INT inc = 1;
    DataType* temp2 = &temp[L];
    // const DataType two = 2;
    const char mklNoTrans = 'N';
    DataType norm2, arg0, k, alpha;
    DataType *x, *subA;
    double normx0, sqrtnorm2, normu;
    // DataType pfaffian = 1;
    int kcount = 0;

    for(uint i = 0; i<L-2; i++) {
        x = &(A[ind(i+1, i, L)]);
        MKL_INT len = L - i - 1;
        // zcopy(&len, &(A[ind(i+1, i, L)]), &inc, x, &inc);
        complexNorm2(x, len, &norm2);
        sqrtnorm2 = sqrt(norm2.real());
        // printVec(len, x);
        // std::cout << "sqrtnorm=" << sqrtnorm2 << "\n";

        // in the case norm2 = 0, pfaffian is just 0
        if ((i%2 == 0) && (norm2.real() < thresholdDBL)) {
            // pfaffian = 0;
            kVec[kcount] = 0;
            // std::cout << "pfaffian zero\n";
            return;
        } else {
            normx0 = std::norm(x[0]);
            arg0 = 1;
            if (normx0 > thresholdDBL) {
                arg0 = x[0]/sqrt(normx0);
            }

            k = arg0 * sqrtnorm2;
            x[0] += k;
            A[ind(i, i+1, L)] = k;

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
            // pfaffian *=  double(i%2) + double(1-i%2) * k;

            // printMat(L, A);
        }

        if (i%2 == 0) {
            kVec[kcount] = k;
            kcount++;
        }
        // std::cout << "=========================loop end \n";
    }
    kVec[kcount] = A[ind(L-2, L-1, L)];
    // *pfaf = pfaffian * A[ind(L-2, L-1, L)];
}

// Wrapper function for pfaffian calculation
// A is 2N*2N matrix, temp is a 4N vector
DataType pfaf(const int N, MatType& A) {
    DataType r = 1;
    cVecType kVec = cVecType::Zero(N);
    cVecType temp(4*N);
    SkewMatHouseholder_PureMKL(N, A.data(), temp.data(), kVec.data());
    for(int i=0; i<N; i++) {
        r *= kVec(i);
    }
    return r;
}

DataType signOfPfaf(MatType& A) {
    const int N = (A.cols() / 2);
    DataType r = 1.0;
    cVecType kVec = cVecType::Zero(N);
    cVecType temp(4*N);
    MatType B = A;
    SkewMatHouseholder_PureMKL(N, B.data(), temp.data(), kVec.data());
    for(int i=0; i<N; i++) {
        // std::cout << "kVec i=" << i << " " << kVec(i) << "\n";
        if (kVec(i) == 0.0) return 0.0; // TODO: DBL check?
        r *= kVec(i) / std::abs(kVec(i));
    }
    return r;
}


// x^\dagger . x
void complexNorm2(const DataType* x, MKL_INT len, DataType* res) {
    MKL_INT inc = 1;
    zdotc(res, &len, x, &inc, x, &inc);
}

void generateMatForEta(const MatType& H, MatType& A) {
    int n = H.cols() / 2;
    Eigen::SelfAdjointEigenSolver<MatType> es(H);
    MatType V = es.eigenvectors();
    cVecType D = es.eigenvalues();
    int N = H.rows();
    MatType sinhH(N, N);
    D = D * 0.25;
    cVecType expd = (D.array().exp() - (-D).array().exp())*0.5 * sqrt(2);
    sinhH.noalias() = V * expd.asDiagonal() * V.adjoint();

    A = MatType::Zero(4*n, 4*n);
    A.block(0, 0, 2*n, 2*n) = sinhH;
    A.block(2*n, 2*n, 2*n, 2*n) = sinhH;
    A.block(0, 2*n, 2*n, 2*n) = - MatType::Identity(2*n, 2*n);
    A.block(2*n, 0, 2*n, 2*n) = MatType::Identity(2*n, 2*n);
}


 DataType pfaffianForEta(const MatType &H) {
    int n = H.cols() / 2;
    MatType A;
    generateMatForEta(H, A);

    // cVecType tmp(8*n);
    DataType r = pfaf(2*n, A);
    r = r * pow(2, n);
    if (n%2 == 1) r = -r;
    return r;
}

 DataType pfaffianForSignOfEta(const MatType &H) {
    int n = H.cols() / 2;
    MatType A;
    generateMatForEta(H, A);

    // cVecType tmp(8*n);
    DataType r = signOfPfaf(A);
    if (n%2 == 1) r = -r;
    return r;
}

DataType pfaffianForSignOfProduct(const MatType &G1, const MatType &G2, bool diagno) {
    int n = G1.cols() / 2;
    MatType A = MatType::Zero(4*n, 4*n);
    A.block(0, 0, 2*n, 2*n) = G1;
    A.block(2*n, 2*n, 2*n, 2*n) = G2;
    A.block(0, 2*n, 2*n, 2*n) = - MatType::Identity(2*n, 2*n);
    A.block(2*n, 0, 2*n, 2*n) = MatType::Identity(2*n, 2*n);
    
    MatType Acopy = A;

    DataType r = signOfPfaf(A);

    if (diagno) {
        if (std::abs(r.real()) < 0.99 ) {
            std::cout << " ===== begin writing ==== \n";
            std::cout << "pf(Acopy) = " << pfaf(2*n, Acopy) << "\n";
            std::fstream myfile;
            myfile.open("1.dat", std::fstream::out);
            std::cout << " r = " << r << "\n";
            myfile << " r = " << r << "\n";
            for(int i=0; i<2*n; i++) {
                for(int j=0; j<2*n; j++) {
                    myfile << G1(i, j) << " ";
                    std::cout << ".";
                }
                myfile << "\n";
            }

            myfile << "\n\n\n";

            for(int i=0; i<2*n; i++) {
                for(int j=0; j<2*n; j++) {
                    myfile << G2(i, j) << " ";
                }
                myfile << "\n";
            }
            std::cout << " ===== end writing ==== \n";
            myfile.close();
        }
    }
    return r;
}

DataType signOfHamiltonian(const MatType &H) {
    int n = H.cols() / 2;
    MatType Hcopy = H;
    MatType sinhK = sinhHQuarterSqrt2(Hcopy);
    DataType sign = (n % 2 == 0) ? 1 : -1;
    sign *= pfaffianForSignOfProduct(sinhK, sinhK);
    return sign;
}