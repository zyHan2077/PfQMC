#ifndef QR_UDT_H
#define QR_UDT_H

#include "types.h"

class UDT
{
public:
    MatType U;
    dVecType D;
    MatType T;
    UDT(MatType& _U, dVecType& _D, MatType& _T) {
        U = _U;
        D = _D;
        T = _T;
    }

    // use qr to get UDT decomposition
    UDT(MatType& A, const int nDim) {
        T = MatType::Zero(nDim, nDim);
        D = dVecType(nDim);
        int jpvt[nDim];
        DataType tau[nDim];
        for(int j=0; j<nDim; j++) {
            jpvt[j] = 0;
        }

        LAPACKE_zgeqp3(LAPACK_COL_MAJOR, nDim, nDim, A.data(), nDim, jpvt, tau);

        double alpha;
        for (int i=0; i<nDim; i++) {
            D(i) = std::abs(A(i, i).real()); // A(i, i)'s are all real
            alpha = 1.0 / D(i);
            for (int j = i; j < nDim; j++) {
			    A(i, j) = A(i, j) * alpha;
		    }
        }
        
        int j;
        for (int i=0; i<nDim; i++) {
            j = jpvt[i] - 1;
            T(i, j) = 1;
        }

        cblas_ztrmm(CblasColMajor, CblasLeft, CblasUpper,   CblasNoTrans, CblasNonUnit, nDim, nDim, &one, A.data(), nDim, T.data(), nDim);

        LAPACKE_zungqr(LAPACK_COL_MAJOR, nDim, nDim, nDim, A.data(), nDim, tau);

        U = A;
    }

    inline UDT factorizedMult(UDT& F, int nDim) {
        MatType mat = T * F.U;
        mat = D.asDiagonal() * mat;
        mat = mat * F.D.asDiagonal();
        UDT tmp(mat, nDim);
        tmp.U = U * tmp.U;
        tmp.T = tmp.T * F.T;
        return tmp;
    }
};

#endif