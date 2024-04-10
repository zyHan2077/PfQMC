#ifndef QR_UDT_H
#define QR_UDT_H

#include "types.h"

class UDT
{
public:
    MatType U;
    dVecType D;
    MatType T;

    UDT(int nDim) {
        U = MatType::Identity(nDim, nDim);
        D = dVecType::Ones(nDim);
        T = MatType::Identity(nDim, nDim);
    }

    UDT(const MatType& _U, const dVecType& _D, const MatType& _T) {
        U = _U;
        D = _D;
        T = _T;
    }

    // UDT& operator=(const UDT& other) {
    //     UDT F(other.U, other.D, other.T);
    //     return F;
    // }

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

    // F = (*this) * F
    inline void factorizedMultUpdate(UDT& F, int nDim) {
        MatType mat = T * F.U;
        mat = D.asDiagonal() * mat;
        mat = mat * F.D.asDiagonal();
        UDT tmp(mat, nDim);
        F.U = U * tmp.U;
        F.T = tmp.T * F.T;
        F.D = tmp.D;
    }

    // (*this) = B * (*this)
    inline void bMultUpdate(const MatType& B, int nDim) {
        // MatType r = ((B * U) * D.asDiagonal()) * T;
        MatType tmp = B * U;
        MatType tmp2 = tmp * D.asDiagonal();
        UDT F(tmp2, nDim);
        U = F.U;
        D = F.D;
        // std::cout << " r - UDTT " << (r - U * D.asDiagonal() * F.T * T).squaredNorm() << "\n";
        tmp = (F.T) * T;
        T = tmp;
    }

    // g = 2 * (1 + UDT)^{-1}
    inline void onePlusInv(int nDim, MatType& g) const {
        MatType Xinv = T;
        MatType tmp1, tmp2;
        // std::cout << Tn << "\n === \n";
        int ipiv[nDim];
        LAPACKE_zgetrf(LAPACK_COL_MAJOR, nDim, nDim, Xinv.data(), nDim, ipiv);
        // std::cout << "ipiv= " << ipiv << "\n";
        // std::cout << Tn << "\n";
        LAPACKE_zgetri(LAPACK_COL_MAJOR, nDim, Xinv.data(), nDim, ipiv);

        // std::cout << "Xinv * X - I = " << (Xinv*T - MatType::Identity(nDim, nDim)).squaredNorm() << "\n";

        dVecType Dpinv(nDim);
        dVecType Dm(nDim);
        for(int i=0; i<nDim; i++) {
            Dpinv(i) = 1.0 / std::max(D(i), 1.0);
            Dm(i) = std::min(D(i), 1.0);
        }

        tmp1 = Xinv * Dpinv.asDiagonal();
        tmp2 = U * Dm.asDiagonal();
        tmp1 = tmp1 + tmp2;

        UDT f = UDT(tmp1, nDim);
        LAPACKE_zgetrf(LAPACK_COL_MAJOR, nDim, nDim, f.T.data(), nDim, ipiv);
        LAPACKE_zgetri(LAPACK_COL_MAJOR, nDim, f.T.data(), nDim, ipiv);

        // std::cout << "Uinv * U - I = " << ((f.U)*(f.U.adjoint()) - MatType::Identity(nDim, nDim)).squaredNorm() << "\n";

        // MatType identity = MatType::Identity(nDim, nDim);
        // g = Xinv * Dpinv.asDiagonal() * f.T * (f.U * f.D.asDiagonal()).inverse();
        // return;

        tmp1 = (f.T) * (f.D.cwiseInverse()).asDiagonal();
        tmp2 = tmp1 * f.U.adjoint();
        tmp1 = Dpinv.asDiagonal() * tmp2;

        f = UDT(tmp1, nDim);

        f.D *= 2.0;

        tmp2 = Xinv * f.U;
        tmp1 = tmp2 * f.D.asDiagonal();
        g = tmp1 * f.T;
    }
};

#endif