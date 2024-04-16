#ifndef QR_UDT_H
#define QR_UDT_H

#include "types.h"

class UDT
{
public:
    int nDim;
    MatType U;
    dVecType D;
    MatType T;
    UDT() = default;
    UDT(int nDim)
    {
        this->nDim = nDim;
        U = MatType::Identity(nDim, nDim);
        D = dVecType::Ones(nDim);
        T = MatType::Identity(nDim, nDim);
    }

    UDT(const MatType &_U, const dVecType &_D, const MatType &_T)
    {
        this->nDim = _U.rows();
        U = _U;
        D = _D;
        T = _T;
    }

    UDT &operator=(const UDT &other)
    {
        nDim = other.nDim;
        U = other.U;
        D = other.D;
        T = other.T;
        return *this;
    }

    UDT(const UDT &other)
    {
        *this = other;
    }

    UDT &operator=(UDT &&other)
    {
        if (this != &other)
        {
            nDim = other.nDim;
            U = std::move(other.U);
            D = std::move(other.D);
            T = std::move(other.T);
        }
    }

    UDT(UDT &&other)
    {
        *this = std::move(other);
    }

    // use qr to get UDT decomposition
    explicit UDT(MatType &A)
    {
        nDim = A.rows();
        T = MatType::Zero(nDim, nDim);
        D = dVecType(nDim);
        int jpvt[nDim];
        DataType tau[nDim];
        for (int j = 0; j < nDim; j++)
        {
            jpvt[j] = 0;
        }

        LAPACKE_zgeqp3(LAPACK_COL_MAJOR, nDim, nDim, A.data(), nDim, jpvt, tau);

        double alpha;
        for (int i = 0; i < nDim; i++)
        {
            D(i) = std::abs(A(i, i).real()); // A(i, i)'s are all real
            alpha = 1.0 / D(i);
            for (int j = i; j < nDim; j++)
            {
                A(i, j) = A(i, j) * alpha;
            }
        }

        int j;
        for (int i = 0; i < nDim; i++)
        {
            j = jpvt[i] - 1;
            T(i, j) = 1;
        }

        cblas_ztrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, nDim, nDim, &one, A.data(), nDim, T.data(), nDim);

        LAPACKE_zungqr(LAPACK_COL_MAJOR, nDim, nDim, nDim, A.data(), nDim, tau);

        U = A;
    }

    // // F = (*this) * F
    // inline void factorizedMultUpdate(UDT &F)
    // {
    //     MatType mat = T * F.U;
    //     mat = D.asDiagonal() * mat;
    //     mat = mat * F.D.asDiagonal();
    //     UDT tmp(mat);
    //     F.U = U * tmp.U;
    //     F.T = tmp.T * F.T;
    //     F.D = tmp.D;
    // }

    // // (*this) = B * (*this)
    // inline void bMultUpdate(const MatType &B)
    // {
    //     // MatType r = ((B * U) * D.asDiagonal()) * T;
    //     MatType tmp = B * U;
    //     MatType tmp2 = tmp * D.asDiagonal();
    //     UDT F(tmp2);
    //     U = F.U;
    //     D = F.D;
    //     // std::cout << " r - UDTT " << (r - U * D.asDiagonal() * F.T * T).squaredNorm() << "\n";
    //     tmp = (F.T) * T;
    //     T = tmp;
    // }

    // // g = 2 * (1 + UDT)^{-1}
    // inline void onePlusInv(MatType &g) const
    // {
    //     MatType Xinv = T;
    //     MatType tmp1, tmp2;
    //     // std::cout << Tn << "\n === \n";
    //     int ipiv[nDim];
    //     LAPACKE_zgetrf(LAPACK_COL_MAJOR, nDim, nDim, Xinv.data(), nDim, ipiv);
    //     // std::cout << "ipiv= " << ipiv << "\n";
    //     // std::cout << Tn << "\n";
    //     LAPACKE_zgetri(LAPACK_COL_MAJOR, nDim, Xinv.data(), nDim, ipiv);

    //     // std::cout << "Xinv * X - I = " << (Xinv*T - MatType::Identity(nDim, nDim)).squaredNorm() << "\n";

    //     dVecType Dpinv(nDim);
    //     dVecType Dm(nDim);
    //     for (int i = 0; i < nDim; i++)
    //     {
    //         Dpinv(i) = 1.0 / std::max(D(i), 1.0);
    //         Dm(i) = std::min(D(i), 1.0);
    //     }

    //     tmp1 = Xinv * Dpinv.asDiagonal();
    //     tmp2 = U * Dm.asDiagonal();
    //     tmp1 = tmp1 + tmp2;

    //     UDT f = UDT(tmp1);
    //     LAPACKE_zgetrf(LAPACK_COL_MAJOR, nDim, nDim, f.T.data(), nDim, ipiv);
    //     LAPACKE_zgetri(LAPACK_COL_MAJOR, nDim, f.T.data(), nDim, ipiv);

    //     // std::cout << "Uinv * U - I = " << ((f.U)*(f.U.adjoint()) - MatType::Identity(nDim, nDim)).squaredNorm() << "\n";

    //     // MatType identity = MatType::Identity(nDim, nDim);
    //     // g = Xinv * Dpinv.asDiagonal() * f.T * (f.U * f.D.asDiagonal()).inverse();
    //     // return;

    //     tmp1 = (f.T) * (f.D.cwiseInverse()).asDiagonal();
    //     tmp2 = tmp1 * f.U.adjoint();
    //     tmp1 = Dpinv.asDiagonal() * tmp2;

    //     f = UDT(tmp1);

    //     f.D *= 2.0;

    //     tmp2 = Xinv * f.U;
    //     tmp1 = tmp2 * f.D.asDiagonal();
    //     g = tmp1 * f.T;
    // }


    // g = 2 * (1 + UDT)^{-1}
    inline void onePlusInv(MatType &g) const
    {
        MatType Xinv = T.inverse();
        dVecType Dpinv(nDim);
        dVecType Dm(nDim);
        for (int i = 0; i < nDim; i++)
        {
            Dpinv(i) = 1.0 / std::max(D(i), 1.0);
            Dm(i) = std::min(D(i), 1.0);
        }
        MatType tmp1 = Xinv * Dpinv.asDiagonal();
        MatType tmp2 = (tmp1 + U * Dm.asDiagonal()).inverse();
        g = 2.0 * tmp1 * tmp2;
    }
};

inline UDT operator*(const UDT &udtL, const UDT &udtR)
{
    MatType mat = udtL.T * udtR.U;
    mat = udtL.D.asDiagonal() * mat;
    mat = mat * udtR.D.asDiagonal();
    UDT tmp(mat);
    tmp.U = udtL.U * tmp.U;
    tmp.T = tmp.T * udtR.T;
    return tmp;
}

// TO DO: lazy evaluation, check udta = B * udta
inline UDT operator*(const MatType &B, const UDT &udtR)
{
    MatType mat = B * udtR.U;
    mat = mat * udtR.D.asDiagonal();
    UDT tmp(mat);
    tmp.T = tmp.T * udtR.T;
    return tmp;
}

// return (1+udtR@(udtL).adjoint)^{-1}
inline MatType onePlusInv(UDT &udtL, UDT &udtR)
{
    int n = udtR.U.cols();
    MatType tem1 = udtR.U.adjoint() * udtL.U;
    MatType tem2 = udtR.T * udtL.T.adjoint();
    auto DrPinv = dVecType(n);
    auto DrM = dVecType(n);
    auto DlPinv = dVecType(n);
    auto DlM = dVecType(n);
    for (int i = 0; i < n; i++)
    {
        if (std::abs(udtR.D(i)) > 1.0)
        {
            DrPinv(i) = 1.0 / udtR.D(i);
            DrM(i) = 1.0;
        }
        else
        {
            DrPinv(i) = 1.0;
            DrM(i) = udtR.D(i);
        }
        if (std::abs(udtL.D(i)) > 1.0)
        {
            DlPinv(i) = 1.0 / udtL.D(i);
            DlM(i) = 1.0;
        }
        else
        {
            DlPinv(i) = 1.0;
            DlM(i) = udtL.D(i);
        }
    }

    tem2 = DrPinv.asDiagonal() * tem1 * DlPinv.asDiagonal() + DrM.asDiagonal() * tem2 * DlM.asDiagonal();
    tem1 = DlPinv.asDiagonal() * tem2.inverse() * DrPinv.asDiagonal();

    return 2. * udtL.U * tem1 * (udtR.U.adjoint());
}
#endif