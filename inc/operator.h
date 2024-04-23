#ifndef Operator_H
#define Operator_H

#include "types.h"
#include "qr_udt.h"
#include "skewMatUtils.h"

class Operator
{
public:
    // long long N;
    // long long para_rows;
    // long long para_cols;
    // MatType *para;
    virtual void left_multiply(const MatType &A, MatType &B){};
    virtual void inv_left_multiply(const MatType &A, MatType &B){};
    virtual void adjoint_left_multiply(const MatType &A, MatType &B){};
    virtual void right_multiply(const MatType &A, MatType &B){};
    virtual void inv_right_multiply(const MatType &A, MatType &B){};
    virtual void adjoint_inv_right_multiply(const MatType &A, MatType &B){};
    virtual void left_propagate(MatType &A, MatType &B){};
    virtual void right_propagate(MatType &A, MatType &B){};
    virtual void update(MatType &g){};
    virtual DataType getSignOfWeight() { return 1.0; };
    // virtual inline DataType getSignPfGInv() { return 1.0; };
    virtual DataType signOfUpdatedWeight(const MatType& g) {return 1.0; };
    virtual void getGreensMat(MatType &g0){};
    virtual void getGreensMatInv(MatType& g) {};

    // F = B * F
    virtual void stabilizedLeftMultiply(UDT &F){};

    virtual iVecType *getAuxField() { return NULL; };
    virtual int getType() { return -1; };
    virtual ~Operator(){};
};

class DenseOperator : public Operator
{
public:
    MatType mat;
    MatType mat_inv;
    MatType g0;
    MatType g0_inv;
    DataType signPf_g0_inv;
    DataType signOfWeight;
    DenseOperator(const MatType &mat_, DataType _s)
    {
        // this->N = Mat.rows();
        // this->para_cols = this->N;
        // this->para_rows = this->N;
        // this->para = &mat;
        mat = mat_;
        mat_inv = mat_.inverse();
        signOfWeight = _s;

        int nDim = mat.cols();
        g0 = (MatType::Identity(nDim, nDim) + mat).inverse() * 2.0 - MatType::Identity(nDim, nDim);
        // std::cout << "det of g0 = " << g0.determinant() <<"\n";
        g0_inv = g0.inverse();
        // std::cout << g0_inv << "===g0_inv======= \n\n";
        MatType tmp = g0_inv;
        signPf_g0_inv = signOfPfaf(tmp);
        // std::cout << g0_inv << "g0_inv \n";
    }

    void left_multiply(const MatType &A, MatType &B) override
    {
        B = mat * A;
    }
    void inv_left_multiply(const MatType &A, MatType &B) override
    {
        B = mat_inv * A;
    }
    void adjoint_left_multiply(const MatType &A, MatType &B) override
    {
        B = mat.adjoint() * A;
    }
    void right_multiply(const MatType &A, MatType &B) override
    {
        B = A * mat;
    }
    void inv_right_multiply(const MatType &A, MatType &B) override
    {
        B = A * mat_inv;
    }
    void adjoint_inv_right_multiply(const MatType &A, MatType &B) override
    {
        B = A * mat_inv.adjoint();
    }
    void left_propagate(MatType &A, MatType &B) override
    {
        B = mat * A;
        A = B * mat_inv;
    }

    void right_propagate(MatType &A, MatType &B) override
    {
        B = mat_inv * A;
        A = B * mat;
    }
    DataType getSignOfWeight() override
    {
        return signOfWeight;
    }

    // inline DataType getSignPfGInv() override {
    //     return signPf_g0_inv;
    // }

    void getGreensMat(MatType &g) override
    {
        g = g0;
    }

    void getGreensMatInv(MatType &g) override
    {
        g = g0_inv;
    }

    void stabilizedLeftMultiply(UDT &F) override
    {
        F = mat * F;
        // F.bMultUpdate(mat);
    }
};

// class SpinlessTvHoneycombUtils;

// class SpinlessVHoneycombOperator : public Operator
// {
// private:
//     const SpinlessTvHoneycombUtils *config;
//     const double etaM;
//     // delayed update for spinless t-V requires additional diagonalization
//     // const int delay_max = 32;
//     // cannot be larger than 64 (unless you change the threads value in kernel_linalg.cpp)
// public:
//     const int nUnitcell;
//     const int bondType;
//     // const int Naux; // number of auxillary fields (live on bonds)
//     const int nDim; // dimension of hamiltonian, number of sites × 2 (number of majorana species)
//     iVecType *s;
//     MatType B;
//     MatType B_inv;
//     rdGenerator *rd;

//     // _s: aux fields, Z_2 variable, length = nUnitcell
//     SpinlessVHoneycombOperator(const SpinlessTvHoneycombUtils *_config, iVecType *_s, int _bondType, rdGenerator *_rd);

//     ~SpinlessVHoneycombOperator() override;

//     void reCalcInv();

//     bool singleFlip(MatType &g, int idxCell, double rand) override;

//     void update(MatType &g) override;

//     void right_multiply(const MatType &AIn, MatType &AOut) override;

//     void left_multiply(const MatType &AIn, MatType &Aout) override
//     {
//         Aout = B * AIn;
//     }

//     void left_propagate(MatType &g, MatType &gTmp) override
//     {
//         reCalcInv();
//         gTmp = B * g;
//         g = gTmp * B_inv;
//     }
//     void right_propagate(MatType &g, MatType &gTmp) override
//     {
//         reCalcInv();
//         gTmp = B_inv * g;
//         g = gTmp * B;
//     }

//     iVecType *getAuxField() override
//     {
//         return s;
//     }

//     int getType() override
//     {
//         return bondType;
//     }

//     void getGreensMat(MatType &g) override;
//     void getGreensMatInv(MatType& g) override;
//     inline DataType getSignPfGInv() override;

//     void stabilizedLeftMultiply(UDT &F) override
//     {
//         // F.bMultUpdate(B);
//         F = B * F;
//     }
// };

#endif