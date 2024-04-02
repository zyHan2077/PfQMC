#ifndef Operator_H
#define Operator_H

#include "types.h"

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
    virtual void update(MatType &d_gbar, MatType &d_gbaru, MatType &d_gbarv){};

    virtual iVecType* getAuxField(){return NULL;};
    virtual int getType(){return -1;};
};

class DenseOperator : public Operator
{
public:
    MatType mat;
    MatType mat_inv;
    DenseOperator(MatType &mat_)
    {
        // this->N = Mat.rows();
        // this->para_cols = this->N;
        // this->para_rows = this->N;
        // this->para = &mat;
        mat = mat_;
        mat_inv = mat_.inverse();
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
};

class SpinlessTvHoneycombUtils;

class SpinlessVOperator : public Operator
{
private:
    const SpinlessTvHoneycombUtils* config;
    const double etaM;
    // delayed update for spinless t-V requires additional diagonalization
    // const int delay_max = 32; 
    // cannot be larger than 64 (unless you change the threads value in kernel_linalg.cpp)
public:
    const int nUnitcell;
    const int bondType;
    // const int Naux; // number of auxillary fields (live on bonds)
    const int nDim; // dimension of hamiltonian, number of sites × 2 (number of majorana species)
    iVecType* s;
    MatType B;
    MatType B_inv;
    rdGenerator* rd;

    // _s: aux fields, Z_2 variable, length = nUnitcell
    SpinlessVOperator(const SpinlessTvHoneycombUtils* _config, iVecType* _s, int _bondType, rdGenerator* _rd);

    ~SpinlessVOperator() {
        delete s;
    }

    void reCalc();


    bool singleFlip(MatType &g, int idxCell, double rand);

    void update(MatType &g);

    void right_multiply(const MatType &AIn, MatType &AOut) override;

    void left_multiply(const MatType &AIn, MatType &Aout) override
    {
        Aout = B * AIn;
    }

    iVecType* getAuxField() override {
        return s;
    }

    int getType() override {
        return bondType;
    }
};

#endif