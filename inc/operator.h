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

class SpinlessV_honeycombOperator : public Operator
{
private:
    // delayed update for spinless t-V requires additional diagonalization
    // const int delay_max = 32; 
    // cannot be larger than 64 (unless you change the threads value in kernel_linalg.cpp)
public:
    int Nunitcell;
    int Naux; // number of auxillary fields (live on bonds)
    int Ndim; // dimension of hamiltonian, number of sites × 2 (number of majorana species)
    Eigen::MatrixXi *auxFields; // aux fields, 3 × Nunitcell, with value  +1 or -1
    MatType ham;
    MatType ham_inv;

    SpinlessV_honeycombOperator(double V, )


}