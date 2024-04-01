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

class SpinlesstvHoneycombUtils;

class SpinlessVOperator : public Operator
{
private:
    const SpinlesstvHoneycombUtils* config;
    const int bondType;
    const double etaM;
    std::uniform_real_distribution<double> uniform01;
    // delayed update for spinless t-V requires additional diagonalization
    // const int delay_max = 32; 
    // cannot be larger than 64 (unless you change the threads value in kernel_linalg.cpp)
public:
    const int nUnitcell;
    // const int Naux; // number of auxillary fields (live on bonds)
    const int nDim; // dimension of hamiltonian, number of sites × 2 (number of majorana species)
    iVecType* s;
    MatType B;
    MatType B_inv;

    // _s: aux fields, Z_2 variable, length = nUnitcell
    SpinlessVOperator(SpinlesstvHoneycombUtils* _config, iVecType* _s, int _bondType)
        :config(_config), bondType(_bondType), etaM(_config->etaM), 
            nDim(config->nsites * 2), nUnitcell(config->nUnitcell) {
        s = _s;
        uniform01 = std::uniform_real_distribution<double>(0.0, 1.0);

        B = MatType::Identity(nDim, nDim);
        config->InteractionBGenerator(B, *s, bondType);
    }

    void reCalc() {
        B = MatType::Identity(nDim, nDim);
        config->InteractionBGenerator(B, *s, bondType);
    }


    bool singleFlip(MatType &g, int idxCell, double rand) {
        DataType r = etaM;
        auto m = config->idxCell2Coord(idxCell);
        int auxCur = (*s)(idxCell);
        int idx1, idx2;
        DataType tmp[2];
        const int inc =  1;
        DataType alpha;
        

        for (int imaj = 0; imaj < 2; imaj ++) {
            idx1 = config->majoranaCoord2Idx(m.ix, m.iy, 0, imaj);
            idx2 = config->neighborSiteIdx(m.ix, m.iy, imaj, bondType);
            // tmp = [1 + i \sigma_{12} \tanh(\lambda / 2) G_{12}]
            tmp[imaj] = ( 1.0 + ( (1.0i) * (config->thHalflV) * double(auxCur) * g(idx1, idx2) ) );
            r *= tmp[imaj];
        }

        bool flag = rand < std::abs(r);

        if (flag) {
            for (int imaj = 0; imaj < 2; imaj ++) {
                idx1 = config->majoranaCoord2Idx(m.ix, m.iy, 0, imaj);
                idx2 = config->neighborSiteIdx(m.ix, m.iy, imaj, bondType);

                // update aux field and B matrix
                (*s)(idxCell) = -auxCur;
                B(idx1, idx2) = -B(idx1, idx2);
                B(idx2, idx1) = -B(idx2, idx1);

                // update Green's function
                cVecType x1 = B.col(idx1);
                cVecType x2 = B.col(idx2);
                alpha = (-1.0i) * (config->thHalflV) / tmp[imaj];
                zgeru(&nDim, &nDim, &alpha, x1.data(), &inc, x2.data(), &inc, B.data(), &nDim);
                alpha = -alpha;
                zgeru(&nDim, &nDim, &alpha, x2.data(), &inc, x1.data(), &inc, B.data(), &nDim);
            }
        }

        return flag;
    }

    void update(MatType &g) {
        for (int i=0; i<nUnitcell; i++) {
            // random real between (0, 1)
            double rand = uniform01(config->rdEng);
            bool flag = singleFlip(g, i, rand);
        }
    }
};