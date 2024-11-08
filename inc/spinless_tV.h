#ifndef Spinless_tV_H
#define Spinless_tV_H

#include "operator.h"
#include "types.h"

class SpinlessTvUtils {
   public:
    int Lx, Ly;
    // int nsites;
    int nDim;
    // std::vector<int> nBond;
    double dt;
    double V;
    int l;
    double lambdaV, chlV, shlV, thlV, etaM;

    const int hsScheme;  

    const bool singleMaj;

    SpinlessTvUtils(int _Lx, int _Ly, double _dt, double _V, int _l, int _nDim,
                    bool _singleMaj = false, int _hsScheme = 0)
        : singleMaj(_singleMaj), hsScheme(_hsScheme) {
        // model configuration
        Lx = _Lx;
        Ly = _Ly;
        dt = _dt;
        V = _V;
        l = _l;  // imaginary time slices
        nDim = _nDim;

        lambdaV = acosh(exp(0.5 * V * dt));
        chlV = cosh(lambdaV);
        shlV = sinh(lambdaV);
        thlV = tanh(lambdaV);
        etaM = chlV * chlV;
    }

    inline int unitCellCoord2Idx(int ix, int iy) const { return ix * Ly + iy; }

    // virtual inline void KineticGenerator(MatType &H, DataType t) const {};

    // virtual inline DataType energyFromGreensFunc(const MatType &g) {};
    virtual inline void aux2MajoranaIdx(int idxAux, int imaj, int bType,
                                        int &idx1, int &idx2) const = 0;

    // Generate the Greens function for single slice
    inline void InteractionTanhGenerator(MatType &H, const iVecType &s,
                                         const int bondType,
                                         bool inv = false) const {
        DataType tmp = (1.0i) * tanh(0.5 * lambdaV);
        if (inv) {
            tmp = (-1.0) / tmp;
        }

        if (hsScheme == 0) {
            int idx1, idx2;
            if (!singleMaj) {
                for (int i = 0; i < s.size(); i++) {
                    for (int k = 0; k < 2; k++) {
                        aux2MajoranaIdx(i, k, bondType, idx1, idx2);
                        H(idx1, idx2) += -tmp * double(s(i));
                        H(idx2, idx1) += +tmp * double(s(i));
                    }
                }
            } else {
                for (int i = 0; i < s.size(); i++) {
                    aux2MajoranaIdx(i, 0, bondType, idx1, idx2);
                    H(idx1, idx2) += -tmp * double(s(i));
                    H(idx2, idx1) += +tmp * double(s(i));
                }
            }
        } else if (hsScheme == 1) {
            int idxi1, idxi2, idxj1, idxj2;
            assert(!singleMaj);
            for (int i = 0; i < s.size(); i++) {
                aux2MajoranaIdx(i, 0, bondType, idxi1, idxj1);
                aux2MajoranaIdx(i, 1, bondType, idxi2, idxj2);
                H(idxi1, idxi2) += -tmp * double(s(i));
                H(idxi2, idxi1) += +tmp * double(s(i));
                H(idxj1, idxj2) += +tmp * double(s(i));
                H(idxj2, idxj1) += -tmp * double(s(i));
            }
        }
    }

    // Directly generate B by directly writing each 2*2 block
    // B should be initialized as Identity
    inline void InteractionBGenerator(MatType &B, const iVecType &s,
                                      const int bondType,
                                      bool inv = false) const {
        DataType ch = chlV;
        DataType ish = (1.0i) * shlV;
        if (inv) ish = -ish;

        if (hsScheme == 0) {
            // B = exp(λ/2 * s * (i γi1 γj1 + i γi2 γj2))
            int idx1, idx2;
            if (!singleMaj) {
                for (int i = 0; i < s.size(); i++) {
                    for (int k = 0; k < 2; k++) {
                        aux2MajoranaIdx(i, k, bondType, idx1, idx2);
                        B(idx1, idx1) = ch;
                        B(idx2, idx2) = ch;
                        B(idx1, idx2) = +ish * double(s(i));
                        B(idx2, idx1) = -ish * double(s(i));
                    }
                }
            } else {
                for (int i = 0; i < s.size(); i++) {
                    aux2MajoranaIdx(i, 0, bondType, idx1, idx2);
                    B(idx1, idx1) = ch;
                    B(idx2, idx2) = ch;
                    B(idx1, idx2) = +ish * double(s(i));
                    B(idx2, idx1) = -ish * double(s(i));
                }
            }
        } else if (hsScheme == 1) {
            // B = exp(λ/2 * s * (i γi1 γi2 - i γj1 γj2))
            int idxi1, idxi2, idxj1, idxj2;
            assert(!singleMaj);
            for (int i = 0; i < s.size(); i++) {
                aux2MajoranaIdx(i, 0, bondType, idxi1, idxj1);
                aux2MajoranaIdx(i, 1, bondType, idxi2, idxj2);
                B(idxi1, idxi1) = ch;
                B(idxi2, idxi2) = ch;
                B(idxi1, idxi2) = +ish * double(s(i));
                B(idxi2, idxi1) = -ish * double(s(i));

                B(idxj1, idxj1) = ch;
                B(idxj2, idxj2) = ch;
                B(idxj1, idxj2) = -ish * double(s(i));
                B(idxj2, idxj1) = +ish * double(s(i));
            }
        }
    }
};

class SpinlessVOperator : public Operator {
   protected:
    const SpinlessTvUtils *config;
    const double etaM;
    // delayed update for spinless t-V requires additional diagonalization
    // const int delay_max = 32;
    // cannot be larger than 64 (unless you change the threads value in
    // kernel_linalg.cpp)
   public:
    // const int nUnitcell;
    const int bondType;
    // const int Naux; // number of auxillary fields (live on bonds)
    const int nDim;  // dimension of hamiltonian, number of sites × 2 (number of
                     // majorana species)

    const int hsScheme; // 0: hopping channel, 1: density channel
    iVecType *s;
    MatType B;
    MatType B_inv;
    rdGenerator *rd;

    const bool singleMaj;

    // _s: aux fields, Z_2 variable, length = nUnitcell
    SpinlessVOperator(const SpinlessTvUtils *_config, iVecType *_s,
                      int _bondType, rdGenerator *_rd)
        : config(_config),
          etaM(_config->etaM),
          bondType(_bondType),
          nDim(config->nDim),
          singleMaj(_config->singleMaj),
          hsScheme(_config->hsScheme) {
        s = _s;
        rd = _rd;
        B = MatType::Identity(nDim, nDim);
        config->InteractionBGenerator(B, *s, bondType, false);
    }

    ~SpinlessVOperator() { delete s; }

    void reCalcInv() {
        B_inv = MatType::Identity(nDim, nDim);
        config->InteractionBGenerator(B_inv, *s, bondType, true);
    }

    // virtual inline void aux2MajoranaIdx(int idxAux, int imaj, int& idx1, int&
    // idx2) {};

    void singleFlip(MatType &g, int idxAux, double rand, bool &flag,
                    DataType &signCur) {
        DataType r;
        // auto m = mConfig->idxCell2Coord(idxCell);
        int auxCur = (*s)(idxAux);
        int idx1, idx2, idx3, idx4;
        DataType tmp[2];
        const int inc = 1;
        DataType alpha;

        config->aux2MajoranaIdx(idxAux, 0, bondType, idx1, idx2); // 1j, 1k
        config->aux2MajoranaIdx(idxAux, 1, bondType, idx3, idx4); // 2j, 2k

        if (hsScheme == 0) {
            tmp[0] =
                (1.0 - ((1.0i) * (config->thlV) * double(auxCur) * g(idx1, idx2)));
            tmp[1] =
                (1.0 - ((1.0i) * (config->thlV) * double(auxCur) * g(idx3, idx4)));
            r = tmp[0] * tmp[1];
            // std::cout << "r1" << r << " ";
            r += (config->thlV * config->thlV) * ((g(idx1, idx3) * g(idx2, idx4)) -
                                                (g(idx2, idx3) * g(idx1, idx4)));
            // std::cout << " r2" << r << "\n";
            r *= etaM;
        } else if (hsScheme == 1) {
            tmp[0] =
                (1.0 - ((1.0i) * (config->thlV) * double(auxCur) * g(idx1, idx3)));
            tmp[1] =
                (1.0 + ((1.0i) * (config->thlV) * double(auxCur) * g(idx2, idx4)));
            r = tmp[0] * tmp[1];
            // std::cout << "r1" << r << " ";
            r -= (config->thlV * config->thlV) * ((g(idx1, idx2) * g(idx3, idx4)) -
                                                (g(idx3, idx2) * g(idx1, idx4)));
            // std::cout << " r2" << r << "\n";
            r *= etaM;
        }
        // for (int imaj = 0; imaj < 2; imaj ++) {
        //     config->aux2MajoranaIdx(idxAux, imaj, bondType, idx1, idx2);
        //     // idx1 = mConfig->majoranaCoord2Idx(m.ix, m.iy, 0, imaj);
        //     // idx2 = mConfig->neighborSiteIdx(m.ix, m.iy, imaj, bondType);
        //     // tmp = [1 + i \sigma_{12} \tanh(\lambda / 2) G_{12}]
        //     tmp[imaj] = ( 1.0 - ( (1.0i) * (config->thlV) * double(auxCur) *
        //     g(idx1, idx2) ) ); r *= tmp[imaj];
        // }

        flag = rand < std::abs(r);
        // std::cout << rand << "=rand " << "r = " << r << "\n";

        if (flag) {
            // DataType t = (r / std::abs(r));
            // if (std::abs(t.imag()) > 1e-5) {
            //     std::cout << "r= " << r << " sign(r)= " << t << "\n";
            // }
            signCur *= (r / std::abs(r));

            if (hsScheme == 0) {
                for (int imaj = 0; imaj < 2; imaj++) {
                    config->aux2MajoranaIdx(idxAux, imaj, bondType, idx1, idx2);
                    if (imaj == 1) {
                        tmp[1] = (1.0 - ((1.0i) * (config->thlV) * double(auxCur) *
                                        g(idx1, idx2)));
                    }
                    // update aux field and B matrix
                    (*s)(idxAux) = -auxCur;
                    B(idx1, idx2) = -B(idx1, idx2);
                    B(idx2, idx1) = -B(idx2, idx1);

                    // update Green's function
                    cVecType x1 = -g.col(idx1);
                    cVecType x2 = -g.col(idx2);
                    x1(idx1) += 2;
                    x2(idx2) += 2;
                    alpha = (+1.0i) * double(auxCur) * (config->thlV) / tmp[imaj];
                    zgeru(&nDim, &nDim, &alpha, x1.data(), &inc, x2.data(), &inc,
                        g.data(), &nDim);
                    alpha = -alpha;
                    zgeru(&nDim, &nDim, &alpha, x2.data(), &inc, x1.data(), &inc,
                        g.data(), &nDim);
                }
            } else if (hsScheme == 1) {
                int idxj1, idxk1, idxj2, idxk2;
                config->aux2MajoranaIdx(idxAux, 0, bondType, idxj1, idxk1);
                config->aux2MajoranaIdx(idxAux, 1, bondType, idxj2, idxk2);
                for (int iaux = 0; iaux < 2; iaux++) {
                    if (iaux == 0) {
                        idx1 = idxj1;
                        idx2 = idxj2;
                    } else {
                        idx1 = idxk1;
                        idx2 = idxk2;
                        tmp[1] = (1.0 + ((1.0i) * (config->thlV) * double(auxCur) *
                                        g(idx1, idx2)));
                    }

                    // update aux field and B matrix
                    (*s)(idxAux) = -auxCur;
                    B(idx1, idx2) = -B(idx1, idx2);
                    B(idx2, idx1) = -B(idx2, idx1);

                    // update Green's function
                    cVecType x1 = -g.col(idx1);
                    cVecType x2 = -g.col(idx2);
                    x1(idx1) += 2;
                    x2(idx2) += 2;
                    alpha = (+1.0i) * double(auxCur) * (config->thlV) / tmp[iaux];
                    if (iaux == 1) alpha = -alpha;
                    zgeru(&nDim, &nDim, &alpha, x1.data(), &inc, x2.data(), &inc,
                        g.data(), &nDim);
                    alpha = -alpha;
                    zgeru(&nDim, &nDim, &alpha, x2.data(), &inc, x1.data(), &inc,
                        g.data(), &nDim);
                }
            }
        }
    };

    void singleFlipSingleMajorana(MatType &g, int idxAux, double rand,
                                  bool &flag, DataType &signCur) {
        DataType r;
        // auto m = mConfig->idxCell2Coord(idxCell);
        int auxCur = (*s)(idxAux);
        int idx1, idx2;
        DataType tmp;
        const int inc = 1;
        DataType alpha;

        config->aux2MajoranaIdx(idxAux, 0, bondType, idx1, idx2);
        tmp =
            (1.0 - ((1.0i) * (config->thlV) * double(auxCur) * g(idx1, idx2)));
        r = tmp * tmp * etaM;

        flag = rand < std::abs(r);
        // std::cout << rand << "=rand " << "r = " << r << "\n";

        if (flag) {
            // std::cout << "tmp = " << tmp << " r = " << r << "\n";
            signCur *= (tmp / std::abs(tmp));
            config->aux2MajoranaIdx(idxAux, 0, bondType, idx1, idx2);
            // update aux field and B matrix
            (*s)(idxAux) = -auxCur;
            B(idx1, idx2) = -B(idx1, idx2);
            B(idx2, idx1) = -B(idx2, idx1);

            // update Green's function
            cVecType x1 = -g.col(idx1);
            cVecType x2 = -g.col(idx2);
            x1(idx1) += 2;
            x2(idx2) += 2;
            alpha = (+1.0i) * double(auxCur) * (config->thlV) / tmp;
            zgeru(&nDim, &nDim, &alpha, x1.data(), &inc, x2.data(), &inc,
                  g.data(), &nDim);
            alpha = -alpha;
            zgeru(&nDim, &nDim, &alpha, x2.data(), &inc, x1.data(), &inc,
                  g.data(), &nDim);
        }
    };

    DataType update(MatType &g) override {
        double rand;
        bool flag;
        DataType signCur = 1.0;
        if (singleMaj) {
            for (int i = 0; i < s->size(); i++) {
                // random real between (0, 1)
                rand = rd->rdUniform01();
                singleFlipSingleMajorana(g, i, rand, flag, signCur);
            }
        } else {
            for (int i = 0; i < s->size(); i++) {
                // random real between (0, 1)
                rand = rd->rdUniform01();
                singleFlip(g, i, rand, flag, signCur);
            }
        }
        return signCur;
    }

    void right_multiply(const MatType &AIn, MatType &AOut) override {
        AOut = AIn * B;
    }

    void left_multiply(const MatType &AIn, MatType &Aout) override {
        Aout = B * AIn;
    }

    void left_propagate(MatType &g, MatType &gTmp) override {
        reCalcInv();
        gTmp = B * g;
        g = gTmp * B_inv;
    }
    void right_propagate(MatType &g, MatType &gTmp) override {
        reCalcInv();
        gTmp = B_inv * g;
        g = gTmp * B;
    }

    iVecType *getAuxField() override { return s; }

    int getType() override { return bondType; }

    void getGreensMat(MatType &g) override {
        g = MatType::Zero(nDim, nDim);
        config->InteractionTanhGenerator(g, *s, bondType, false);
    }

    void getGreensMatInv(MatType &g) override {
        g = MatType::Zero(nDim, nDim);
        config->InteractionTanhGenerator(g, *s, bondType, true);
    }

    // inline DataType getSignPfGInv() override;

    void stabilizedLeftMultiply(UDT &F) override {
        // F.bMultUpdate(B);
        F = B * F;
    }
};

class Spinless_tV {
   public:
    std::vector<Operator *> op_array;
    int nDim;
    ~Spinless_tV() {
        for (int i = 0; i < op_array.size(); i++) {
            delete op_array[i];
        }
    }
};

#endif