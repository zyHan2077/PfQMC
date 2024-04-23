#ifndef SquareLattice_H
#define SquareLattice_H

#include "operator.h"
#include "skewMatUtils.h"
#include "spinless_tV.h"
#include "types.h"

// a square lattice is defined by
// - Lx * Ly unit cells
// - 1 sites per unit cell
// - 2 bonds per unit cell
// - bonds are divided into 4 mutually commuting group
// |---1---|---2---|---1---|---2---|
// 3       3       3       3       3
// |---1---|---2---|---1---|---2---|
// 4       4       4       4       4
// |---1---|---2---|---1---|---2---|
// so
// - 2 majorana species
class SpinlessTvSquareUtils : public SpinlessTvUtils {
   public:
    int nsites;

    SpinlessTvSquareUtils(int _Lx, int _Ly, double _dt, double _V, int _l)
        : SpinlessTvUtils(_Lx, _Ly, _dt, _V, _l, _Lx * _Ly * 2) {
        // model configuration
        nsites = Lx * Ly;
    }

    inline int majoranaCoord2Idx(int ix, int iy, int imaj) const {
        return nsites * imaj + unitCellCoord2Idx(ix, iy);
    }

    inline void KineticGenerator(MatType &H, DataType t) const {
        H.setZero();
        int idx1, idx2;
        DataType tmp = (1.0i) * t;
        for (int i = 0; i < Lx; i++) {
            for (int j = 0; j < Ly; j++) {
                for (int k = 0; k < 2; k++) {
                    idx1 = majoranaCoord2Idx(i, j, k);
                    idx2 = majoranaCoord2Idx((i + 1) % Lx, j, k);
                    H(idx1, idx2) = tmp;
                    H(idx2, idx1) = -tmp;
                    idx2 = majoranaCoord2Idx(i, (j + 1) % Ly, k);
                    H(idx1, idx2) = tmp;
                    H(idx2, idx1) = -tmp;
                }
            }
        }
    }

    inline DataType energyFromGreensFunc(const MatType &g) {
        DataType r = 0.0;
        DataType tmp = (0.5i);
        int idx1, idx2;
        for (int i = 0; i < Lx; i++) {
            for (int j = 0; j < Ly; j++) {
                for (int k = 0; k < 2; k++) {
                    idx1 = majoranaCoord2Idx(i, j, k);
                    idx2 = majoranaCoord2Idx((i + 1) % Lx, j, k);
                    r += tmp * g(idx1, idx2);
                    idx2 = majoranaCoord2Idx(i, (j + 1) % Ly, k);
                    r += tmp * g(idx1, idx2);
                }
            }
        }

        int idxi1, idxi2, idxj1, idxj2;
        tmp = (0.25) * V;
        for (int i = 0; i < Lx; i++) {
            for (int j = 0; j < Ly; j++) {
                idxi1 = majoranaCoord2Idx(i, j, 0);             // i1
                idxi2 = majoranaCoord2Idx(i, j, 1);             // i2
                idxj1 = majoranaCoord2Idx((i + 1) % Lx, j, 0);  // j1
                idxj2 = majoranaCoord2Idx((i + 1) % Lx, j, 1);  // j2
                r += tmp * g(idxi1, idxj1) * g(idxi2, idxj2);
                r += tmp * g(idxi1, idxj2) * g(idxj1, idxi2);
                r -= tmp * g(idxi1, idxi2) * g(idxj1, idxj2);

                idxj1 = majoranaCoord2Idx(i, (j + 1) % Ly, 0);  // j1
                idxj2 = majoranaCoord2Idx(i, (j + 1) % Ly, 1);  // j2
                r += tmp * g(idxi1, idxj1) * g(idxi2, idxj2);
                r += tmp * g(idxi1, idxj2) * g(idxj1, idxi2);
                r -= tmp * g(idxi1, idxi2) * g(idxj1, idxj2);
            }
        }
        return r;
    }

    // Directly generate B by directly writing each 2*2 block
    // B should be initialized as Identity
    inline void InteractionBGenerator(MatType &B, const iVecType &s,
                                      const int bondType,
                                      bool inv = false) const override {
        DataType ch = chlV;
        DataType ish = (1.0i) * shlV;
        if (inv) ish = -ish;
        int idx1, idx2;
        int auxCount = 0;

        if ((bondType == 0) || (bondType == 1)) {
            // std::cout << "here!" << "\n";
            for (int i = bondType; i < Lx; i += 2) {
                for (int j = 0; j < Ly; j++) {
                    for (int k = 0; k < 2; k++) {
                        idx1 = majoranaCoord2Idx(i, j, k);
                        idx2 = majoranaCoord2Idx((i + 1) % Lx, j, k);
                        //  |   \cosh(\lambda)           ,  -i \sinh(\lambda)
                        //  \sigma | |  +i \sinh(\lambda) \sigma  ,
                        //  \cosh(\lambda)          |
                        B(idx1, idx1) = ch;
                        B(idx2, idx2) = ch;
                        // std::cout << idx1 << " " << idx2 << " " << ish <<
                        // "\n";
                        B(idx1, idx2) = +ish * double(s(auxCount));
                        B(idx2, idx1) = -ish * double(s(auxCount));
                    }
                    auxCount++;
                }
            }
        } else {
            for (int i = 0; i < Lx; i++) {
                for (int j = bondType - 2; j < Ly; j += 2) {
                    for (int k = 0; k < 2; k++) {
                        idx1 = majoranaCoord2Idx(i, j, k);
                        idx2 = majoranaCoord2Idx(i, (j + 1) % Ly, k);
                        B(idx1, idx1) = ch;
                        B(idx2, idx2) = ch;
                        B(idx1, idx2) = +ish * double(s(auxCount));
                        B(idx2, idx1) = -ish * double(s(auxCount));
                    }
                    auxCount++;
                }
            }
        }
    }
};

class SpinlessVSquareOperator : public SpinlessVOperator {
   public:
    const SpinlessTvSquareUtils *mConfig;

    SpinlessVSquareOperator(const SpinlessTvSquareUtils *_config, iVecType *_s,
                            int _bondType, rdGenerator *_rd)
        : SpinlessVOperator(_config, _s, _bondType, _rd), mConfig(_config) {}

    void aux2MajoranaIdx(int idAux, int imaj, int &idx1, int &idx2) override {
        int ix, iy;
        int Lx = mConfig->Lx;
        int Ly = mConfig->Ly;
        if (bondType == 0 || bondType == 1) {
            iy = idAux % Ly;
            ix = (idAux / Ly) * 2 + bondType;
            idx1 = mConfig->majoranaCoord2Idx(ix, iy, imaj);
            idx2 = mConfig->majoranaCoord2Idx((ix + 1) % Lx, iy, imaj);
        } else {
            ix = idAux % Lx;
            iy = (idAux / Lx) * 2 + bondType - 2;
            idx1 = mConfig->majoranaCoord2Idx(ix, iy, imaj);
            idx2 = mConfig->majoranaCoord2Idx(ix, (iy + 1) / Ly, imaj);
        }
    }
};

class Square_tV : public Spinless_tV {
   public:
    const SpinlessTvSquareUtils *modelConfig;
    double dt;
    int l;
    int nSites;
    rdGenerator *rd;
    int nBond[4];

    Square_tV(SpinlessTvSquareUtils *_config, rdGenerator *_rd) {
        modelConfig = _config;
        dt = modelConfig->dt;
        l = modelConfig->l;
        nSites = modelConfig->nsites;
        rd = _rd;

        nDim = nSites * 2;
        MatType Ht(nDim, nDim);
        const MatType identity = MatType::Identity(nDim, nDim);
        modelConfig->KineticGenerator(Ht, 1.0);
        MatType expK = expm(Ht, -dt);
        MatType expKhalf = expm(Ht, -dt / 2.0);

        // number of bonds, therefore auxillary fields
        // for each species
        nBond[1] = (modelConfig->Lx / 2) * modelConfig->Ly;
        nBond[0] = nSites - nBond[1];
        nBond[3] = (modelConfig->Ly / 2) * modelConfig->Lx;
        nBond[2] = nSites - nBond[3];

        // for (int i = 0; i<4; i++)
        // std::cout << nBond[i] << "\n";

        op_array = std::vector<Operator *>(5 * l + 1);
        iVecType *s;
        for (int i = 0; i < l; i++) {
            if (i == 0) {
                op_array[0] = new DenseOperator(expKhalf, 1.0);
            } else {
                op_array[5 * i] = new DenseOperator(expK, 1.0);
            }

            for (int j = 0; j < 4; j++) {
                s = new iVecType(nBond[j]);
                for (int k = 0; k < nBond[j]; k++) (*s)(k) = rd->rdZ2();
                op_array[5 * i + j + 1] =
                    new SpinlessVSquareOperator(modelConfig, s, j, rd);
            }
        }
        op_array[5 * l] = new DenseOperator(expKhalf, 1.0);
    }

    ~Square_tV() {
        for (int i = 0; i < 5 * l + 1; i++) {
            delete op_array[i];
        }
    }
};

#endif