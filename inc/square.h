#ifndef SquareLattice_H
#define SquareLattice_H

#include "operator.h"
#include "skewMatUtils.h"
#include "spinless_tV.h"
#include "types.h"
/* 
 * This macro is purely for fun.
 * For 4n*4m lattice, rearanging the
 * hopping signs introduces no changes
 * to the final result. However for e.g 4*6
 * lattice, this macro defines a model which
 * is different from the original t-V model
 */
// #define INVERSEBOND

// a square lattice is defined by
// - Lx * Ly unit cells
// - 1 sites per unit cell
// - 2 bonds per unit cell
// - bonds are divided into 4 mutually commuting group
// |---0---|---1---|---0---|---1---|
// 3       3       3       3       3
// |---0---|---1---|---0---|---1---|
// 2       2       2       2       2
// |---0---|---1---|---0---|---1---|
// so
// - 2 majorana species
class SpinlessTvSquareUtils : public SpinlessTvUtils {
   public:
    int nsites;
    double delta;

    SpinlessTvSquareUtils(int _Lx, int _Ly, double _dt, double _V, int _l, double _delta=0.0)
        : SpinlessTvUtils(_Lx, _Ly, _dt, _V, _l, _Lx * _Ly * 2) {
        // model configuration
        nsites = Lx * Ly;
        delta = _delta;
    }

    inline int majoranaCoord2Idx(int ix, int iy, int imaj) const {
        return nsites * imaj + unitCellCoord2Idx(ix, iy);
    }

    inline void aux2MajoranaIdx(int idAux, int imaj, int bType, int &idx1, int &idx2) const override {
        int ix, iy;
        if (bType == 0 || bType == 1) {
            iy = idAux % Ly;
            ix = (idAux / Ly) * 2 + bType;
            idx1 = majoranaCoord2Idx(ix, iy, imaj);
            idx2 = majoranaCoord2Idx((ix + 1) % Lx, iy, imaj);
        } else {
            ix = idAux % Lx;
            iy = (idAux / Lx) * 2 + bType - 2;
            idx1 = majoranaCoord2Idx(ix, iy, imaj);
            idx2 = majoranaCoord2Idx(ix, (iy + 1) % Ly, imaj);
        }

        #ifndef INVERSEBOND
        if ( (ix + iy) % 2 == 1) {
            std::swap(idx1, idx2);
        }   
        #endif
        // std::cout << "idAux=" << idAux << " btype=" << bondType << " " << idx1 << idx2 << "\n";
    }

    inline void KineticGenerator(MatType &H, DataType t=1.0) const {
        H.setZero();
        int idx1, idx2, idx2p;
        DataType tmp = (1.0i) * t;
        DataType tmpDelta = (1.0i) * delta;
        for (int i = 0; i < Lx; i++) {
            for (int j = 0; j < Ly; j++) {
                for (int k = 0; k < 2; k++) {
                    idx1 = majoranaCoord2Idx(i, j, k);
                    idx2 = majoranaCoord2Idx((i + 1) % Lx, j, k);
                    idx2p = majoranaCoord2Idx(i, (j + 1) % Ly, k);

                    // hopping part
                    #ifndef INVERSEBOND
                        if((i+j)% 2 == 0) {
                            H(idx1, idx2) = tmp;
                            H(idx2, idx1) = -tmp;
                            H(idx1, idx2p) = tmp;
                            H(idx2p, idx1) = -tmp;
                        } else {
                            H(idx2, idx1) = tmp;
                            H(idx1, idx2) = -tmp;
                            H(idx2p, idx1) = tmp;
                            H(idx1, idx2p) = -tmp;
                        }
                    #else
                        H(idx1, idx2) = tmp;
                        H(idx2, idx1) = -tmp;
                        H(idx1, idx2p) = tmp;
                        H(idx2p, idx1) = -tmp;
                    #endif
                }
                // pairing part, i \to i+x
                idx1 = majoranaCoord2Idx(i, j, 0);
                idx2 = majoranaCoord2Idx((i + 1) % Lx, j, 0);
                H(idx1, idx2) += tmpDelta;
                H(idx2, idx1) += -tmpDelta;
                idx1 = majoranaCoord2Idx(i, j, 1);
                idx2 = majoranaCoord2Idx((i + 1) % Lx, j, 1);
                H(idx1, idx2) += -tmpDelta;
                H(idx2, idx1) += tmpDelta;

                // pairing part, i \to i+y
                idx1 = majoranaCoord2Idx(i, j, 0);
                idx2 = majoranaCoord2Idx(i, (j + 1) % Ly, 1);
                H(idx1, idx2) += tmpDelta;
                H(idx2, idx1) += -tmpDelta;
                idx1 = majoranaCoord2Idx(i, j, 1);
                idx2 = majoranaCoord2Idx(i, (j + 1) % Ly, 0);
                H(idx1, idx2) += tmpDelta;
                H(idx2, idx1) += -tmpDelta;
            }
        }
    }

    inline DataType energyFromGreensFunc(const MatType &g) {
        DataType r = 0.0;
        DataType tmp = (0.5i);
        DataType tmpDelta = (0.5i) * delta;

        int idx1, idx2, idx2p;
        DataType r0;
        for (int i = 0; i < Lx; i++) {
            for (int j = 0; j < Ly; j++) {
                for (int k = 0; k < 2; k++) {
                    r0 = 0.0;
                    idx1 = majoranaCoord2Idx(i, j, k);
                    idx2 = majoranaCoord2Idx((i + 1) % Lx, j, k);
                    idx2p = majoranaCoord2Idx(i, (j + 1) % Ly, k);
                    r0 += tmp * g(idx1, idx2);
                    r0 += tmp * g(idx1, idx2p);
                    #ifndef INVERSEBOND
                        if((i+j)% 2 == 1) {
                            r0 = -r0;
                        }
                    #endif
                    r += r0;
                }

                // pairing part
                idx1 = majoranaCoord2Idx(i, j, 0);
                idx2 = majoranaCoord2Idx((i + 1) % Lx, j, 0);
                r += tmpDelta * g(idx1, idx2);
                idx2 = majoranaCoord2Idx(i, (j+1)%Ly, 1);
                r += tmpDelta * g(idx1, idx2);
                idx1 = majoranaCoord2Idx(i, j, 1);
                idx2 = majoranaCoord2Idx((i + 1) % Lx, j, 1);
                r -= tmpDelta * g(idx1, idx2);
                idx2 = majoranaCoord2Idx(i, (j+1)%Ly, 0);
                r += tmpDelta * g(idx1, idx2);
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

        //TODO: pairing energy
        return r;
    }
};

// class SpinlessVSquareOperator : public SpinlessVOperator {
//    public:
//     const SpinlessTvSquareUtils *mConfig;

//     SpinlessVSquareOperator(const SpinlessTvSquareUtils *_config, iVecType *_s,
//                             int _bondType, rdGenerator *_rd)
//         : SpinlessVOperator(_config, _s, _bondType, _rd), mConfig(_config) {}
// };

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
        modelConfig->KineticGenerator(Ht);
        MatType Hcopy(Ht);
        MatType expK = expm(Hcopy, -dt);
        Hcopy = Ht;
        MatType expKhalf = expm(Hcopy, -dt / 2.0);

        // sign associated with K and Khalf (should be 1.0 if dt is small enough)
        Hcopy = dt * Ht;
        DataType signK = signOfHamiltonian(Hcopy);
        Hcopy = (0.5 * dt ) * Ht;
        DataType signKHalf = signOfHamiltonian(Hcopy);
        // MatType sinhKQuater = expKQuater 

        // number of bonds, therefore auxillary fields
        // for each species
        nBond[1] = (modelConfig->Lx / 2) * modelConfig->Ly;
        nBond[0] = nSites - nBond[1];
        nBond[3] = (modelConfig->Ly / 2) * modelConfig->Lx;
        nBond[2] = nSites - nBond[3];

        // for (int i = 0; i<4; i++) std::cout << "nBond " << i << " " << nBond[i] << "\n";

        op_array = std::vector<Operator *>(5 * l + 1);
        iVecType *s;
        for (int i = 0; i < l; i++) {
            if (i == 0) {
                op_array[0] = new DenseOperator(expKhalf, signKHalf);
            } else {
                op_array[5 * i] = new DenseOperator(expK, signK);
            }

            for (int j = 0; j < 4; j++) {
                s = new iVecType(nBond[j]);
                for (int k = 0; k < nBond[j]; k++) (*s)(k) = rd->rdZ2();
                op_array[5 * i + j + 1] =
                    new SpinlessVOperator(modelConfig, s, j, rd);
            }
        }
        op_array[5 * l] = new DenseOperator(expKhalf, signKHalf);
    }
};

#endif