#ifndef KITAEVCHAIN_H
#define KITAEVCHAIN_H

#include "operator.h"
#include "skewMatUtils.h"
#include "spinless_tV.h"
#include "types.h"
#include <iostream>

class  SpinlessTvChainUtils : public SpinlessTvUtils {
  public:
    int nsites;
    double delta;
    double mu;
    int boundaryType; // 0: PBC, 1: OBC

    int hsScheme;

    SpinlessTvChainUtils(int _L, double _dt, double _V, int _l, int _boundary, double _delta=0.0, double _mu=0.0, int _hsScheme=0) 
        : SpinlessTvUtils(_L, 1, _dt, _V, _l, _L*2, false, _hsScheme) {
        
        boundaryType = _boundary;
        nsites = _L;
        delta = _delta;
        mu = _mu;
        hsScheme = _hsScheme;
    }

    inline int majoranaCoord2Idx(int ix, int imaj) const {
        return nsites * imaj + ix;
    }

    // for 1d, btype = 0 or 1
    inline void aux2MajoranaIdx(int idAux, int imaj, int bType, int &idx1, int &idx2) const override {
        int ix;
        ix = idAux * 2 + bType;
        if (bType == 0) {
            idx1 = majoranaCoord2Idx(ix, imaj);
            idx2 = majoranaCoord2Idx((ix+1)%Lx, imaj);
        } else {
            idx2 = majoranaCoord2Idx(ix, imaj);
            idx1 = majoranaCoord2Idx((ix+1)%Lx, imaj);
        }
    }

    // boundaryType = 0 : PBC, boundaryType = 1 : OBC
    inline void KineticGenerator(MatType &H) const {
        H.setZero();
        int idx1, idx2;
        DataType tmp = (1.0i);
        DataType tmpDelta = (1.0i)*delta;
        DataType tmpMu = (1.0i)*mu;
        // int Li = Lx;
        // if (boundaryType == 1) {
        //     Li = Lx - 1;
        // }
        for (int i=0; i<Lx-1; i++) {
            for (int k=0; k<2; k++) {
                idx1 = majoranaCoord2Idx(i, k);
                idx2 = majoranaCoord2Idx((i+1)%Lx, k);

                if((i % 2) == 0) {
                    H(idx1, idx2) = tmp;
                    H(idx2, idx1) = -tmp;
                } else {
                    H(idx1, idx2) = -tmp;
                    H(idx2, idx1) = +tmp;
                }

            }

            // pairing part
            idx1 = majoranaCoord2Idx(i, 0);
            idx2 = majoranaCoord2Idx((i+1)%Lx, 0);
            H(idx1, idx2) += tmpDelta;
            H(idx2, idx1) += -tmpDelta;

            idx1 = majoranaCoord2Idx(i, 1);
            idx2 = majoranaCoord2Idx((i+1)%Lx, 1);
            H(idx1, idx2) += -tmpDelta;
            H(idx2, idx1) += tmpDelta;
        }

        if (boundaryType == 0) { // PBC
            int idxL1 = majoranaCoord2Idx(Lx-1, 0);
            int idxL2 = majoranaCoord2Idx(Lx-1, 1);
            int idx11 = majoranaCoord2Idx(0, 0);
            int idx12 = majoranaCoord2Idx(0, 1);
            if (Lx % 2 == 0) {
                // idx1 = majoranaCoord2Idx(Lx-1, 0);
                // idx2 = majoranaCoord2Idx(0, 0);
                H(idxL1, idx11) += tmpDelta;
                H(idx11, idxL1) += -tmpDelta;
                H(idxL1, idx11) += -tmp;
                H(idx11, idxL1) += +tmp;

                // idx1 = majoranaCoord2Idx(Lx-1, 1);
                // idx2 = majoranaCoord2Idx(0, 1);
                H(idxL2, idx12) += -tmpDelta;
                H(idx12, idxL2) += tmpDelta;
                H(idxL2, idx12) += -tmp;
                H(idx12, idxL2) += +tmp;
            } else {
                H(idxL1, idx12) += -tmp;
                H(idx12, idxL1) += +tmp;
                H(idxL2, idx11) += +tmp;
                H(idx11, idxL2) += -tmp;

                H(idxL1, idx12) += -tmpDelta;
                H(idx12, idxL1) += +tmpDelta;
                H(idxL2, idx11) += -tmpDelta;
                H(idx11, idxL2) += +tmpDelta;    
            }
        }

        // chemical potential
        for (int i=0; i<Lx; i++) {
            idx1 = majoranaCoord2Idx(i, 0);
            idx2 = majoranaCoord2Idx(i, 1);
            H(idx1, idx2) += -tmpMu;
            H(idx2, idx1) += +tmpMu;
        }
    }

    // this defines the correlation function between 
    // the two majorana "on the edge"
    // in principle only be used when Lx = 2n and Delta > 0
    // <i \gamma^2_1 \gamma^2_L>
    inline DataType EdgeCorrelator(const MatType &g) {
        int idx1, idx2;
        idx1 = majoranaCoord2Idx(0, 1);
        idx2 = majoranaCoord2Idx(Lx-1, 1);
        return (1.0i) * g(idx1, idx2);
    }

    inline DataType StructureFactorCDW(const MatType &g) const {
        DataType r = 0.0;
        int idxi1, idxi2, idxj1, idxj2;
        DataType tmp;
        for (int i = 0; i < Lx; i++) {
            for (int j = 0; j < Lx; j++) {
                idxi1 = majoranaCoord2Idx(i, 0);
                idxi2 = majoranaCoord2Idx(i, 1);
                idxj1 = majoranaCoord2Idx(j, 0);
                idxj2 = majoranaCoord2Idx(j, 1);
                tmp = (-g(idxi1, idxj1) * g(idxi2, idxj2) + g(idxi1, idxi2) * g(idxj1, idxj2) + g(idxi1, idxj2) * g(idxi2, idxj1));
                if ( (i+j) % 2 == 0) {
                    r += tmp;
                } else {
                    r -= tmp;
                }
            }
        }

        return r / (4.0 * Lx * Lx);
    }

    inline DataType StructureFactorCDWOffset(const MatType &g) const {
        DataType r = 0.0;
        int idxi1, idxi2, idxj1, idxj2;
        DataType tmp;
        DataType deltaq = 2.0 * M_PI / double(Lx);
        for (int i = 0; i < Lx; i++) {
            for (int j = 0; j < Lx; j++) {
                idxi1 = majoranaCoord2Idx(i, 0);
                idxi2 = majoranaCoord2Idx(i, 1);
                idxj1 = majoranaCoord2Idx(j, 0);
                idxj2 = majoranaCoord2Idx(j, 1);
                tmp = (-g(idxi1, idxj1) * g(idxi2, idxj2) + g(idxi1, idxi2) * g(idxj1, idxj2) + g(idxi1, idxj2) * g(idxi2, idxj1));
                if ( (i+j) % 2 == 0) {
                    r += tmp * exp(1.0i * deltaq * double(i - j) );
                } else {
                    r -= tmp * exp(1.0i * deltaq * double(i - j) );
                }
            }
        }

        return r / (4.0 * Lx * Lx);
    }

    inline void StructureFactorCDWFFT(const MatType &g, cVecType &out) const {
        int idxi1, idxi2, idxj1, idxj2;
        // double r;
        DataType deltaq = 2.0 * M_PI / double(Lx);
        DataType curQ;
        cVecType r = cVecType::Zero(Lx);
        DataType tmp;
        for (int i = 0; i < Lx; i++) {
            for (int j = 0; j < Lx; j++) {
                idxi1 = majoranaCoord2Idx(i, 0);
                idxi2 = majoranaCoord2Idx(i, 1);
                idxj1 = majoranaCoord2Idx(j, 0);
                idxj2 = majoranaCoord2Idx(j, 1);
                tmp = (-g(idxi1, idxj1) * g(idxi2, idxj2) + g(idxi1, idxi2) * g(idxj1, idxj2) + g(idxi1, idxj2) * g(idxi2, idxj1));
                for (int k=0; k<Lx; k++) {
                    curQ = deltaq * double(k);
                    r(k) += tmp * exp(1.0i * curQ * double(i - j) );
                }
            }
        }
        out = r / (4.0 * Lx * Lx);
    }

    // inline DataType StructureFactorCDWM4(const MatType &g) const {
    //     DataType r = 0.0;
    //     int idxi1, idxi2, idxj1, idxj2, idxk1, idxk2, idxl1, idxl2;
    //     DataType tmp;
    //     for (int i = 0; i < Lx; i++) {
    //         for (int j = 0; j < Lx; j++) {
    //             for (int k = 0; k < Lx; k++) {
    //                 for (int l = 0; l < Lx; l++) {
    //                     idxi1 = majoranaCoord2Idx(i, 0);
    //                     idxi2 = majoranaCoord2Idx(i, 1);
    //                     idxj1 = majoranaCoord2Idx(j, 0);
    //                     idxj2 = majoranaCoord2Idx(j, 1);
    //                     idxk1 = majoranaCoord2Idx(k, 0);
    //                     idxk2 = majoranaCoord2Idx(k, 1);
    //                     idxl1 = majoranaCoord2Idx(l, 0);
    //                     idxl2 = majoranaCoord2Idx(l, 1);
                        
    //                     tmp = (g(idxi1, idxj1) * g(idxk1, idxl1) +  g(idxi1, idxl1) * g(idxj1, idxk1) - g(idxi1, idxk1) * g(idxj1, idxl1));
                        
    //                     tmp *= (g(idxi2, idxj2) * g(idxk2, idxl2) +  g(idxi2, idxl2) * g(idxj2, idxk2) - g(idxi2, idxk2) * g(idxj2, idxl2));

    //                     if ( (i+j+k+l) % 2 == 0) {
    //                         r += tmp;
    //                     } else {
    //                         r -= tmp;
    //                     }
    //                 }
    //             }
    //         }
    //     }

    //     return r / (16.0 * Lx * Lx * Lx * Lx);
    // }

    inline DataType Z2FermionParity(const MatType &g) const {
        DataType additionalSign;
        switch (Lx % 4) {
            case 0:
                additionalSign = 1.0;
                break;
            case 1:
                additionalSign = -(1.0i);
                break;
            case 2:
                additionalSign = -1.0;
                break;
            case 3:
                additionalSign = (1.0i);
                break;
        }

        int tmp = Lx * (Lx -1) /2;
        if (tmp % 2 == 1) {
            additionalSign *= -1.0;
        }
        MatType gcopy = g;
        DataType r = additionalSign * pfaf(Lx, gcopy);
        // std::cout << "Z2FermionParity = " << r << std::endl;
        return r;
    }

    inline DataType Z2FermionParityEdgeCorrelator(const MatType &g) const {
        DataType additionalSign;
        switch (Lx % 4) {
            case 0:
                additionalSign = 1.0;
                break;
            case 1:
                additionalSign = -(1.0i);
                break;
            case 2:
                additionalSign = -1.0;
                break;
            case 3:
                additionalSign = (1.0i);
                break;
        }

        int idx1, idx2;
        idx1 = majoranaCoord2Idx(0, 1);
        idx2 = majoranaCoord2Idx(Lx-1, 1);
        iVecType curInd(2*Lx - 2);
        int count = 0;
        for (int i = 0; i < 2*Lx; i++) {
            if (i != idx1 && i != idx2) {
                curInd(count) = i;
                count++;
            }
        }

        // sign = (-i)^(Lx) * (-1)^(Lx(Lx-1)/2) * (-1)^(Lx-1) * i
        int tmp = Lx * (Lx -1) /2;
        if (tmp % 2 == 1) {
            additionalSign *= -1.0;
        }

        if(Lx % 2 == 0) {
            additionalSign *= -1.0;
        }

        MatType gcopy = g(curInd, curInd);

        return (1.0i) * additionalSign * pfaf(Lx-1, gcopy);
    }

    inline DataType energyFromGreensFunc(const MatType & g) {
        DataType r = 0.0;
        DataType tmp = (0.5i);
        DataType tmpDelta = (0.5i) * delta;
        DataType tmpMu = (0.5i) * mu;
        int idx1, idx2;
        int Li = Lx;
        
        // if (boundaryType == 1) {
        //     Li = Lx - 1;
        // }

        for (int i=0; i<Lx-1; i++) {
            for (int k=0; k<2; k++) {
                idx1 = majoranaCoord2Idx(i, k);
                idx2 = majoranaCoord2Idx((i + 1) % Lx, k);
                if ((i % 2) == 0) {
                    r += tmp * g(idx1, idx2);
                } else {
                    r -= tmp * g(idx1, idx2);
                }
            }

             // pairing part
            idx1 = majoranaCoord2Idx(i, 0);
            idx2 = majoranaCoord2Idx((i + 1) % Lx, 0);
            r += tmpDelta * g(idx1, idx2);

            idx1 = majoranaCoord2Idx(i, 1);
            idx2 = majoranaCoord2Idx((i + 1) % Lx, 1);
            r -= tmpDelta * g(idx1, idx2);

        }

        if (boundaryType == 0) { // PBC
            int idxL1 = majoranaCoord2Idx(Lx-1, 0);
            int idxL2 = majoranaCoord2Idx(Lx-1, 1);
            int idx11 = majoranaCoord2Idx(0, 0);
            int idx12 = majoranaCoord2Idx(0, 1);
            if (Lx % 2 == 0) {
                r += tmpDelta * g(idxL1, idx11);
                r += -tmp * g(idxL1, idx11);

                r += -tmpDelta * g(idxL2, idx12);
                r += -tmp * g(idxL2, idx12);
                // H(idxL1, idx11) += tmpDelta;
                // H(idx11, idxL1) += -tmpDelta;
                // H(idxL1, idx11) += -tmp;
                // H(idx11, idxL1) += +tmp;

                // H(idxL2, idx12) += -tmpDelta;
                // H(idx12, idxL2) += tmpDelta;
                // H(idxL2, idx12) += -tmp;
                // H(idx12, idxL2) += +tmp;
            } else {
                r += -tmp * g(idxL1, idx12);
                r += tmp * g(idxL2, idx11);
                // H(idxL1, idx12) += -tmp;
                // H(idx12, idxL1) += +tmp;
                // H(idxL2, idx11) += +tmp;
                // H(idx11, idxL2) += -tmp;
                r += -tmpDelta * g(idxL1, idx12);
                r += -tmpDelta * g(idxL2, idx11);
                // H(idxL1, idx12) += -tmpDelta;
                // H(idx12, idxL1) += +tmpDelta;
                // H(idxL2, idx11) += -tmpDelta;
                // H(idx11, idxL2) += +tmpDelta;    
            }
        }

        int idxi1, idxi2, idxj1, idxj2;
        tmp = (0.25) * V;
        for (int i = 0; i < Lx-1; i++) {
            idxi1 = majoranaCoord2Idx(i, 0);             // i1
            idxi2 = majoranaCoord2Idx(i, 1);             // i2
            idxj1 = majoranaCoord2Idx((i + 1) % Lx, 0);  // j1
            idxj2 = majoranaCoord2Idx((i + 1) % Lx, 1);  // j2
            r += tmp * g(idxi1, idxj1) * g(idxi2, idxj2);
            r += tmp * g(idxi1, idxj2) * g(idxj1, idxi2);
            r -= tmp * g(idxi1, idxi2) * g(idxj1, idxj2);
        }

        //chemical potential part
        for (int i=0; i<Lx; i++) {
            idx1 = majoranaCoord2Idx(i, 0);
            idx2 = majoranaCoord2Idx(i, 1);
            r -= tmpMu * g(idx1, idx2);
        }

       return r;

    }

};

class Chain_tV : public Spinless_tV {
  public:
    const SpinlessTvChainUtils *modelConfig;
    double dt; 
    int l, nSites;
    rdGenerator *rd;
    int nBond[2];

    Chain_tV(SpinlessTvChainUtils *_config, rdGenerator *_rd) {
        modelConfig = _config;
        dt = _config->dt;
        l = _config->l;
        nSites = _config->nsites;
        rd = _rd;

        nDim = 2 * nSites;

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

        if (modelConfig->boundaryType == 0) { // PBC
            if (nSites % 2 == 0) {
                nBond[0] = nSites / 2;
                nBond[1] = nSites / 2;
            } else {
                nBond[0] = (nSites + 1) / 2;
                nBond[1] = (nSites - 1) / 2;
            }
        } else { // OBC
            if (nSites % 2 == 0) {
                nBond[0] = nSites / 2;
                nBond[1] = (nSites / 2) - 1;
            } else {
                nBond[0] = (nSites - 1) / 2;
                nBond[1] = (nSites - 1) / 2;
            }
        }

        // only OBC is used for interaction terms
        // if (nSites % 2 == 0) {
        //     nBond[0] = nSites / 2;
        //     nBond[1] = (nSites / 2) - 1;
        // } else {
        //     nBond[0] = (nSites - 1) / 2;
        //     nBond[1] = (nSites - 1) / 2;
        // }

        // std::cout << "nBond = " << nBond[0] << " " << nBond[1] << std::endl;

        op_array = std::vector<Operator *>(3 * l + 1);
        iVecType *s;
        for (int i = 0; i < l; i++) {
            if (i == 0) {
                op_array[0] = new DenseOperator(expKhalf, signKHalf);
            } else {
                op_array[3 * i] = new DenseOperator(expK, signK);
            }

            for (int j = 0; j < 2; j++) {
                s = new iVecType(nBond[j]);
                for (int k = 0; k < nBond[j]; k++) (*s)(k) = rd->rdZ2();
                op_array[3 * i + j + 1] =
                    new SpinlessVOperator(modelConfig, s, j, rd);
            }
        }
        op_array[3 * l] = new DenseOperator(expKhalf, signKHalf);
    }
};

#endif
