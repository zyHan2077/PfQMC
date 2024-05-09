#ifndef HONEYCOMB_SINGLE_MAJORANA_H
#define HONEYCOMB_SINGLE_MAJORANA_H

#include "operator.h"
#include "spinless_tV.h"
#include "types.h"

class SpinlessTvHoneycombSingleMajoranaUtils : public SpinlessTvUtils {
   public:
    struct majoranaCoord {
        int ix, iy, iSubcell;
        majoranaCoord(int _ix = 0, int _iy = 0, int _iSubcell = 0)
            : ix(_ix), iy(_iy), iSubcell(_iSubcell){};
    };

    int nsites;
    int nUnitcell;

    SpinlessTvHoneycombSingleMajoranaUtils(int _Lx, int _Ly, double _dt,
                                           double _V, int _l)
        : SpinlessTvUtils(_Lx, _Ly, _dt, _V, _l, _Lx * _Ly * 2, true) {
        nUnitcell = Lx * Ly;
        nsites = nUnitcell * 2;
    }

    inline majoranaCoord idxCell2Coord(int idx) const {
        return majoranaCoord(idx / Ly, idx % Ly, 0);
    }

    inline int majoranaCoord2Idx(majoranaCoord s) const {
        return unitCellCoord2Idx(s.ix, s.iy)*2 + s.iSubcell;
    }

    inline int majoranaCoord2Idx(int ix, int iy, int isubcell) const {
        return unitCellCoord2Idx(ix, iy)*2 + isubcell;
    }

    inline int neighborSiteIdx(int ix, int iy, int bondType) const {
        switch (bondType) {
            case 0:
                return majoranaCoord2Idx(ix, iy, 1);
            case 1:
                return majoranaCoord2Idx((ix+1) % Lx, iy, 1);
            default:
                return majoranaCoord2Idx(ix, (iy+1) % Ly, 1);
        }
    }

    inline void aux2MajoranaIdx(int idAux, int , int bType, int& idx1, int& idx2) const override {
        majoranaCoord m = idxCell2Coord(idAux);
        idx1 = majoranaCoord2Idx(m.ix, m.iy, 0);
        idx2 = neighborSiteIdx(m.ix, m.iy, bType);
    }

    inline void KineticGenerator(MatType &H, DataType t) const {
        H.setZero();
        int idx1, idx2;
        DataType tmp = (1.0i) * t;
        for (int i=0; i<Lx; i++) {
            for (int j=0; j<Ly; j++) {
                idx1 = majoranaCoord2Idx(i, j, 0);
                idx2 = majoranaCoord2Idx(i, j, 1);
                H(idx1, idx2) = tmp;
                H(idx2, idx1) = -tmp;
                idx2 = majoranaCoord2Idx((i+1) % Lx, j, 1);
                H(idx1, idx2) = tmp;
                H(idx2, idx1) = -tmp;
                idx2 = majoranaCoord2Idx(i, (j+1) % Ly, 1);
                H(idx1, idx2) = tmp;
                H(idx2, idx1) = -tmp;
            }
        }
    }

    inline DataType energyFromGreensFunc(const MatType &g) {
        DataType r = 0.0;
        DataType tmp = (0.5i) * 2.0;
        int idx1, idx2;
        for (int i=0; i<Lx; i++) {
            for (int j=0; j<Ly; j++) {
                idx1 = majoranaCoord2Idx(i, j, 0);
                idx2 = majoranaCoord2Idx(i, j, 1);
                r += tmp * g(idx1, idx2);
                idx2 = majoranaCoord2Idx((i+1) % Lx, j, 1);
                r += tmp * g(idx1, idx2);
                idx2 = majoranaCoord2Idx(i, (j+1) % Ly, 1);
                r += tmp * g(idx1, idx2);
            }
        }

        int idxi, idxj;
        tmp = (0.25) * V;
        for (int i=0; i<Lx; i++) {
            for (int j=0; j<Ly; j++) {
                for (int btype=0; btype<3; btype++) {
                    idxi = majoranaCoord2Idx(i, j, 0); // i1
                    idxj = neighborSiteIdx(i, j, btype); // j1
                    r += tmp * g(idxi, idxj) * g(idxi, idxj);
                }
            }
        }
        return r;
    }

    // The structure factor is defined by:
    // S(k) = 1/N^2 \sum_{i,j} s_{ij} <(n_i - 1/2)(n_j - 1/2)>
    //      = 1/(4 N^2) \sum_{i,j} s_{ij} <\gamma^1_i \gamma^1_j> <\gamma^2_i \gamma^2_j>
    // where s_{ij} is :
    // +1 if i and j are in the same unit cell
    // -1 if i and j are in different unit cells
    inline DataType structureFactorCDW(const MatType &g) {
        DataType r = 0.0;
        for (int i=0; i<nsites; i++) {
            for (int j=0; j<nsites; j++) {
                if ( (i+j) % 2 == 0 ) {
                    r += g(i, j) * g(i, j);
                } else {
                    r -= g(i, j) * g(i, j);
                }
            }
        }
        return r / (4.0 * nsites * nsites);
    }
};

class HoneycombSingleMajorana_tV : public Spinless_tV {
public:
    const SpinlessTvHoneycombSingleMajoranaUtils* modelConfig;
    double dt;
    int l;
    int nSites, nUnitcell;
    rdGenerator* rd;

    HoneycombSingleMajorana_tV(SpinlessTvHoneycombSingleMajoranaUtils* _modelConfig,  rdGenerator* _rd) {
        modelConfig = _modelConfig;
        dt = _modelConfig->dt;
        l = _modelConfig->l;
        nUnitcell = _modelConfig->nUnitcell;
        nSites = _modelConfig->nsites;
        rd = _rd;

        nDim = nSites;
        MatType Ht(nDim, nDim);
        const MatType identity = MatType::Identity(nDim, nDim);
        modelConfig->KineticGenerator(Ht, 1.0);
        MatType expK = expm(Ht, -dt);
        MatType expKhalf = expm(Ht, -0.5*dt);
        op_array = std::vector<Operator*>(4*l + 1);
        iVecType* s;
        for (int i=0; i<l; i++) {
            if (i == 0) {
                op_array[0] = new DenseOperator(expKhalf, 1.0);
            } else {
                op_array[4*i] = new DenseOperator(expK, 1.0);
            }

            for (int j=0; j<3; j++) {
                s = new iVecType(nUnitcell);
                for (int k=0; k<nUnitcell; k++) (*s)(k) = rd->rdZ2();
                op_array[4*i + j + 1] = new SpinlessVOperator(modelConfig, s, j, rd);
            }
        }
        op_array[4*l] = new DenseOperator(expKhalf, 1.0);
    }
};

#endif