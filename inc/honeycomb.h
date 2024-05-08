#ifndef Honeycomb_H
#define Honeycomb_H

#include "types.h"
#include "operator.h"
#include "skewMatUtils.h"
#include "spinless_tV.h"

// a honeycomb lattice is defined by 
// - Lx * Ly unit cells
// - 2 sites per unit cell
// - 3 bonds per unit cell
// - 2 majorana species
class SpinlessTvHoneycombUtils: public SpinlessTvUtils
{
public:
    struct majoranaCoord
    {
        int ix, iy, iSubcell, imaj;
        majoranaCoord(int _ix=0, int _iy=0, int _iSubcell=0, int _imaj=0)
            : ix(_ix), iy(_iy), iSubcell(_iSubcell), imaj(_imaj) {};
    };

    int nsites;
    int nUnitcell;

    SpinlessTvHoneycombUtils(int _Lx, int _Ly, double _dt, 
            double _V, int _l) 
            :SpinlessTvUtils(_Lx, _Ly, _dt, _V, _l, _Lx*_Ly*4) {
    
        // model configuration
        nUnitcell = Lx*Ly;
        nsites = nUnitcell*2;
    }

    inline majoranaCoord idxCell2Coord(int idx) const {
        return majoranaCoord(idx / Ly, idx % Ly, 0, 0);
    }
    
    inline int majoranaCoord2Idx(majoranaCoord s) const {
        return nsites*s.imaj + unitCellCoord2Idx(s.ix, s.iy)*2 + s.iSubcell;
    }

    inline int majoranaCoord2Idx(int ix, int iy, int isubcell, int imaj) const {
        return nsites*imaj + unitCellCoord2Idx(ix, iy)*2 + isubcell;
    }

    inline int neighborSiteIdx(int ix, int iy, int imaj, int bondType) const {
        switch (bondType) {
            case 0:
                return majoranaCoord2Idx(ix, iy, 1, imaj);
            case 1:
                return majoranaCoord2Idx((ix+1) % Lx, iy, 1, imaj);
            default:
                return majoranaCoord2Idx(ix, (iy+1) % Ly, 1, imaj);
        }
    }

    inline majoranaCoord idx2MajoranaCoord(int idx) const {
        majoranaCoord s;
        s.imaj = idx / nsites;
        int tmp = idx % nsites;
        s.iSubcell = tmp % 2;
        tmp = tmp / 2;
        s.iy = tmp % Ly;
        s.ix = tmp / Ly;
        return s;
    }

    inline void aux2MajoranaIdx(int idAux, int imaj, int bType, int& idx1, int& idx2) const override {
        SpinlessTvHoneycombUtils::majoranaCoord m = idxCell2Coord(idAux);
        idx1 = majoranaCoord2Idx(m.ix, m.iy, 0, imaj);
        idx2 = neighborSiteIdx(m.ix, m.iy, imaj, bType);
    }

    inline void KineticGenerator(MatType &H, DataType t) const {
        H.setZero();
        int idx1, idx2;
        DataType tmp = (1.0i) * t;
        for (int i=0; i<Lx; i++) {
            for (int j=0; j<Ly; j++) {
                for (int k=0; k<2; k++) {
                    idx1 = majoranaCoord2Idx(i, j, 0, k);
                    idx2 = majoranaCoord2Idx(i, j, 1, k);
                    H(idx1, idx2) = tmp;
                    H(idx2, idx1) = -tmp;
                    idx2 = majoranaCoord2Idx((i+1) % Lx, j, 1, k);
                    H(idx1, idx2) = tmp;
                    H(idx2, idx1) = -tmp;
                    idx2 = majoranaCoord2Idx(i, (j+1) % Ly, 1, k);
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
        for (int i=0; i<Lx; i++) {
            for (int j=0; j<Ly; j++) {
                for (int k=0; k<2; k++) {
                    idx1 = majoranaCoord2Idx(i, j, 0, k);
                    idx2 = majoranaCoord2Idx(i, j, 1, k);
                    r += tmp * g(idx1, idx2);
                    idx2 = majoranaCoord2Idx((i+1) % Lx, j, 1, k);
                    r += tmp * g(idx1, idx2);
                    idx2 = majoranaCoord2Idx(i, (j+1) % Ly, 1, k);
                    r += tmp * g(idx1, idx2);
                }
            }
        }

        int idxi1, idxi2, idxj1, idxj2;
        tmp = (0.25) * V;
        for (int i=0; i<Lx; i++) {
            for (int j=0; j<Ly; j++) {
                for (int btype=0; btype<3; btype++) {
                    idxi1 = majoranaCoord2Idx(i, j, 0, 0); // i1
                    idxi2 = majoranaCoord2Idx(i, j, 0, 1); // i2
                    idxj1 = neighborSiteIdx(i, j, 0, btype); // j1
                    idxj2 = neighborSiteIdx(i, j, 1, btype); // j2
                    r += tmp * g(idxi1, idxj1) * g(idxi2, idxj2);
                    r += tmp * g(idxi1, idxj2) * g(idxj1, idxi2);
                    r -= tmp * g(idxi1, idxi2) * g(idxj1, idxj2);
                }
            }
        }
        return r;
    }

    // Generate the H-S transformed Hamiltonian
    // in principle ONLY be used for testing
    inline void InteractionHGenerator(MatType &H, const iVecType &s, const int bondType) const {
        DataType tmp = (1.0i) * lambdaV;

        int idx1, idx2, idUnitcell;
        for (int i=0; i<Lx; i++) {
            for (int j=0; j<Ly; j++) {
                idUnitcell = unitCellCoord2Idx(i, j);
                for (int k=0; k<2; k++) {
                    idx1 = majoranaCoord2Idx(i, j, 0, k);
                    idx2 = neighborSiteIdx(i, j, k, bondType);
                    // std::cout << idx1 << " " << idx2 << "\n";
                    H(idx1, idx2) += -tmp * double(s(idUnitcell));
                    H(idx2, idx1) += +tmp * double(s(idUnitcell));
                }
            }
        }
    }
};

// class SpinlessTvHoneycombUtils;

// class SpinlessVHoneycombOperator : public SpinlessVOperator
// {
// public:
//     // _s: aux fields, Z_2 variable, length = nUnitcell
//     const SpinlessTvHoneycombUtils *mConfig;
//     const int nUnitcell;

//     SpinlessVHoneycombOperator(const SpinlessTvHoneycombUtils *_config, iVecType *_s, int _bondType, rdGenerator *_rd)
//         :SpinlessVOperator(_config, _s, _bondType, _rd), mConfig(_config), nUnitcell(mConfig->nUnitcell) {
//     }

//     // ~SpinlessVHoneycombOperator() override;
// };


class Honeycomb_tV: public Spinless_tV
{
public:
    const SpinlessTvHoneycombUtils* modelConfig;
    double dt;
    int l;
    int nSites, nUnitcell;
    rdGenerator* rd;

    Honeycomb_tV(SpinlessTvHoneycombUtils* _config, rdGenerator* _rd) {
        modelConfig = _config;
        dt = modelConfig->dt;
        l = modelConfig->l;
        nUnitcell = modelConfig->Lx * modelConfig->Ly;
        nSites = nUnitcell * 2;
        rd = _rd;

        nDim = nSites * 2;
        MatType Ht(nDim, nDim);
        const MatType identity = MatType::Identity(nDim, nDim);
        modelConfig->KineticGenerator(Ht, 1.0);
        MatType expK = expm(Ht, -dt);
        MatType expKhalf = expm(Ht, -dt / 2.0);
        
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