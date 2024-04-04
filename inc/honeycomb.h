#ifndef Honeycomb_H
#define Honeycomb_H

#include "types.h"
#include "operator.h"
#include "skewMatUtils.h"

// a honeycomb lattice is defined by 
// - Lx * Ly unit cells
// - 2 sites per unit cell
// - 3 bonds per unit cell
// - 2 majorana species
class SpinlessTvHoneycombUtils
{
public:
    struct majoranaCoord
    {
        int ix, iy, iSubcell, imaj;
        majoranaCoord(int _ix=0, int _iy=0, int _iSubcell=0, int _imaj=0)
            : ix(_ix), iy(_iy), iSubcell(_iSubcell), imaj(_imaj) {};
    };
    
    int Lx, Ly;
    int nsites;
    int nUnitcell; // also number of auxillary fields per direction
    double dt;
    double V;
    int l;
    double lambdaV, chlV, shlV, thlV, etaM;

    SpinlessTvHoneycombUtils(int _Lx, int _Ly, double _dt, 
            double _V, int _l) {
    
        // model configuration
        Lx = _Lx;
        Ly = _Ly;
        nUnitcell = Lx*Ly;
        nsites = nUnitcell*2;
        dt = _dt;
        V = _V;
        l = _l; // imaginary time slices

        lambdaV = acosh(exp(0.5*V*dt));
        chlV = cosh(lambdaV);
        shlV = sinh(lambdaV);
        thlV = tanh(lambdaV);
        etaM = chlV * chlV;
        
    }

    inline int unitCellCoord2Idx(int ix, int iy) const {
        return ix*Ly + iy;
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

    // Directly generate B by directly writing each 2*2 block
    // B should be initialized as Identity
    inline void InteractionBGenerator(MatType &B, const iVecType &s, const int bondType) const {
        DataType ch = chlV;
        DataType ish = (1.0i) * shlV;

        int idx1, idx2, idUnitcell;
        for (int i=0; i<Lx; i++) {
            for (int j=0; j<Ly; j++) {
                idUnitcell = unitCellCoord2Idx(i, j);
                for (int k=0; k<2; k++) {
                    idx1 = majoranaCoord2Idx(i, j, 0, k);
                    idx2 = neighborSiteIdx(i, j, k, bondType);
                    // std::cout << "(" << idx1 << "," << idx2 << "\n";
                    //  |   \cosh(\lambda)           ,  -i \sinh(\lambda) \sigma |
                    //  |  +i \sinh(\lambda) \sigma  ,   \cosh(\lambda)          |
                    B(idx1, idx1) = ch;
                    B(idx2, idx2) = ch;
                    B(idx1, idx2) = +ish * double(s(idUnitcell));
                    B(idx2, idx1) = -ish * double(s(idUnitcell));
                }
            }
        }
    }
};



class Honeycomb_tV 
{
public:
    const SpinlessTvHoneycombUtils* modelConfig;
    
    std::vector<Operator *> op_array;
    double dt;
    int l;
    int nSites, nUnitcell;
    rdGenerator* rd;

    Honeycomb_tV(SpinlessTvHoneycombUtils* _config, rdGenerator* _rd) {
        modelConfig = _config;
        dt = modelConfig->dt;
        l = modelConfig->l;
        nSites = modelConfig->nsites;
        nUnitcell = modelConfig->nUnitcell;
        rd = _rd;

        int ndim = nSites * 2;
        MatType Ht(ndim, ndim);
        modelConfig->KineticGenerator(Ht, 1.0);
        MatType expK = expm(Ht, -dt);
        MatType expKhalf = expm(Ht, -dt / 2.0);
        auto s = new iVecType[3*l];
        for (int i = 0; i < 3 * l; i++)
        {
            s[i] = iVecType(nUnitcell);
            for (int j=0; j<nUnitcell; j++) s[i](j) = rd->rdZ2();
        }
        
        op_array = std::vector<Operator*>(4*l + 1);
        for (int i=0; i<l; i++) {
            if (i == 0) {
                op_array[0] = new DenseOperator(expKhalf);
            } else {
                op_array[4*i] = new DenseOperator(expK);
            }
            op_array[4*i + 1] = new SpinlessVOperator(modelConfig, &(s[3*i]), 0, rd);
            op_array[4*i + 2] = new SpinlessVOperator(modelConfig, &(s[3*i+1]), 1, rd);
            op_array[4*i + 3] = new SpinlessVOperator(modelConfig, &(s[3*i+2]), 2, rd);
        }
        op_array[4*l] = new DenseOperator(expKhalf);
    }

    ~Honeycomb_tV() {
        for (int i=0; i<4l+1; i++) {
            delete op_array[i];
        }
    }
    
};

#endif