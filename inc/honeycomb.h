#include "types.h"
#include "operator.h"

// a honeycomb lattice is defined by 
// - Lx * Ly unit cells
// - 2 sites per unit cell
// - 3 bonds per unit cell
// - 2 majorana species
class HoneycombUtils
{
public:
    struct majoranaCoord
    {
        int ix, iy, isubcell, imaj;
        majoranaCoord(int _ix=0, int _iy=0, int _isubcell=0, int _imaj=0)
            : ix(_ix), iy(_iy), isubcell(_isubcell), imaj(_imaj) {};
    };
    
    int Lx, Ly;
    int nsites;

    HoneycombUtils(int _Lx, int _Ly) {
        Lx = _Lx;
        Ly = _Ly;
        nsites = Lx*Ly*2;
    }

    inline int unitCellCoord2Idx(int ix, int iy) {
        return ix*Ly + iy;
    }
    
    inline int majoranaCoord2Idx(majoranaCoord s) {
        return nsites*s.imaj + unitCellCoord2Idx(s.ix, s.iy)*2 + s.isubcell;
    }

    inline majoranaCoord idx2MajoranaCoord(int idx) {
        majoranaCoord s;
        s.imaj = idx / nsites;
        int tmp = idx % nsites;
        s.isubcell = tmp % 2;
        tmp = tmp / 2;
        s.iy = tmp % Ly;
        s.ix = tmp / Ly;
        return s;
    }

    inline void KineticGenerator(MatType &H, DataType t) {
        H.setZero();
        int idx1, idx2;
        DataType tmp = (0.25i) * t;
        for (int i=0; i<Lx; i++) {
            for (int j=0; j<Ly; j++) {
                for (int k=0; k<2; k++) {
                    idx1 = majoranaCoord2Idx(majoranaCoord(i, j, 0, k));
                    idx2 = majoranaCoord2Idx(majoranaCoord(i, j, 1, k));
                    H(idx1, idx2) = tmp;
                    H(idx2, idx1) = -tmp;
                    idx2 = majoranaCoord2Idx(majoranaCoord((i+1) % Lx, j, 1, k));
                    H(idx1, idx2) = tmp;
                    H(idx2, idx1) = -tmp;
                    idx2 = majoranaCoord2Idx(majoranaCoord(i, (j+1) % Ly, 1, k));
                    H(idx1, idx2) = tmp;
                    H(idx2, idx1) = -tmp;
                }
            }
        }
    }

    inline void InteractionGenerator(MatType &H, const Eigen::MatrixXi &s, DataType lambda, const int bondType) {
        DataType tmp = (0.25i) * lambda;
        // int idxGen(int, int, int);
        // switch (bondType)
        // {
        // case 0:
        //     idxGen = [&] (int i, int j, int k) 
        //     { return majoranaCoord2Idx(majoranaCoord(i, j, 1, k)); };
        //     break;
        // case 1:
        //     idxGen = [&] (int i, int j, int k) 
        //     { return majoranaCoord2Idx(majoranaCoord((i+1) % Lx, j, 1, k)); };
        //     break;
        // default:
        //     idxGen = [&] (int i, int j, int k) 
        //     { return majoranaCoord2Idx(majoranaCoord(i, (j+1) % Ly, 1, k)); };
        //     break;
        // }
        auto idxGen = [&] (int i, int j, int k) {
            if (bondType == 0) {
                return majoranaCoord2Idx(majoranaCoord(i, j, 1, k));
            } else if (bondType == 1) {
                return majoranaCoord2Idx(majoranaCoord((i+1) % Lx, j, 1, k));
            } else {
                return majoranaCoord2Idx(majoranaCoord(i, (j+1) % Ly, 1, k));
            }
        };
        int idx1, idx2, idUnitcell;
        for (int i=0; i<Lx; i++) {
            for (int j=0; j<Ly; j++) {
                idUnitcell = unitCellCoord2Idx(i, j);
                for (int k=0; k<2; k++) {
                    idx1 = majoranaCoord2Idx(majoranaCoord(i, j, 0, k));
                    idx2 = idxGen(i, j, k);
                    H(idx1, idx2) = tmp * double(s(idUnitcell, bondType));
                    H(idx2, idx1) = -tmp * double(s(idUnitcell, bondType));
                }
            }
        }
    }
};



class Honeycomb_tV 
{
    HoneycombUtils latticeConfig;
    double dt;
    double V;
    int l; // time slices
    
    std::vector<Operator *> op_array;
    
    
};