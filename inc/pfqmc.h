#ifndef PFQMC_H
#define PFQMC_H

#include "honeycomb.h"
#include "qr_udt.h"

class PfQMC {
   public:
    int stb;
    int nDim;
    MatType g;
    std::vector<Operator *> *op_array;
    int op_length;
    std::vector<bool> need_stabilization;
    int checkpoints;
    std::vector<UDT*> udtR;
    UDT* Al;

    PfQMC(Honeycomb_tV *walker, int _stb = 10);

    void rightInit() {
        int id = 0;
        MatType Aseg, Atmp;
        for (int i=0; i<checkpoints; i++) {
            Aseg = MatType::Identity(nDim, nDim);
            for (int j=0; j<stb; j++) {
                op_array->at(id)->left_multiply(Aseg, Atmp);
                std::swap(Aseg, Atmp);
                id++; 
                if (id >= op_length) break;
            }
            udtR[i] = new UDT(Aseg, nDim);
        }
        
        // initialize green function g
        Al = new UDT(*udtR[0]);
        for (int i=1; i<checkpoints; i++) {
            udtR[i]->factorizedMultUpdate(*Al, nDim);
        }
        Al->onePlusInv(nDim, g);
    }

    // after each sweep
    // the greens function g
    // is automatically updated
    void sweep() {
        MatType tmp, Aseg;
        int curSeg = 0;
        for(int i=0; i<op_length; i++) {
            op_array->at(i)->update(g);
            if (need_stabilization[(i+1) % op_length]) {
                
                // re-evaluate the UDT of current segment
                Aseg = MatType::Identity(nDim, nDim);
                for (int j = (i - stb + 1); (j<=i) && j<op_length; j++) {
                    op_array->at(j)->left_multiply(Aseg, tmp);
                    std::swap(Aseg, tmp);
                }

                delete udtR[curSeg];
                udtR[curSeg] = new UDT(Aseg, nDim); // TODO: performance check

                // re-evaluate the Green's function of next time slice
                curSeg = (curSeg + 1) % checkpoints;

                delete Al;
                Al = new UDT(*udtR[curSeg]);
                for (int k= curSeg + 1; k<checkpoints; k++) {
                    udtR[k]->factorizedMultUpdate(*Al, nDim);
                }
                for (int k = 0; k<curSeg; k++) {
                    udtR[k]->factorizedMultUpdate(*Al, nDim);
                }
                Al->onePlusInv(nDim, g);

            } else {
                // no need for stabilization
                // direct propagate
                op_array->at(i)->left_propagate(g, tmp);
            }
        }
    }

    DataType getSign() {
        return 1.0;
    }

    ~PfQMC() {
        for (int i=0; i<udtR.size(); i++) {
            delete udtR[i];
        }
        delete Al;
    }
};



#endif