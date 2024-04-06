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
                (*op_array)[id]->left_multiply(Aseg, Atmp);
                std::swap(Aseg, Atmp);
                id++; 
                if (id >= op_length) break;
            }
            udtR[i] = new UDT(Aseg, nDim);
        }
        
        Al = udtR[0];
        for (int i=1; i<checkpoints; i++) {
            udtR[i]->factorizedMult(*Al, nDim);
        }

    }

    void sweep() {

    }

    //
    void greenFunctionUpdate() {

    }
};



#endif