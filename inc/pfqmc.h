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
    void sweep();

    DataType getSign() {
        const MatType identity = MatType::Identity(nDim, nDim);
        const DataType extraSign = ((nDim / 2) % 2 == 0) ? 1.0 : -1.0;
        UDT A(nDim);
        op_array->at(0)->stabilizedLeftMultiply(A);
        MatType gNext, gCur;
        DataType signCur, signNext, signPfaf;
        signCur = op_array->at(0)->getSign();
        A.onePlusInv(nDim, gCur); gCur -= identity;
        for (int i=1; i<op_length; i++) {
            op_array->at(i)->getGreensMat(gNext);
            signNext = op_array->at(i)->getSign();
            // std::cout << gNext << "==== gnext ====\n \n";
            // std::cout << gCur << "==== gcur ====\n \n";
            signPfaf = pfaffianForSignOfProduct(gNext, gCur);
            signCur = (signCur * signNext * signPfaf * extraSign);
            
            // std::cout << signCur << " sign cur\n";
            if (i == op_length-1) break;
            op_array->at(i)->stabilizedLeftMultiply(A);
            A.onePlusInv(nDim, gCur); gCur -= identity;
        }
        return signCur;
    }

    ~PfQMC() {
        for (int i=0; i<udtR.size(); i++) {
            delete udtR[i];
        }
        delete Al;
    }
};



#endif