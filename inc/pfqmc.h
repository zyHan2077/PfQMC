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

    PfQMC(Honeycomb_tV *walker, int _stb = 10);
    ~PfQMC();
};

PfQMC::PfQMC(Honeycomb_tV *walker, int _stb) {
    stb = _stb;
    nDim = (walker->nSites) * 2;
    g = MatType::Identity(nDim, nDim);
    op_array = &(walker->op_array);
    op_length = op_array->size();

    need_stabilization = std::vector<bool>(op_length + 1);
    checkpoints = 0;
    bool flag;
    for (int i = 0; i < (op_length + 1); i++) {
        flag = ((i % stb) == 0) | (i == op_length);
        need_stabilization[i] = flag;
        if (flag) {
            checkpoints++;
        }
    }


}

PfQMC::~PfQMC() {}

#endif