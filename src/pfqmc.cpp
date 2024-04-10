#include "pfqmc.h"

PfQMC::PfQMC(Honeycomb_tV *walker, int _stb) {
    stb = _stb;
    nDim = (walker->nSites) * 2;
    g = MatType::Identity(nDim, nDim);
    op_array = &(walker->op_array);
    op_length = op_array->size();

    need_stabilization = std::vector<bool>(op_length);
    checkpoints = 0;
    bool flag;

    // if a checkpoint is reached
    // stabilized Green's function is re-evaluated
    for (int i = 0; i < op_length; i++) {
        flag = ((i % stb) == 0);
        need_stabilization[i] = flag;
        if (flag) {
            checkpoints++;
        }
    }
    udtR = std::vector<UDT*>(checkpoints);
    rightInit();
}