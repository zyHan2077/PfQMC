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

void PfQMC::sweep() {
    MatType tmp, Aseg;
    int curSeg = 0;
    int j;
    for(int i=0; i<op_length; i++) {
        op_array->at(i)->update(g);
        if (need_stabilization[(i+1) % op_length]) {
            // re-evaluate the UDT of current segment
            Aseg = MatType::Identity(nDim, nDim);
            for (j = i; !need_stabilization[j]; j--) {
                op_array->at(j)->right_multiply(Aseg, tmp);
                std::swap(Aseg, tmp);
            }
            op_array->at(j)->right_multiply(Aseg, tmp);
            std::swap(Aseg, tmp);

            delete udtR[curSeg];
            udtR[curSeg] = new UDT(Aseg, nDim); // TODO: performance check

            // re-evaluate the Green's function of next time slice
            curSeg = (curSeg + 1) % checkpoints;
            // std::cout << i << "=i, " << curSeg << "=curSeg \n";

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