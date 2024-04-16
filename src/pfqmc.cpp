#include "pfqmc.h"

PfQMC::PfQMC(Honeycomb_tV *walker, int _stb)
{
    stb = _stb;
    nDim = (walker->nSites) * 2;
    g = MatType::Identity(nDim, nDim);
    op_array = walker->op_array;
    op_length = op_array.size();
    need_stabilization = std::vector<bool>(op_length);
    checkpoints = 0;
    // if a checkpoint is reached
    // stabilized Green's function is re-evaluated
    for (int l = 0; l < op_length; l++)
    {
        bool flag = ((l % stb) == 0);
        need_stabilization[l] = flag;
        if (flag)
        {
            checkpoints++;
        }
    }
    udtL = std::vector<UDT>(checkpoints);
    udtR = std::vector<UDT>(checkpoints);
    leftInit();
    rightInit();
}

void PfQMC::rightSweep()
{
    MatType tmp = MatType::Identity(nDim, nDim);
    MatType Aseg = MatType::Identity(nDim, nDim);
    int curSeg = 0;
    for (int l = 0; l < op_length; l++)
    {
        op_array[l]->update(g);
        op_array[l]->left_multiply(Aseg, tmp);
        std::swap(Aseg, tmp);
        //%op_length is important, cannot be remove. or else the last segment will not be calculated
        if (need_stabilization[(l + 1) % op_length])
        {
            // auto g2 = g;
            // op_array[i]->left_propagate(g2, tmp);
            // re-evaluate the UDT of current segment
            if (curSeg == 0)
            {
                udtR[curSeg] = UDT(Aseg); // TODO: performance check
            }
            else
            {
                udtR[curSeg] = Aseg * udtR[curSeg - 1];
            }
            Aseg = MatType::Identity(nDim, nDim);
            // re-evaluate the Green's function of next time slice
            if (curSeg == (checkpoints - 1))
            {
                udtR[curSeg].onePlusInv(g);
            }
            else
            {
                g = onePlusInv(udtL[curSeg + 1], udtR[curSeg]);
            }
            // std::cout<<"right g recal "<<(g2-g).norm()<<std::endl;
            curSeg++;
        }
        else
        {
            // no need for stabilization
            // direct propagate
            op_array[l]->left_propagate(g, tmp);
        }
    }
}

void PfQMC::leftSweep()
{
    MatType tmp = MatType::Identity(nDim, nDim);
    MatType Aseg = MatType::Identity(nDim, nDim);
    int curSeg = checkpoints - 1;
    for (int l = op_length - 1; l > -1; l--)
    {
        op_array[l]->right_propagate(g, tmp);
        op_array[l]->update(g);
        op_array[l]->right_multiply(Aseg, tmp);
        std::swap(Aseg, tmp);
        if (need_stabilization[l])
        {
            // auto g2 = g;
            Aseg.adjointInPlace();
            // re-evaluate the UDT of current segment
            if (curSeg == (checkpoints - 1))
            {
                udtL[curSeg] = UDT(Aseg); // TODO: performance check
            }
            else
            {
                udtL[curSeg] = Aseg * udtL[curSeg + 1];
            }
            Aseg = MatType::Identity(nDim, nDim);
            // re-evaluate the Green's function of next time slice
            if (curSeg == 0)
            {
                udtL[curSeg].onePlusInv(g);
                g.adjointInPlace();
            }
            else
            {
                g = onePlusInv(udtL[curSeg], udtR[curSeg - 1]);
            }
            curSeg--;
            // std::cout<<"left g recal "<<(g2-g).norm()<<std::endl;
        }
    }
}