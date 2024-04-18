#ifndef PFQMC_H
#define PFQMC_H
#include "honeycomb.h"
#include "qr_udt.h"

class PfQMC
{
public:
    int stb;
    int nDim;
    MatType g;
    std::vector<Operator *> op_array;
    int op_length;
    std::vector<bool> need_stabilization;
    int checkpoints;
    std::vector<UDT> udtL;
    std::vector<UDT> udtR;

    PfQMC(Honeycomb_tV *walker, int _stb = 10);

    void rightInit()
    {
        MatType tmp = MatType::Identity(nDim, nDim);
        MatType Aseg = MatType::Identity(nDim, nDim);
        int curSeg = 0;
        for (int l = 0; l < op_length; l++)
        {
            op_array[l]->left_multiply(Aseg, tmp);
            std::swap(Aseg, tmp);
            // %op_length is important, cannot be remove. or else the last segment will not be calculated
            if (need_stabilization[(l + 1) % op_length])
            {
                if (curSeg == 0)
                {
                    udtR[curSeg] = UDT(Aseg); // TODO: performance check
                }
                else
                {
                    udtR[curSeg] = Aseg * udtR[curSeg - 1];
                }
                Aseg = MatType::Identity(nDim, nDim);
                curSeg++;
            }
        }
        udtR[checkpoints - 1].onePlusInv(g);
    }

    void leftInit()
    {
        MatType tmp = MatType::Identity(nDim, nDim);
        MatType Aseg = MatType::Identity(nDim, nDim);
        int curSeg = checkpoints - 1;
        for (int l = op_length - 1; l > -1; l--)
        {
            op_array[l]->right_multiply(Aseg, tmp);
            std::swap(Aseg, tmp);
            if (need_stabilization[l])
            {
                Aseg.adjointInPlace();
                if (curSeg == (checkpoints - 1))
                {
                    udtL[curSeg] = UDT(Aseg); // TODO: performance check
                }
                else
                {
                    udtL[curSeg] = Aseg * udtL[curSeg + 1];
                }
                Aseg = MatType::Identity(nDim, nDim);
                curSeg--;
            }
        }
        udtL[0].onePlusInv(g);
        g.adjointInPlace();
    }

    // after each sweep
    // the greens function g
    // is automatically updated
    void rightSweep();
    void leftSweep();

    // get sign by computing the pfaffian of
    // a 4N * 4N matrix, for testing purpose
    DataType getSignRaw();

    // should provide same result as getSignRaw
    // but by computing the pfaffian of a 2N * 2N matrix
    DataType getSign();
    
    // ~PfQMC()
    // {
        // for (int i = 0; i < udtR.size(); i++)
        // {
        //     delete udtR[i];
        // }
        // delete Al;
    // }
};

#endif