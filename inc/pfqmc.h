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

    DataType getSign()
    {
        const MatType identity = MatType::Identity(nDim, nDim);
        const DataType extraSign = ((nDim / 2) % 2 == 0) ? 1.0 : -1.0;
        UDT A(nDim);
        op_array[0]->stabilizedLeftMultiply(A);
        MatType gNext, gCur;
        DataType signCur, signNext, signPfaf;
        signCur = op_array[0]->getSign();
        A.onePlusInv(gCur);
        gCur -= identity;
        for (int i = 1; i < op_length; i++)
        {
            op_array[i]->getGreensMat(gNext);
            signNext = op_array[i]->getSign();
            // std::cout << gNext << "==== gnext ====\n \n";
            // std::cout << gCur << "==== gcur ====\n \n";
            signPfaf = pfaffianForSignOfProduct(gNext, gCur);
            signCur = (signCur * signNext * signPfaf * extraSign);

            // std::cout << signCur << " sign cur\n";
            if (i == op_length - 1)
                break;
            op_array[i]->stabilizedLeftMultiply(A);
            A.onePlusInv(gCur);
            gCur -= identity;
        }
        return signCur;
    }
    
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