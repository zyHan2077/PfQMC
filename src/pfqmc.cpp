#include "pfqmc.h"

PfQMC::PfQMC(Spinless_tV *walker, int _stb)
{
    stb = _stb;
    nDim = walker->nDim;
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
    sign = getSignRaw();
}

void PfQMC::rightSweep()
{
    MatType tmp = MatType::Identity(nDim, nDim);
    MatType Aseg = MatType::Identity(nDim, nDim);
    int curSeg = 0;
    DataType signCur;
    for (int l = 0; l < op_length; l++)
    {
        signCur = op_array[l]->update(g);
        this->sign *= signCur;

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
    DataType signCur;
    for (int l = op_length - 1; l > -1; l--)
    {
        op_array[l]->right_propagate(g, tmp);
        signCur = op_array[l]->update(g);
        sign *= signCur;
        
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

DataType PfQMC::getSignRaw()
{
    const MatType identity = MatType::Identity(nDim, nDim);
    const DataType extraSign = ((nDim / 2) % 2 == 0) ? 1.0 : -1.0;
    UDT A(nDim);
    op_array[0]->stabilizedLeftMultiply(A);
    MatType gNext, gCur;
    DataType signCur, signNext, signPfaf;
    signCur = op_array[0]->getSignOfWeight();
    A.onePlusInv(gCur);
    gCur -= identity;
    for (int i = 1; i < op_length; i++)
    {
        op_array[i]->getGreensMat(gNext);
        signNext = op_array[i]->getSignOfWeight();
        // std::cout << gNext << "==== gnext ====\n \n";
        // std::cout << gCur << "==== gcur ====\n \n";
        signPfaf = pfaffianForSignOfProduct(gNext, gCur);
        // std::cout << "Raw sCur, sNext, sPfaf=" << signCur << " " << signNext << " " << signPfaf << "\n"; 
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

// DataType PfQMC::getSign() {
//     //TODO: this current method has fundamental flaws
//     const MatType identity = MatType::Identity(nDim, nDim);
//     UDT A(nDim);
//     op_array[0]->stabilizedLeftMultiply(A);
//     MatType gNext, gCur, gInv, gTemp, t;
//     DataType signCur, signNext, signPfaf;
//     signCur = op_array[0]->getSignOfWeight();
//     A.onePlusInv(gCur);
//     gCur -= identity;
//     for (int i = 1; i < op_length; i++)
//     {
//         op_array[i]->getGreensMat(gNext);
//         // std::cout << gNext << "==== gnext ====\n \n";
//         // std::cout << "det of g=" << gNext.determinant() << "\n";
//         op_array[i]->getGreensMatInv(gInv);
//         MatType t = gNext * gInv;
//         // std::cout << (t - identity).squaredNorm() << "gNext * gInv, bType" << op_array[i]->getType()<<  " \n";
//         // EXPECT_NEAR(, 0.0, 1e-10);
//         signNext = op_array[i]->getSignOfWeight();
//         // std::cout << gNext << "==== gnext ====\n \n";
//         // std::cout << gCur << "==== gcur ====\n \n";
//         // std::cout << gInv << "=gInv\n";
//         gTemp = gInv + gCur;
//         signPfaf = signOfPfaf(gTemp) / op_array[i]->getSignPfGInv();
//         // signPfaf = 1.0;
//         // std::cout << "sCur, sNext, sPfaf=" << signCur << " " << signNext << " " << signPfaf << "\n\n"; 
//         signCur = (signCur * signNext * signPfaf);
//         // std::cout << signCur << " sign cur\n";
//         if (i == op_length - 1)
//             break;
//         op_array[i]->stabilizedLeftMultiply(A);
//         A.onePlusInv(gCur);
//         gCur -= identity;
//     }
//     return signCur;
// }