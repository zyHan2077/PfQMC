#include <gtest/gtest.h>

#include <ctime>

#include "../inc/square.h"
#include "../inc/pfqmc.h"

TEST(SquareTest, Config) {
    int Lx = 4;
    int Ly = 4;
    int LTau = 9;
    int nDim = Lx * Ly * 2;
    const MatType identity = MatType::Identity(nDim, nDim);
    SpinlessTvSquareUtils config(Lx, Ly, 0.1, 0.7, LTau);
    MatType B = identity;
    iVecType s = iVecType::Constant(Lx*Ly/2, +1);
    config.InteractionBGenerator(B, s, 0);
    rdGenerator rd(42);
    Square_tV walker(&config, &rd);
    PfQMC pfqmc(&walker, 10);
    for (int i = 0; i < 10; i++)
    {
        pfqmc.rightSweep();
        pfqmc.leftSweep();
    }

}