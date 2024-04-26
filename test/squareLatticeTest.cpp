#include <gtest/gtest.h>

#include <ctime>

#include "../inc/square.h"
#include "../inc/pfqmc.h"

TEST(SinhTest, SinhH) {
    int nDim = 20;
    MatType H = MatType::Random(nDim, nDim);
    MatType s1 = sinhHQuarterSqrt2(H);
    MatType s2 = (expm(H, 0.25) - expm(H, -0.25)) * sqrt(2) / 2.0;
    EXPECT_NEAR((s1 - s2).squaredNorm(), 0.0, 1e-10);
}

TEST(SquareTest, Config) {
    int Lx = 10;
    int Ly = 10;
    int LTau = 10;
    int nDim = Lx * Ly * 2;
    const MatType identity = MatType::Identity(nDim, nDim);
    SpinlessTvSquareUtils config(Lx, Ly, 0.1, 0.7, LTau);
    MatType B = identity;
    // iVecType s = iVecType::Constant(Lx*Ly/2, +1);
    // config.InteractionBGenerator(B, s, 0);
    rdGenerator rd(42);
    Square_tV walker(&config, &rd);
    PfQMC pfqmc(&walker, 10);
    DataType sign = pfqmc.getSignRaw();
    EXPECT_NEAR(std::abs(sign-1.0), 0.0, 1e-10);
}

TEST(SquareTest, SignTest) {
    int Lx = 12;
    int Ly = 10;
    int LTau = 10;
    double dt = 0.2;
    int nDim = Lx * Ly * 2;
    SpinlessTvSquareUtils config(Lx, Ly, 0.1, 0.7, LTau, 0.5);
    MatType Ht(nDim, nDim);
    config.KineticGenerator(Ht);
    MatType Hcopy = dt * Ht;
    DataType signK = signOfHamiltonian(Hcopy);
    Hcopy = (0.5 * dt ) * Ht;
    DataType signKHalf = signOfHamiltonian(Hcopy);
    // std::cout << "signK = " << signK << " signKhalf = " << signKHalf << "\n";
    EXPECT_NEAR(std::abs(signK - 1.0), 0.0, 1e-10);
    EXPECT_NEAR(std::abs(signKHalf - 1.0), 0.0, 1e-10);
}