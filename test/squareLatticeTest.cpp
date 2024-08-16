#include <gtest/gtest.h>

#include <ctime>
#include <fstream>

#include "../inc/square.h"
#include "../inc/pfqmc.h"
#include "../inc/pfapack/c_interface/pfapack.h"

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

// TEST(SquareTest, errorSignCatcher) {
//     std::fstream myfile;
//     myfile.open("../doc/1.dat", std::fstream::in);
//     int L = 32;
//     MatType G1(L, L);
//     MatType G2(L, L);
//     DataType tmp;
//     for(int i=0; i<L; i++) {
//         for(int j=0; j<L; j++) {
//             myfile >> G1(i, j);
//         }
//     }

//     for(int i=0; i<L; i++) {
//         for(int j=0; j<L; j++) {
//             myfile >> G2(i, j);
//         }
//         // std::cout << G2(i, 0) << "\n";
//     }
//     myfile.close();

//     G1 = 0.5 * (G1 - G1.transpose());
//     // std::cout << G2 << "\n====G2Raw===\n\n";
//     // std::cout << G2.transpose() << "\n====G2RawTrans===\n\n";
//     // std::cout << 0.5 * (G2 - G2.transpose()) << "\n====G2Raw - Trans===\n\n";
//     MatType G3 = 0.5 * (G2 - G2.transpose());
//     std::cout << G3 << "\n====G3===\n\n";
//     G2 = 0.5 * (G2 - G2.transpose().eval());
//     std::cout << G2 << "\n====G2Skew===\n\n";
//     // G2 = G3;

//     MatType G1copy = G1;
//     MatType G2copy = G2;
//     std::cout << G1 << "\n====G1===\n";
//     std::cout << G2 << "\n====G2Skew===\n";

//     DataType detg2 = G2copy.determinant();
//     DataType pfg2 = pfaf(L/2, G2copy);

//     std::cout << "pfaffian of G1 " << pfaf(L/2, G1copy) << "\n";
//     std::cout << "pfaffian of G2 " << pfg2 << " " << detg2 << " error=" << std::abs(pfg2 * pfg2 - detg2) << "\n";
//     std::cout << "pfaffian of G1+G2 " << pfaffianForSignOfProduct(G1, G2, false) << "\n";
//     // std::cout << G2 << "\n=======\n";
// }

// TEST(PFAPACKTEST, SKTRD) {
//     std::fstream myfile;
//     myfile.open("../test/1.dat", std::fstream::in);
//     int L = 32;
//     MatType G1(L, L);
//     MatType G2(L, L);
//     DataType tmp;
//     for(int i=0; i<L; i++) {
//         for(int j=0; j<L; j++) {
//             myfile >> G1(i, j);
//         }
//     }

//     for(int i=0; i<L; i++) {
//         for(int j=0; j<L; j++) {
//             myfile >> G2(i, j);
//         }
//         // std::cout << G2(i, 0) << "\n";
//     }
//     myfile.close();
//     MatType G2copy = G2;
//     DataType pf = pfaf(L/2, G2copy);
//     // std::cout << "pf G2 = " << pf << "\n";
//     EXPECT_NEAR(std::abs(pf - (831715.0+287.442i)), 0.0, 1.0);

//     DataType signPf = pfaffianForSignOfProduct(G1, G2);
//     EXPECT_NEAR(std::abs(signPf - 1.0), 0.0, 1e-4);
// }