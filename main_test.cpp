#include <gtest/gtest.h>
#include <ctime>
#include "inc/skewMatUtils.h"
#include "inc/honeycomb.h"

#define MPRINTF(...) {printf("[   INFO   ] "); printf(__VA_ARGS__); }

TEST(PfaffianTest, SimpleIntegerMatrix)
{
    int L = 20;
    long long result[] = { -1, 8, -144, 4512, -218912, 
        15247232, -1444779392, 178933241088,
        -28079775960320, 5447517439766528 };

    MatType A(L, L);
    cVecType temp(2 * L);

    // initialize
    int count = 1;
    for (int j = 1; j < L; j++)
    {
        for (int i = 0; i < j; i++)
        {
            A(j, i) = DataType{(double)count, -0.0};
            A(i, j) = DataType{-(double)count, 0.0};
            count++;
        }
    }

    for (int l = 2; l <= L; l += 2) {
        MatType ABlock = A(Eigen::seq(0, l-1), Eigen::seq(0, l-1));
        DataType pf = pfaf(l/2, ABlock, temp);
        double pf0 = result[(l/2)-1];
        EXPECT_NEAR((pf.real() - pf0) / pf0, 0.0, 1e-10) << "failed with l = " << l << "and pf = " << pf << "\n";    
    }
}

TEST(PfaffianTest, RandomEntries)
{
    mkl_set_num_threads(1);

    int N = 100;
    int Nrounds = 100;
    int L = 2 * N;
    srand(114514);
    cVecType temp(4*N);
    DataType pf, det, pf2; 
    std::vector<MatType> matList;
    std::vector<DataType> detList;
    std::vector<DataType> pfList;

    clock_t t1, t2;

    for (int rounds = 0; rounds < Nrounds; rounds++) {
        MatType A = MatType::Random(L, L);
        A = A - A.transpose().eval();
        matList.push_back(A);
    }

    t1 = clock();
    for(auto A: matList) {
        detList.push_back(A.determinant());
    }
    t1 = clock() - t1;

    t2 = clock();
    for (auto A: matList) {
        pfList.push_back( pfaf(N, A, temp) );
    }
    t2 = clock() - t2;

    for (int rounds = 0; rounds < 100; rounds++) {
        pf = pfList.at(rounds);
        det = detList.at(rounds);
        pf2 = pf * pf;
        double r = std::abs(pf2 - det) / std::abs(det);
        // std::cout << pf2 << " " << det << " " << r << "\n";
        EXPECT_NEAR( r, 0.0, 1e-10);
    }

    MPRINTF("complete %d tests with matrix size %d, pfaffian %ld clicks, det %ld clicks\n",
        Nrounds, L, t2, t1);
    MPRINTF("Do not take this result seriously if you are multi-threading\n");
}

TEST(HamiltonianGeneratorTest, honeycombKineticAndInteraction) {
    // Lx, Ly should be >= 2
    // so that the PBC does not double counts 
    int Lx = 11;
    int Ly = 13;
    int nUnitcell = Lx*Ly;
    int hamiltonianDim = Lx*Ly*4;
    SpinlessTvHoneycombUtils config(Lx, Ly, 0.1, 0.7, 10);
    

    // test if Hv1 + Hv2 + Hv3 \propto Ht
    MatType Ht(hamiltonianDim, hamiltonianDim);
    config.KineticGenerator(Ht, config.lambdaV);

    MatType Hv = MatType::Zero(hamiltonianDim, hamiltonianDim);
    iVecType s = iVecType::Constant(nUnitcell, +1);
     // set all +1
    config.InteractionHGenerator(Hv, s, 0);
    config.InteractionHGenerator(Hv, s, 1);
    config.InteractionHGenerator(Hv, s, 2);
    // std::cout << Ht << "\n =1= \n";
    // std::cout << Hv << "\n =2= \n";
    double r = (Ht + Hv).squaredNorm();
    EXPECT_NEAR(r, 0.0, 1e-20);

    
    // test if exp(-Hv) = B
    MatType B;
    rdGenerator rd(114514);
    for(int i=0; i<nUnitcell; i++) s[i] = rd.rdZ2();

    for (int bType=0; bType<3; bType++) {
        Hv = MatType::Zero(hamiltonianDim, hamiltonianDim);
        config.InteractionHGenerator(Hv, s, bType);
        Hv = expm(Hv, -1.0);
        // std::cout << Hv << "\n =2= \n";
        B = MatType::Identity(hamiltonianDim, hamiltonianDim);
        config.InteractionBGenerator(B, s, bType);
        // std::cout << B << "\n == \n";
        r = (Hv - B).squaredNorm();
        // std::cout << r << "\n";
        EXPECT_NEAR(r, 0.0, 1e-20);
    }
}

// test if G = 2 (I + B_1 \cdots B_l)^{-1}
// is updated correctly
// TEST(FastUpdateTest, GreenFunction) {
//     int Lx = 11;
//     int Ly = 13;
//     int nUnitcell = Lx*Ly;
//     int hamiltonianDim = Lx*Ly*4;
//     SpinlessTvHoneycombUtils config(Lx, Ly, 0.1, 0.7, 10);
//     rdGenerator rd(114514);
//     Honeycomb_tV walker(&config, 10, 0.1, 0.7, &rd);
//     MatType g = MatType::Identity(hamiltonianDim, hamiltonianDim);
// }

// test if G_l = 2 (I + B_{l+1} \cdots B_L B_1 \cdots B_l)^{-1}
// is updated correctly
TEST(FastUpdateTest, RatioSquare) {
    // for 8 core, nDim~1600, LTau=20
    // typically ~15 seconds
    mkl_set_num_threads(8);

    int Lx = 19;
    int Ly = 21;
    int LTau = 20;
    int hamiltonianDim = Lx*Ly*4;
    const MatType identity = MatType::Identity(hamiltonianDim, hamiltonianDim);
    SpinlessTvHoneycombUtils config(Lx, Ly, 0.1, 0.7, LTau);
    rdGenerator rd(114514);
    Honeycomb_tV walker(&config, 10, 0.1, 0.7, &rd);
    MatType g = identity;
    MatType A = identity;
    
    int l = 4*LTau - 1;
    for (int i=0; i<4*LTau; i++) {
        // std::cout << i << " here!\n";
        walker.op_array[i]->right_multiply(A, A);
    }
    walker.op_array[4*LTau]->left_multiply(A, A);
    g = 2 * ((g+A).inverse()) - identity;
    
    EXPECT_EQ(walker.op_array[l]->getType(), 2);

    // attempt to flip the first aux field
    iVecType s = *(walker.op_array[l]->getAuxField());
    int auxCur = s(0);
    s(0) = -2 * auxCur;
    for (int i=1; i<s.size(); i++) s(i) = 0;
    MatType Bm = MatType::Zero(hamiltonianDim, hamiltonianDim);
    config.InteractionHGenerator(Bm, s, 2);
    Bm = expm(Bm, -1.0);
    // std::cout << Bm << "\n";


    MatType tmp = (identity + (A*Bm));
    // MatType tmpLU = tmp;
    DataType r2 = logDet(tmp);
    // std::cout << r2 << " " << log(tmp.determinant()) << "==expr2==\n";
    // EXPECT_NEAR(std::abs(r2 - log(tmp.determinant())), 0.0, 1e-10);
    tmp = (identity + A);
    r2 -= logDet(tmp);
    r2 = exp(r2);
    

    DataType r2_1 = (identity + (identity - g)*(Bm - identity)*0.5).determinant();

    // std::cout << r2 << " " << r2_1 << " == r2 and r2_1 ===\n";
    EXPECT_NEAR(std::abs(r2 - r2_1), 0.0, 1e-10);

    
    auto m = config.idxCell2Coord(0);
    // std::cout << config.lambdaV << "\n";
    DataType r = config.etaM;
    for (int imaj = 0; imaj < 2; imaj ++) {
        int idx1 = config.majoranaCoord2Idx(m.ix, m.iy, 0, imaj);
        int idx2 = config.neighborSiteIdx(m.ix, m.iy, imaj, 2);
        // std::cout << idx1 << " " << idx2 << " idx\n";
        // tmp = [1 + i \sigma_{12} \tanh(\lambda / 2) G_{12}]
        r *= ( 1.0 - ( (1.0i) * (config.thlV) * double(auxCur) * g(idx1, idx2) ) );
    }

    // std::cout << r2 << " " << r << " " << r*r << "\n";
    EXPECT_NEAR(std::abs(r*r - r2), 0.0, 1e-10);
}

// test if G = 2 (I + B_1 \cdots B_l)^{-1}
// is updated correctly
TEST(FastUpdateTest, Ratio) {

}