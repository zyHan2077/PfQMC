#include <gtest/gtest.h>
#include <ctime>
#include "inc/skewMatUtils.h"
#include "inc/honeycomb.h"

#define MPRINTF(...) {printf("[   INFO   ] "); printf(__VA_ARGS__); }

TEST(PfafTest, SimpleIntegerMatrix)
{
    int L = 20;
    long long result[] = { -1, 8, -144, 4512, -218912, 
        15247232, -1444779392, 178933241088,
        -28079775960320, 5447517439766528 };

    MatType A(L, L);
    VecType temp(2 * L);

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
        MatType Ablock = A(Eigen::seq(0, l-1), Eigen::seq(0, l-1));
        DataType pf = pfaf(l/2, Ablock, temp);
        double pf0 = result[(l/2)-1];
        EXPECT_NEAR((pf.real() - pf0) / pf0, 0.0, 1e-10) << "failed with l = " << l << "and pf = " << pf << "\n";    
    }
}

TEST(PfafTest, RandomEntries)
{
    mkl_set_num_threads(1);
    
    int N = 100;
    int Nrounds = 100;
    int L = 2 * N;
    srand(114514);
    VecType temp(4*N);
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
    int Lx = 20;
    int Ly = 21;
    int nUnitcell = Lx*Ly;
    int hamiltonianDim = Lx*Ly*4;
    honeycombUtils honeycombLatticeConfig(Lx, Ly);
    
    MatType Ht(hamiltonianDim, hamiltonianDim);
    honeycombLatticeConfig.KineticGenerator(Ht, 1.0);

    MatType Hv = MatType::Zero(hamiltonianDim, hamiltonianDim);
    Eigen::MatrixXi s = Eigen::MatrixXi::Constant(nUnitcell, 3, -1); // set all -1
    honeycombLatticeConfig.InteractionGenerator(Hv, s, 1.0, 0);
    honeycombLatticeConfig.InteractionGenerator(Hv, s, 1.0, 1);
    honeycombLatticeConfig.InteractionGenerator(Hv, s, 1.0, 2);
    double r = (Ht + Hv).squaredNorm();
    EXPECT_NEAR(r, 0.0, 1e-20);
    // std::cout << r << std::endl;
}