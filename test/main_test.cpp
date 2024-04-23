#include <gtest/gtest.h>

#include <ctime>

#include "../inc/honeycomb.h"
#include "../inc/pfqmc.h"
#include "../inc/qr_udt.h"
#include "../inc/skewMatUtils.h"

#define MPRINTF(...)             \
    {                            \
        printf("[   INFO   ] "); \
        printf(__VA_ARGS__);     \
    }

TEST(PfaffianTest, SimpleIntegerMatrix) {
    int L = 20;
    long long result[] = {-1, 8, -144, 4512, -218912,
                          15247232, -1444779392, 178933241088,
                          -28079775960320, 5447517439766528};

    MatType A(L, L);
    // cVecType temp(2 * L);

    // initialize
    int count = 1;
    for (int j = 1; j < L; j++) {
        for (int i = 0; i < j; i++) {
            A(j, i) = DataType{(double)count, -0.0};
            A(i, j) = DataType{-(double)count, 0.0};
            count++;
        }
    }

    for (int l = 2; l <= L; l += 2) {
        MatType ABlock = A(Eigen::seq(0, l - 1), Eigen::seq(0, l - 1));
        DataType pf = pfaf(l / 2, ABlock);
        double pf0 = result[(l / 2) - 1];
        EXPECT_NEAR((pf.real() - pf0) / pf0, 0.0, 1e-10) << "failed with l = " << l << "and pf = " << pf << "\n";
    }
}

TEST(PfaffianTest, RandomEntries) {
    mkl_set_num_threads(1);

    int N = 100;
    int Nrounds = 100;
    int L = 2 * N;
    srand(114514);
    // cVecType temp(4 * N);
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
    for (auto A : matList) {
        detList.push_back(A.determinant());
    }
    t1 = clock() - t1;

    t2 = clock();
    for (auto A : matList) {
        pfList.push_back(pfaf(N, A));
    }
    t2 = clock() - t2;

    for (int rounds = 0; rounds < 100; rounds++) {
        pf = pfList.at(rounds);
        det = detList.at(rounds);
        pf2 = pf * pf;
        double r = std::abs(pf2 - det) / std::abs(det);
        // std::cout << pf2 << " " << det << " " << r << "\n";
        EXPECT_NEAR(r, 0.0, 1e-10);
    }

    MPRINTF("complete %d tests with matrix size %d, pfaffian %ld clicks, det %ld clicks\n",
            Nrounds, L, t2, t1);
    MPRINTF("Do not take this result seriously if you are multi-threading\n");
}

TEST(PfaffianTest, MonteCarloWeight) {
    int Lx = 3;
    int Ly = 3;
    int nUnitcell = Lx * Ly;
    int nDim = Lx * Ly * 4;
    SpinlessTvHoneycombUtils config(Lx, Ly, 0.1, 0.7, 10);

    // test if Hv1 + Hv2 + Hv3 \propto Ht
    MatType Ht(nDim, nDim);
    config.KineticGenerator(Ht, config.lambdaV);
    MatType B = expm(Ht, -1.0);

    DataType eta = pfaffianForEta(Ht);
    DataType w = (MatType::Identity(nDim, nDim) + B).determinant();
    // std::cout << eta << " eta**2 = " << eta*eta << " w= " << w << "\n";
    EXPECT_NEAR(std::abs((eta*eta / w) - 1.0), 0.0, 1e-10);

    DataType signEta = pfaffianForSignOfEta(Ht);
    // std::cout << eta / std::abs(eta) << " signOfEta = " << signEta << "\n";
    // EXPECT_NEAR(std::abs(signEta - eta / std::abs(eta)), 0.0, 1e-30);
    EXPECT_NEAR(std::abs(signEta - 1.0), 0.0, 1e-10);
}

TEST(HamiltonianGeneratorTest, honeycombKineticAndInteraction) {
    // Lx, Ly should be >= 2
    // so that the PBC does not double counts
    int Lx = 19;
    int Ly = 21;
    int nUnitcell = Lx * Ly;
    int hamiltonianDim = Lx * Ly * 4;
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
    const MatType identity = MatType::Identity(hamiltonianDim, hamiltonianDim);
    for (int i = 0; i < nUnitcell; i++) s[i] = rd.rdZ2();

    for (int bType = 0; bType < 3; bType++) {
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

        MatType g = (identity + B).inverse() * 2.0 - identity;
        MatType g1 = MatType::Zero(hamiltonianDim, hamiltonianDim);
        config.InteractionTanhGenerator(g1, s, bType);
        EXPECT_NEAR((g-g1).squaredNorm(), 0.0, 1e-20);
    }
}

// test if G_l = 2 (I + B_{l-1} \cdots B_1 B_L \cdots B_l)^{-1}
// is updated correctly
TEST(FastUpdateTest, GreenFunction) {
    mkl_set_num_threads(8);

    int Lx = 11;
    int Ly = 9;
    int LTau = 20;
    int hamiltonianDim = Lx * Ly * 4;
    const MatType identity = MatType::Identity(hamiltonianDim, hamiltonianDim);
    SpinlessTvHoneycombUtils config(Lx, Ly, 0.1, 0.7, LTau);
    rdGenerator rd(114514);
    Honeycomb_tV walker(&config, &rd);
    MatType g;
    MatType A = identity;

    int l = 4 * LTau - 1;
    for (int i = 0; i < l; i++) {
        walker.op_array[i]->left_multiply(A, A);
    }
    walker.op_array[4 * LTau]->right_multiply(A, A);
    walker.op_array[l]->right_multiply(A, A);

    g = (2 * ((identity + A).inverse()));
    // MatType gcopy = g;
    EXPECT_EQ(walker.op_array[l]->getType(), 2);

    bool flip = false;
    SpinlessVOperator* p = (SpinlessVOperator*) walker.op_array[l];
    for (int i = 0; i < (Lx*Ly)/2; i++) {
        p->singleFlip(g, i, -0.1, flip);
        EXPECT_EQ(flip, true);
    }
    // std::cout << ( g + g.transpose() - (2*identity) ).squaredNorm() << " test g skew\n";

    iVecType s = *(walker.op_array[l]->getAuxField());

    A = identity;
    for (int i = 0; i < l; i++) {
        walker.op_array[i]->left_multiply(A, A);
    }
    walker.op_array[4 * LTau]->right_multiply(A, A);
    walker.op_array[l]->right_multiply(A, A);
    MatType gBrutal = 2 * ((identity + A).inverse());

    double r = (gBrutal - g).squaredNorm();
    // std::cout << "difference between 2 g = " << r << "\n";
    EXPECT_NEAR(r, 0.0, 1e-15);
    // std::cout << (gcopy - g).squaredNorm() << "\n";
}

// test if acceptance ratio R evaluated in fast update
// is in accordance with R^2 given by the determinant formula
TEST(FastUpdateTest, RatioSquare) {
    // for 8 core, nDim~1600, LTau=20
    // typically ~15 seconds
    mkl_set_num_threads(8);

    int Lx = 9;
    int Ly = 11;
    int LTau = 20;
    int hamiltonianDim = Lx * Ly * 4;
    const MatType identity = MatType::Identity(hamiltonianDim, hamiltonianDim);
    SpinlessTvHoneycombUtils config(Lx, Ly, 0.1, 0.7, LTau);
    rdGenerator rd(114514);
    Honeycomb_tV walker(&config, &rd);
    MatType g = identity;
    MatType A = identity;

    int l = 4 * LTau - 1;
    for (int i = 0; i < l; i++) {
        walker.op_array[i]->left_multiply(A, A);
    }
    walker.op_array[4 * LTau]->right_multiply(A, A);
    walker.op_array[l]->right_multiply(A, A);
    g = 2 * ((identity + A).inverse());

    EXPECT_EQ(walker.op_array[l]->getType(), 2);

    // attempt to flip the first aux field
    iVecType s = *(walker.op_array[l]->getAuxField());
    int auxCur = s(0);
    s(0) = -2 * auxCur;
    for (int i = 1; i < s.size(); i++) s(i) = 0;
    MatType Bm = MatType::Zero(hamiltonianDim, hamiltonianDim);
    config.InteractionHGenerator(Bm, s, 2);
    Bm = expm(Bm, -1.0);
    // std::cout << Bm << "\n";

    MatType tmp = (identity + (A * Bm));
    // MatType tmpLU = tmp;
    DataType r2 = logDet(tmp);
    // std::cout << r2 << " " << log(tmp.determinant()) << "==expr2==\n";
    // EXPECT_NEAR(std::abs(r2 - log(tmp.determinant())), 0.0, 1e-10);
    tmp = (identity + A);
    r2 -= logDet(tmp);
    r2 = exp(r2);

    DataType r2_1 = (identity + ((2.0 * identity) - g) * (Bm - identity) * 0.5).determinant();

    // std::cout << r2 << " " << r2_1 << " == r2 and r2_1 ===\n";
    EXPECT_NEAR(std::abs(r2 - r2_1), 0.0, 1e-10);

    auto m = config.idxCell2Coord(0);
    // std::cout << config.lambdaV << "\n";
    DataType r = config.etaM;
    for (int imaj = 0; imaj < 2; imaj++) {
        int idx1 = config.majoranaCoord2Idx(m.ix, m.iy, 0, imaj);
        int idx2 = config.neighborSiteIdx(m.ix, m.iy, imaj, 2);
        // std::cout << idx1 << " " << idx2 << " idx\n";
        // tmp = [1 + i \sigma_{12} \tanh(\lambda / 2) G_{12}]
        r *= (1.0 - ((1.0i) * (config.thlV) * double(auxCur) * g(idx1, idx2)));
    }

    // std::cout << r2 << " " << r << " " << r*r << "\n";
    EXPECT_NEAR(std::abs(r * r - r2), 0.0, 1e-10);
}

TEST(QR_Factorization, UDT_Decomposition) {

    mkl_set_num_threads(8);

    srand(114514);
    int nDim = 1000;
    MatType A1 = MatType::Random(nDim, nDim);
    MatType A2 = MatType::Random(nDim, nDim);
    // std::cout << "A=\n" << A << "\n";
    MatType A1copy = A1;
    MatType A2copy = A2;
    UDT F1(A1);
    // std::cout << "Aafter=\n" << A << "\n";
    EXPECT_NE(A1.data(), A1copy.data());

    EXPECT_NEAR(((F1.U * F1.D.asDiagonal() * F1.T) - A1copy).squaredNorm(), 0.0, 1e-20);

    UDT F2(A2);
    F2 = F1 * F2;
    // F1.factorizedMultUpdate(F2, nDim);
    double r = ((F2.U * F2.D.asDiagonal() * F2.T) - (A1copy * A2copy)).squaredNorm();
    // std::cout << "r=" << r << "\n";
    EXPECT_NEAR(r, 0.0, 1e-18);

    // Green function test
    MatType g;
    MatType gBrutal = MatType::Identity(nDim, nDim);
    gBrutal = (gBrutal + ((F2.U * F2.D.asDiagonal()) * F2.T)).inverse() * 2.0;
    F2.onePlusInv(g);
    // std::cout << g << "\n====\n";
    // std::cout << gBrutal << "\n";
    // std::cout << (g - gBrutal).squaredNorm() << "\n";
    EXPECT_NEAR((g - gBrutal).squaredNorm(), 0.0, 1e-15);
}

TEST(QR_Factorization, BMult) {
    srand(114514);
    int nDim = 1000;
    MatType A1 = MatType::Random(nDim, nDim);
    UDT id(nDim);
    EXPECT_NEAR( (id.U * id.D.asDiagonal() * id.T - 
                MatType::Identity(nDim, nDim)).squaredNorm(), 0.0, 
                1e-20);
    MatType A1copy = A1;
    UDT F1(A1);
    EXPECT_NEAR(((F1.U * F1.D.asDiagonal() * F1.T) - A1copy).squaredNorm(), 0.0, 1e-20);

    MatType A2 = MatType::Random(nDim, nDim);
    // A2 = MatType::Identity(nDim, nDim);
    F1 = A2 * F1;
    // F1.bMultUpdate(A2);
    EXPECT_NEAR(((F1.U * F1.D.asDiagonal() * F1.T) - (A2 * A1copy)).squaredNorm(), 0.0, 1e-20);
}

TEST(PFQMC, GetSignTest) {
    mkl_set_num_threads(8);

    int Lx = 9;
    int Ly = 11;
    int LTau = 9;
    int hamiltonianDim = Lx * Ly * 4;
    const MatType identity = MatType::Identity(hamiltonianDim, hamiltonianDim);
    SpinlessTvHoneycombUtils config(Lx, Ly, 0.1, 0.7, LTau);
    rdGenerator rd(42);
    Honeycomb_tV walker(&config, &rd);
    PfQMC pfqmc(&walker, 10);
    //TODO: getSign optimization
    DataType sRaw = pfqmc.getSignRaw();
    // DataType s = pfqmc.getSign();
    // std::cout << s << "\n";
    EXPECT_NEAR(std::abs(sRaw-1.0), 0.0, 1e-10);
    // EXPECT_NEAR(std::abs(s-1.0), 0.0, 1e-10);
}

TEST(PFQMC, udtRTest) {
    mkl_set_num_threads(8);

    int Lx = 2;
    int Ly = 2;
    int LTau = 10;
    int hamiltonianDim = Lx * Ly * 4;
    const MatType identity = MatType::Identity(hamiltonianDim, hamiltonianDim);
    SpinlessTvHoneycombUtils config(Lx, Ly, 0.1, 0.7, LTau);
    rdGenerator rd(114514);
    Honeycomb_tV walker(&config, &rd);
    PfQMC pfqmc(&walker, 10);

    int thermalLength = 1000;

    for (int i = 0; i < thermalLength; i++) {
        pfqmc.leftSweep();
    }

    MatType A = identity;
    for (int i = 0; i < pfqmc.op_length; i++) {
        walker.op_array[i]->left_multiply(A, A);
    }
    MatType gBrutal = 2.0 * (identity + A).inverse();
    EXPECT_NEAR((gBrutal - pfqmc.g).squaredNorm(), 0.0, 1e-10);
}
