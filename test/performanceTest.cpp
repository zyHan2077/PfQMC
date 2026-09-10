#include <gtest/gtest.h>
#include "pfqmc.h"
#include "singleMajoranaHoneycomb.h"
#include "honeycomb.h"
#include "square.h"
#include "kitaevChain.h"

// Keep a dense oracle for the optimized interaction operations. The proposal
// formula, RNG and signed acceptance phase are shared intentionally.
class DenseInteraction : public SpinlessVOperator {
public:
    using SpinlessVOperator::SpinlessVOperator;
    void left_multiply(const MatType& a, MatType& out) override { out = B*a; }
    void right_multiply(const MatType& a, MatType& out) override { out = a*B; }
    void left_propagate(MatType& g, MatType& tmp) override {
        reCalcInv(); tmp = B*g; g = tmp*B_inv;
    }
    void right_propagate(MatType& g, MatType& tmp) override {
        reCalcInv(); tmp = B_inv*g; g = tmp*B;
    }
    void stabilizedLeftMultiply(UDT& f) override { f = B*f; }
};

void checkLayer(SpinlessTvUtils& cfg, int bonds, int kind) {
    rdGenerator rd(193);
    auto s = new iVecType(bonds);
    for (int i = 0; i < bonds; ++i) (*s)(i) = rd.rdZ2();
    SpinlessVOperator op(&cfg, s, kind, &rd);
    MatType a = MatType::Random(cfg.nDim, cfg.nDim), out;
    for (int iteration = 0; iteration < 2; ++iteration) {
        MatType expected = op.B*a;
        op.left_multiply(a, out);
        EXPECT_LT((out-expected).norm(), 1e-11);
        out = a; op.left_multiply(out, out);
        EXPECT_LT((out-expected).norm(), 1e-11);
        expected = a*op.B;
        op.right_multiply(a, out);
        EXPECT_LT((out-expected).norm(), 1e-11);
        out = a; op.right_multiply(out, out);
        EXPECT_LT((out-expected).norm(), 1e-11);
        op.reCalcInv();
        MatType tmp;
        out = a; expected = op.B*a*op.B_inv;
        op.left_propagate(out, tmp);
        EXPECT_LT((out-expected).norm(), 1e-10);
        out = a; expected = op.B_inv*a*op.B;
        op.right_propagate(out, tmp);
        EXPECT_LT((out-expected).norm(), 1e-10);
        MatType input = a;
        UDT f(input);
        expected = op.B*a;
        op.stabilizedLeftMultiply(f);
        EXPECT_LT((f.U*f.D.asDiagonal()*f.T-expected).norm(), 1e-10);
        // Force accepted flips and exercise the changed signs of block entries.
        MatType g = 2.0*(MatType::Identity(cfg.nDim,cfg.nDim)+op.B).inverse();
        bool accepted; DataType sign = 1.;
        if (bonds) {
            if (cfg.singleMaj) op.singleFlipSingleMajorana(g, 0, -1., accepted, sign);
            else op.singleFlip(g, 0, -1., accepted, sign);
            EXPECT_TRUE(accepted);
        }
    }
}

TEST(InteractionBlocks, ModelsBoundariesAndChannels) {
    mkl_set_num_threads(1);
    SpinlessTvHoneycombSingleMajoranaUtils majorana(8,8,.1,.7,2);
    SpinlessTvHoneycombUtils honeycomb(3,4,.1,.7,2);
    SpinlessTvSquareUtils square(4,6,.1,.7,2);
    for (int k=0;k<3;++k) { checkLayer(majorana,64,k); checkLayer(honeycomb,12,k); }
    for (int k=0;k<4;++k) checkLayer(square,12,k);
    for (int hs=0;hs<2;++hs) {
        SpinlessTvChainUtils chain(8,.1,.7,2,1,.3,.2,hs);
        checkLayer(chain,4,0); checkLayer(chain,3,1);
        // Overlapping odd-periodic layers must use the dense fallback.
        SpinlessTvChainUtils odd(7,.1,.7,2,0,.3,.2,hs);
        checkLayer(odd,4,0); checkLayer(odd,3,1);
    }
    SpinlessTvHoneycombSingleMajoranaUtils zeroV(3,3,.1,0.,2);
    checkLayer(zeroV,9,0);
}

template<class Config, class Walker>
void compareWalkers(Config& config, int sweeps) {
    rdGenerator ra(923), rb(923);
    Walker wa(&config,&ra), wb(&config,&rb);
    for (auto& base : wb.op_array) {
        auto op = dynamic_cast<SpinlessVOperator*>(base);
        if (op) {
            auto replacement = new DenseInteraction(&config,new iVecType(*op->s),op->bondType,&rb);
            delete base; base = replacement;
        }
    }
    PfQMC a(&wa,7), b(&wb,7); // incomplete final stabilization segment
    for (int t=0;t<sweeps;++t) {
        if (t%2) { a.leftSweep(); b.leftSweep(); }
        else { a.rightSweep(); b.rightSweep(); }
        EXPECT_LT((a.g-b.g).norm()/std::max(1.,b.g.norm()),1e-9);
        EXPECT_LT(std::abs(a.sign-b.sign),1e-9);
        for (int i=0;i<a.op_length;++i) {
            auto sa = wa.op_array[i]->getAuxField();
            if (sa) EXPECT_EQ((*sa-*wb.op_array[i]->getAuxField()).squaredNorm(),0);
        }
    }
    EXPECT_LT(std::abs(a.sign-a.getSignRaw()),1e-8);
}

TEST(InteractionBlocks, SignedSweepTrajectories) {
    mkl_set_num_threads(1);
    SpinlessTvHoneycombSingleMajoranaUtils majorana(8,8,.1,.7,6);
    compareWalkers<SpinlessTvHoneycombSingleMajoranaUtils,HoneycombSingleMajorana_tV>(majorana,6);
    SpinlessTvHoneycombUtils honeycomb(3,3,.1,.7,20);
    compareWalkers<SpinlessTvHoneycombUtils,Honeycomb_tV>(honeycomb,6);
    SpinlessTvSquareUtils square(4,4,.1,.7,6,.3);
    compareWalkers<SpinlessTvSquareUtils,Square_tV>(square,6);
    for(int hs=0;hs<2;++hs) {
        SpinlessTvChainUtils chain(8,.1,.7,20,1,.3,.2,hs);
        compareWalkers<SpinlessTvChainUtils,Chain_tV>(chain,6);
    }
}

TEST(QR_Factorization, PivotScatterScaledColumns) {
    mkl_set_num_threads(1);
    const int n=128;
    MatType a = MatType::Random(n,n);
    for(int j=0;j<n;++j) a.col(j) *= std::pow(10., (j%13)-6);
    MatType saved = a;
    UDT f(a);
    EXPECT_LT((f.U*f.D.asDiagonal()*f.T-saved).norm()/saved.norm(),1e-12);
    EXPECT_LT((f.U.adjoint()*f.U-MatType::Identity(n,n)).norm(),1e-12);
}

TEST(InteractionBlocks, AcceptedNegativeWeightKeepsSign) {
    mkl_set_num_threads(1);
    for (int hs=0; hs<2; ++hs) {
        SpinlessTvChainUtils cfg(4,.1,.7,2,1,.3,.2,hs);
        rdGenerator rd(4);
        SpinlessVOperator op(&cfg,new iVecType(iVecType::Ones(2)),0,&rd);
        MatType g = MatType::Identity(8,8);
        int i,j,k,l;
        cfg.aux2MajoranaIdx(0,0,0,i,j);
        cfg.aux2MajoranaIdx(0,1,0,k,l);
        if (hs==1) j=k;
        g(i,j) = -2.0i/cfg.thlV; g(j,i) = -g(i,j);
        bool accepted=false; DataType sign=1.;
        op.singleFlip(g,0,0.,accepted,sign);
        EXPECT_TRUE(accepted);
        EXPECT_NEAR(std::abs(sign+1.),0.,1e-12);
    }
    SpinlessTvHoneycombSingleMajoranaUtils cfg(2,2,.1,.7,2);
    rdGenerator rd(4);
    SpinlessVOperator op(&cfg,new iVecType(iVecType::Ones(4)),0,&rd);
    MatType g = MatType::Identity(8,8);
    g(0,1) = -2.0i/cfg.thlV; g(1,0) = -g(0,1);
    bool accepted=false; DataType sign=1.;
    op.singleFlipSingleMajorana(g,0,0.,accepted,sign);
    EXPECT_TRUE(accepted);
    EXPECT_NEAR(std::abs(sign+1.),0.,1e-12);
}
