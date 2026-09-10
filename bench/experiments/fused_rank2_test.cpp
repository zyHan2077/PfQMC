
TEST(InteractionBlocks, FusedRankTwoMatchesUnconjugatedBlas) {
    mkl_set_num_threads(1);
    for (int n : {8,128,288}) {
        MatType expected = MatType::Random(n,n), actual=expected;
        cVecType x=cVecType::Random(n), y=cVecType::Random(n);
        DataType alpha(.3,-.7), minus=-alpha;
        const int inc=1;
        zgeru(&n,&n,&alpha,x.data(),&inc,y.data(),&inc,expected.data(),&n);
        zgeru(&n,&n,&minus,y.data(),&inc,x.data(),&inc,expected.data(),&n);
        skewRank2Update(actual,x,y,alpha);
        EXPECT_LT((actual-expected).norm()/expected.norm(),1e-14);
    }
}
