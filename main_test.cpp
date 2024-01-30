#include <gtest/gtest.h>
#include "src/skewMatUtils.h"


TEST(PfafTest, ResultError) {
    int N = 5;

    dtype *A, *temp1;
    // dtype *temp2;
    A = (dtype *) mkl_malloc((4*N*N)*sizeof(dtype), 64);
    temp1 = (dtype *) mkl_malloc((4*N)*sizeof(dtype), 64);
    // temp2 = (dtype *) mkl_malloc((2*N)*sizeof(dtype), 64);

    // initialize
    int count = 1;
    for(int i=0; i<2*N; i++) {
        for (int j=i+1; j<2*N; j++) {
            A[i*2*N + j] = dtype{(double)count, -0.0};
            A[j*2*N + i] = dtype{-(double)count, 0.0};
            count++;
        }
        // temp1[i] = count;
    }

    // printMat(2*N, A);

    dtype pfaf;
    SkewMatHouseholder_PureMKL(N, A, temp1, &pfaf);
    
    // printMat(2*N, A);

    // std::cout << "pfaffian= " << (pfaf) << "\n";

    EXPECT_NEAR(pfaf.real(), -2359296, 1e-5) << "failed with error 1111 " << (pfaf+2359296.0);

    MKL_free(A);MKL_free(temp1);
    // MKL_free(temp2);
}

TEST(PfafTest, RandomEntries) {
    int N = 10;
    int sizemat = 4*N*N;
    int tempind;
    dtype tempz;
    VSLStreamStatePtr stream;
    double *x1, *x2;

    x1 = (double *) mkl_malloc(sizemat*sizeof(double), 64);
    x2 = (double *) mkl_malloc(sizemat*sizeof(double), 64);

    vslNewStream(&stream, VSL_BRNG_MCG59,  114514);
    vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, sizemat, x1, -1, 1);
    vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, sizemat, x2, -1, 1);

    // printVec(sizemat, x1);
    // printVec(sizemat, x2);
    

    dtype *A, *temp1;
    // dtype *temp2;
    A = (dtype *) mkl_malloc((4*N*N)*sizeof(dtype), 64);
    temp1 = (dtype *) mkl_malloc((4*N)*sizeof(dtype), 64);
    // temp2 = (dtype *) mkl_malloc((2*N)*sizeof(dtype), 64);

    // initialize
    // int count = 1;
    for(int i=0; i<2*N; i++) {
        for (int j=i+1; j<2*N; j++) {
            tempind = i*2*N + j;
            tempz = x1[tempind] + (im*x2[tempind]);
            A[tempind] = tempz;
            A[j*2*N + i] = - tempz;
            // count++;
        }
        A[i*2*N + i] = 0;
        // temp1[i] = count;
    }

    printMat(2*N, A);

    Eigen::MatrixXcd Acopy(2*N, 2*N);
    memcpy(Acopy.data(), A, (4*N*N)*sizeof(dtype));
    std::cout << Acopy << std::endl;

    dtype pfaf;
    SkewMatHouseholder_PureMKL(N, A, temp1, &pfaf);

    std::cout << "pfaf complete! = " << pfaf << std::endl;
    dtype det;
    // lapack_int ipiv[2*N];
    det = Acopy.determinant();
    
    // printMat(2*N, A);

    std::cout << "det= " << (det) << "\n";

    pfaf = pfaf*pfaf;

    EXPECT_NEAR(pfaf.real(), det.real(), 1e-5) << "failed with pfaf= " << (pfaf);

    MKL_free(A);MKL_free(temp1);MKL_free(x1);MKL_free(x2);
}