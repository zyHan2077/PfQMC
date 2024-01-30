#include"src/skewMatUtils.h"
// #include <Eigen/Dense>
 
// using Eigen::MatrixXd;

int main() {

    // MatrixXd m(2,2);
    // m(0,0) = 3;
    // m(1,0) = 2.5;
    // m(0,1) = -1;
    // m(1,1) = m(1,0) + m(0,1);
    // std::cout << m << std::endl;

    int N = 10;
    dtype *A, *temp1, *temp2;
    A = (dtype *) mkl_malloc((4*N*N)*sizeof(dtype), 64);
    temp1 = (dtype *) mkl_malloc((4*N)*sizeof(dtype), 64);
    temp2 = (dtype *) mkl_malloc((2*N)*sizeof(dtype), 64);
    // y = (dtype *) mkl_malloc((2*N)*sizeof(dtype), 64);

    // initialize
    int count = 1;
    for(int i=0; i<2*N; i++) {
        for (int j=i+1; j<2*N; j++) {
            A[i*2*N + j] = dtype{(double)count, -0.0};
            A[j*2*N + i] = dtype{-(double)count, 0.0};
            count++;
        }
        temp1[i] = count;
    }
    
    // printVec(2*N, temp1);
    // const MKL_INT inc=1;
    // // const dtype two = 2;
    // MKL_INT len = 2*N;
    // const char mklNoTrans = 'N';
    // zgemv(&mklNoTrans, &len, &len, &one, A, &len, temp1, &inc, &zero, temp2, &inc);
    // printVec(2*N, temp2);

    // std::cout << "smallest double: " << thresholdDBL << "\n";

    printMat(2*N, A);

    dtype pfaf;
    SkewMatHouseholder_PureMKL(N, A, temp1, &pfaf);
    
    printMat(2*N, A);

    std::cout << "pfaffian= " << pfaf << "\n";

    MKL_free(A);MKL_free(temp1);MKL_free(temp2);

    return 0;
}