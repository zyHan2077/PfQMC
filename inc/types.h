#ifndef TYPES_H
#define TYPES_H

#ifndef EIGEN_USE_MKL_ALL
#define EIGEN_USE_MKL_ALL
#endif

#include<float.h>
#include<complex>
#include<random>
#include<iostream>

#define MKL_Complex16 std::complex<double>

#include "mkl_types.h"
#include "mkl.h"
#include "mkl_vsl.h"
#include <Eigen/Dense>

typedef std::complex<double> DataType;
using namespace std::complex_literals;
const DataType zero{0.0, 0.0};
const DataType one{1.0, 0.0};
// const DataType im{0.0, 1.0};
typedef Eigen::MatrixXcd MatType;
typedef Eigen::VectorXcd cVecType;
typedef Eigen::VectorXd dVecType;
typedef Eigen::VectorXi iVecType;
#define ind(i, j, N) (((j)*N) + i)
const double thresholdDBL = 1.0e-300;

inline DataType logDet(const MatType &H) {
    DataType ld = 0.0;
    Eigen::PartialPivLU<MatType> lu(H);
    auto& LU = lu.matrixLU();
    DataType c = lu.permutationP().determinant(); // -1 or 1
    // std::cout << "c=" << c << "\n";
    for (unsigned i = 0; i < LU.rows(); ++i) {
      const DataType& lii = LU(i,i);
    //   std::cout << lii << " q\n";
      ld += log(lii);
    }
    ld += log(c);
    return ld;
}

class rdGenerator {

private:
    std::uniform_int_distribution<int> rdDistZ2;
    std::uniform_real_distribution<double> rdDistUniform01;
    std::mt19937 rdEng;
public:
    rdGenerator(int seed=114514) {
        // random generator, Z_2 auxillary field
        rdDistZ2 = std::uniform_int_distribution<int>(0, 1);
        rdDistUniform01 = std::uniform_real_distribution<double>(0.0, 1.0);
        rdEng.seed(seed);
    }

    inline int rdZ2() {
        return 2*rdDistZ2(rdEng) - 1;
    }

    inline int rdUniform01() {
        return rdDistUniform01(rdEng);
    }

};

#endif