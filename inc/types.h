#ifndef TYPES_H
#define TYPES_H

#ifndef EIGEN_USE_MKL_ALL
#define EIGEN_USE_MKL_ALL
#endif

#include<complex>
#include<float.h>

#define MKL_Complex16 std::complex<double>

#include "mkl_types.h"
#include "mkl.h"
#include "mkl_vsl.h"
#include <Eigen/Dense>

typedef std::complex<double> dtype;
const dtype zero{0.0, 0.0};
const dtype one{1.0, 0.0};
const dtype im{0.0, 1.0};
typedef Eigen::MatrixXcd cMat;
typedef Eigen::VectorXcd cVec;
#define ind(i, j, N) (((i)*N) + j)
const double thresholdDBL = 1.0e-300;

#endif