#ifndef TYPES_H
#define TYPES_H

#ifndef EIGEN_USE_MKL_ALL
#define EIGEN_USE_MKL_ALL
#endif

#include<float.h>
#include<complex>

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
typedef Eigen::VectorXi iVecType;
#define ind(i, j, N) (((j)*N) + i)
const double thresholdDBL = 1.0e-300;

// template <typename T>
// struct types;

// template <>
// struct types<double>
// {
// 	using MatrixType = Eigen::MatrixXd;
// 	using DataType = double;
// 	using VectorType = Eigen::VectorXd;
// };

#endif