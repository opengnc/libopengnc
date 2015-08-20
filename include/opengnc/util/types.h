#ifndef NUGNC_TYPES_H
#define NUGNC_TYPES_H
#include <Eigen/Core>
#include <limits.h>

namespace GNC
{
#ifdef __USE_SINGLE_PRECISION__
typedef float Scalar;
#else
typedef double Scalar;
#endif

typedef Eigen::Matrix<Scalar,2,1> Vector2s;
typedef Eigen::Matrix<Scalar,2,2> Matrix2s;
}


#endif // TYPES

