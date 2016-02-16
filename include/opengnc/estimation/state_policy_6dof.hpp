#ifndef STATE_POLICY_6DOF_HPP
#define STATE_POLICY_6DOF_HPP

#include "Eigen/Dense"
#include <complex>

namespace opengnc {
namespace estimation {

template <typename scalar, int x_vector_length=13>
class state_policy_6dof
{
protected:
	typedef Eigen::Matrix<scalar, 3, 1> Vector3s;
	typedef Eigen::Matrix<scalar, 4, 1> Vector4s;
	typedef Eigen::Matrix<std::complex<scalar>, 4, 1> Vector4cd;
	typedef Eigen::Matrix<scalar, 3, 3> Matrix3s;
	typedef Eigen::Matrix<scalar, 4, 4> Matrix4s;
	typedef Eigen::Matrix<std::complex<scalar>, 4, 4> Matrix4cd;
	typedef Eigen::Matrix<scalar, 3, 4> Matrix3x4s;
	typedef Eigen::Matrix<scalar, 4, 3> Matrix4x3s;

public:
	enum { state_vector_length = 13 };
        enum { parameter_vector_length = 0 };
	typedef scalar scalar_type;
        typedef Eigen::Matrix<scalar_type, x_vector_length, 1> XVector;
        typedef Eigen::Matrix<scalar_type, x_vector_length, x_vector_length> PMatrix;

	static const Vector3s rBNn(const XVector& x)
	{
		return x.template segment<3>(0);
	}

	static const Vector4s thetanb(const XVector& x)
	{
		return x.template segment<4>(3);
	}

	static const Vector3s vBNb(const XVector& x)
	{
		return x.template segment<3>(7);
	}

	static const Vector3s omegaBNb(const XVector& x)
	{
		return x.template segment<3>(10);
	}

	static void pack_state(XVector& x,
						   const Vector3s& rBNn,
						   const Vector4s& thetanb,
						   const Vector3s& vBNb,
						   const Vector3s& omegaBNb)
	{
		x.template segment<3>(0) = rBNn;
		x.template segment<4>(3) = thetanb;
		x.template segment<3>(7) = vBNb;
		x.template segment<3>(10) = omegaBNb;
	}

	static void apply_constraints(XVector&)
	{
		// TODO: Something
	}
};

}
}

#endif // STATE_POLICY_6DOF_HPP

