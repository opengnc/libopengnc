#ifndef STATE_POLICY_6DOF_HPP
#define STATE_POLICY_6DOF_HPP

#include "Eigen/Dense"
#include <complex>

namespace opengnc {
namespace estimation {

template <typename scalar>
class state_policy_6dof
{
private:
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
	typedef scalar scalar_type;
	typedef Eigen::Matrix<scalar_type, state_vector_length, 1> XVector;

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

	// TODO: Move all following functions to a separate utils class
	static const Matrix3s Skew(const Vector4s& vector)
	{
		Matrix3s s;
		s <<              0, -vector[3],  vector[2],
				 vector[3],           0, -vector[1],
				-vector[2],  vector[1],           0;
		return s;
	}

	static const Matrix3s rotation(const Vector4s& thetanb)
	{
		Matrix3x4s Rot_A;
		Matrix3x4s Rot_B;
		Matrix3s Rot;
		Matrix3s thetanb_eye;
		Matrix3s s = Skew(thetanb);

		thetanb_eye = thetanb[0] * Matrix3s::Identity();
		Rot_A.template block<3,1>(0,0) = -thetanb.template block<3,1>(1,0);
		Rot_A.template block<3,3>(0,1) = thetanb_eye + s;

		Rot_B.template block<3,1>(0,0) = -thetanb.template block<3,1>(1,0);
		Rot_B.template block<3,3>(0,1) = thetanb_eye - s;

		Rot = Rot_A * Rot_B.transpose();

		return Rot;
	}



	static const Matrix4x3s transform(const Vector4s& thetanb)
	{
		Matrix4x3s T;
		Matrix3s s = Skew(thetanb);
		Matrix3s thetanb_eye = thetanb[0] * Matrix3s::Identity();

		T.template block<1,3>(0,0) = -(thetanb.template block<3,1>(1,0)).transpose();
		T.template block<3,3>(1,0) = thetanb_eye + s;
		T = 0.5*T;

		return T;
	}

	static const Vector4s quaternionRotation(const Matrix3s& R)
	{
		Matrix4s K;
		K << R(0,0) - R(1,1) - R(2,2), R(1,0) + R(0,1), R(2,0) + R(0,2), R(1,2) - R(2,1),
				R(1,0) + R(0,1), R(1,1) - R(0,0) - R(2,2), R(2,1) + R(1,2), R(2,0) - R(0,2),
				R(2,0) + R(0,2), R(2,1) + R(1,2), R(2,2) - R(0,0) - R(1,1), R(0,1) - R(1,0),
				R(1,2) - R(2,1), R(2,0) - R(0,2), R(0,1) - R(1,0), R(0,0) + R(1,1) + R(2,2);
		K *= 1.0/3.0;

		Eigen::EigenSolver<Matrix4s> es(K);
		Vector4cd D = es.eigenvalues();
		Matrix4cd V = es.eigenvectors();

		int idx;
		D.cwiseAbs().maxCoeff(&idx);

		Vector4s q;
		q << -V(3, idx).real(),
				V(0, idx).real(),
				V(1, idx).real(),
				V(2, idx).real();
		return q;
	}
};

}
}

#endif // STATE_POLICY_6DOF_HPP

