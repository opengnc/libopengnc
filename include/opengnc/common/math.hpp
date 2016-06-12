#ifndef OPENGNC_COMMON_MATH_HPP
#define OPENGNC_COMMON_MATH_HPP

#include "Eigen/Dense"
#include <complex>

namespace opengnc {
namespace common {

struct math
{
    template <typename scalar>
    static const Eigen::Matrix<scalar, 3, 3> Skew(const Eigen::Matrix<scalar, 3, 1>& vector)
    {
        Eigen::Matrix<scalar, 3, 3> s;
        s <<              0, -vector(2),  vector(1),
                vector(2),           0, -vector(0),
                -vector(1),  vector(0),           0;
        return s;
    }

    template <typename scalar>
    static const Eigen::Matrix<scalar, 3, 3> Skew(const Eigen::Matrix<scalar, 4, 1>& vector)
    {
        Eigen::Matrix<scalar, 3, 3> s;
        s <<              0, -vector(3),  vector(2),
                vector(3),           0, -vector(1),
                -vector(2),  vector(1),           0;
        return s;
    }

    template <typename scalar>
    static const Eigen::Matrix<scalar, 3, 3> rotationQuaternion(const Eigen::Matrix<scalar, 4, 1>& thetanb)
    {
        Eigen::Matrix<scalar, 3, 4> Rot_A;
        Eigen::Matrix<scalar, 3, 4> Rot_B;
        Eigen::Matrix<scalar, 3, 3> Rot;
        Eigen::Matrix<scalar, 3, 3> thetanb_eye;
        Eigen::Matrix<scalar, 3, 3> s = Skew(thetanb);

        thetanb_eye = thetanb[0] * Eigen::Matrix<scalar, 3, 3>::Identity();
        Rot_A.template block<3,1>(0,0) = -thetanb.template block<3,1>(1,0);
        Rot_A.template block<3,3>(0,1) = thetanb_eye + s;

        Rot_B.template block<3,1>(0,0) = -thetanb.template block<3,1>(1,0);
        Rot_B.template block<3,3>(0,1) = thetanb_eye - s;

        Rot = Rot_A * Rot_B.transpose();

        return Rot;
    }

    template <typename scalar>
    static const Eigen::Matrix<scalar, 4, 3> transform(const Eigen::Matrix<scalar, 4, 1>& thetanb)
    {
        Eigen::Matrix<scalar, 4, 3> T;
        Eigen::Matrix<scalar, 3, 3> s = Skew(thetanb);
        Eigen::Matrix<scalar, 3, 3> thetanb_eye = thetanb[0] * Eigen::Matrix<scalar, 3, 3>::Identity();

        T.template block<1,3>(0,0) = -(thetanb.template block<3,1>(1,0)).transpose();
        T.template block<3,3>(1,0) = thetanb_eye + s;
        T = 0.5*T;

        return T;
    }

    template <typename scalar>
    static const Eigen::Matrix<scalar, 4, 1> quaternionRotation(const Eigen::Matrix<scalar, 3, 3>& R)
    {
        Eigen::Matrix<scalar, 4, 4> K;
        K << R(0,0) - R(1,1) - R(2,2), R(1,0) + R(0,1), R(2,0) + R(0,2), R(1,2) - R(2,1),
                R(1,0) + R(0,1), R(1,1) - R(0,0) - R(2,2), R(2,1) + R(1,2), R(2,0) - R(0,2),
                R(2,0) + R(0,2), R(2,1) + R(1,2), R(2,2) - R(0,0) - R(1,1), R(0,1) - R(1,0),
                R(1,2) - R(2,1), R(2,0) - R(0,2), R(0,1) - R(1,0), R(0,0) + R(1,1) + R(2,2);
        K *= 1.0/3.0;

        Eigen::EigenSolver<Eigen::Matrix<scalar, 4, 4>> es(K);
        Eigen::Matrix<std::complex<scalar>, 4, 1> D = es.eigenvalues();
        Eigen::Matrix<std::complex<scalar>, 4, 4> V = es.eigenvectors();

        int idx;
        D.cwiseAbs().maxCoeff(&idx);

        Eigen::Matrix<scalar, 4, 1> q;
        q << -V(3, idx).real(),
                V(0, idx).real(),
                V(1, idx).real(),
                V(2, idx).real();
        return q;
    }

    template <typename scalar, int m, int n>
    static Eigen::Matrix<scalar, m*n, 1> vectorise(const Eigen::Matrix<scalar, m, n>& X)
    {
        Eigen::Matrix<scalar, m*n, 1> Y;

        for (int j = 0; j < n; ++j)
        {
            Y.template segment<m>(m*j) = X.template block<m,1>(0,j);
        }

        return Y;
    }

    template<typename scalar>
    static Eigen::Matrix<scalar, 3, 3> rotationEuler(const Eigen::Matrix<scalar, 3, 1>& v)
    {
        scalar sphi = std::sin(v[0]);
        scalar cphi = std::cos(v[0]);
        scalar ctheta = std::cos(v[1]);
        scalar stheta = std::sin(v[1]);
        scalar cpsi = std::cos(v[2]);
        scalar spsi = std::sin(v[2]);

        Eigen::Matrix<scalar, 3, 3> R;
        R << cpsi*ctheta, -spsi*cphi + cpsi*stheta*sphi,  spsi*sphi + cpsi*cphi*stheta,
                spsi*ctheta,  cpsi*cphi + spsi*stheta*sphi, -cpsi*sphi + spsi*cphi*stheta,
                -stheta,                   ctheta*sphi,                   ctheta*cphi;

        return R;
    }

    template<typename scalar>
    static Eigen::Matrix<scalar, 3, 1> eulerRotation(const Eigen::Matrix<scalar, 3, 3>& R)
    {
         Eigen::Matrix<scalar, 3, 1> euler;

        auto& phi = euler(0);
        auto& theta = euler(1);
        auto& psi = euler(2);

        theta = -asin(R(2,0));
        auto cos_theta = cos(theta);
        psi = atan2(R(1,0)/cos_theta, R(0,0)/cos_theta);
        phi = atan2(R(2,1)/cos_theta, R(2,2)/cos_theta);

        return euler;
    }

    template<typename scalar>
    static Eigen::Matrix<scalar, 6, 1> euler6Rotation(const Eigen::Matrix<scalar, 3, 3>& R)
    {
         Eigen::Matrix<scalar, 6, 1> euler6;

        auto& sin_phi = euler6(0);
        auto& cos_phi = euler6(1);
        auto& sin_theta = euler6(2);
        auto& cos_theta = euler6(3);
        auto& sin_psi = euler6(4);
        auto& cos_psi = euler6(5);

        auto euler = eulerRotation(R);

        sin_theta = sin(euler(1));
        cos_theta = cos(euler(1));
        sin_phi = sin(euler(0));
        cos_phi = cos(euler(0));
        sin_psi = sin(euler(2));
        cos_psi = cos(euler(2));

        return euler6;
    }
};

}
}

#endif // OPENGNC_COMMON_MATH_HPP

