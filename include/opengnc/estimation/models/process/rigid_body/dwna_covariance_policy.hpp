#ifndef OPENGNC_ESTIMATION_MODELS_PROCESS_RIGID_BODY_DWNA_COVARIANCE_POLICY_HPP
#define OPENGNC_ESTIMATION_MODELS_PROCESS_RIGID_BODY_DWNA_COVARIANCE_POLICY_HPP

#include <Eigen/Core>

namespace opengnc {
namespace estimation {
namespace models {
namespace process {
namespace rigid_body {

/* Based on the Discrete White Noise Acceleration Model from
 * the book 'Estimation with Applications to Tracking and Navigation',
 * page 272. Assumes state vector x is packed as {pos, vel, parameters}
 */
template<typename _state_policy>
struct dwna_covariance_policy
{
    typedef Eigen::Matrix<typename _state_policy::scalar_type, _state_policy::state_vector_length, 1> x_vec;
    typedef Eigen::Matrix<typename _state_policy::scalar_type, _state_policy::state_vector_length, _state_policy::state_vector_length> x_mat;
    typedef Eigen::Matrix<typename _state_policy::scalar_type, 6, 1> param_type;

    x_mat apply(const x_vec& x, float timestep, const param_type& Q_input_diag)
    {
        typedef Eigen::Matrix<typename _state_policy::scalar_type, 6, 6> Matrix6s;
        typedef Eigen::Matrix<typename _state_policy::scalar_type, _state_policy::state_vector_length, 6> Gamma_mat;

        float T = timestep;
        auto thetanb  = _state_policy::thetanb(x);
        auto Rnb = _state_policy::rotation(thetanb);
        auto Tnb = _state_policy::transform(thetanb);

        typedef Eigen::Matrix<typename _state_policy::scalar_type,
                Rnb.RowsAtCompileTime + Tnb.RowsAtCompileTime,
                Rnb.ColsAtCompileTime + Tnb.ColsAtCompileTime> J_mat;

        J_mat J = J_mat::Zero();
        J.block<Rnb.rows(),Rnb.cols()>(0,0) = Rnb;
        J.block<Tnb.rows(),Tnb.cols()>(Rnb.rows(),Rnb.cols()) = Tnb;


        Gamma_mat Gamma;
        if (_state_policy::state_vector_length == Eigen::Dynamic) Gamma.resize(x.size(), 6);

        Gamma.setZero();
        Gamma.block(0,0,J.rows(),J.cols()) = static_cast<typename _state_policy::scalar_type>(0.5)*T*T*J;
        Gamma.block(J.rows(),0,6,6) = T*Matrix6s::Identity();

        x_mat Q = Gamma*Q_input_diag.asDiagonal()*Gamma.transpose();

        return Q;
    }
};

}
}
}
}
}
#endif // OPENGNC_ESTIMATION_MODELS_PROCESS_RIGID_BODY_DWNA_COVARIANCE_POLICY_HPP
