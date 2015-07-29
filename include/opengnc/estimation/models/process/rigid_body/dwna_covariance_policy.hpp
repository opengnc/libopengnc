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
struct dwna_covariance_policy : private _state_policy
{
    typedef Eigen::Matrix<typename _state_policy::scalar_type, _state_policy::state_vector_length, 1> x_vec;
    typedef Eigen::Matrix<typename _state_policy::scalar_type, _state_policy::state_vector_length, _state_policy::state_vector_length> x_mat;
    typedef Eigen::Matrix<typename _state_policy::scalar_type, 6, 1> param_type;

    x_mat apply(const x_vec& x, float timestep, const param_type& Q_input_diag)
    {
        typedef Eigen::Matrix<typename _state_policy::scalar_type, 6, 6> Matrix6s;
        typedef Eigen::Matrix<typename _state_policy::scalar_type, _state_policy::state_vector_length, 6> Gamma_mat;
        typedef Eigen::Matrix<typename _state_policy::scalar_type, Eigen::Dynamic, Eigen::Dynamic> MatrixXs;

        float T = timestep;
        auto thetanb  = _state_policy::thetanb(x);
        auto Rnb = _state_policy::rotation(thetanb);
        auto Tnb = _state_policy::transform(thetanb);

        MatrixXs J = MatrixXs::Zero(Rnb.rows()+Tnb.rows(),Rnb.cols()+Tnb.cols());
        J.block(0,0,Rnb.rows(),Rnb.cols()) = Rnb;
        J.block(Rnb.rows(),Rnb.cols(),Tnb.rows(),Tnb.cols()) = Tnb;

        Gamma_mat Gamma;
        if (_state_policy::state_vector_length == Eigen::Dynamic) Gamma.resize(x.size(), 6);

        Gamma.setZero();
        Gamma.block(0,0,J.rows(),J.cols()) = 0.5*T*T*J;
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
