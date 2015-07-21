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
template<typename state_policy>
struct dwna_covariance_policy : state_policy
{
    typedef Eigen::Matrix<typename state_policy::scalar_type, state_policy::state_vector_length, 1> x_vec;
    typedef Eigen::Matrix<typename state_policy::scalar_type, state_policy::state_vector_length, state_policy::state_vector_length> x_mat;
    typedef Eigen::Matrix<typename state_policy::scalar_type, 6, 1> param_type;

    x_mat apply(const x_vec& x, float timestep, const param_type& Q_input_diag)
    {
        typedef Eigen::Matrix<typename state_policy::scalar_type, 6, 6> Matrix6s;
        typedef Eigen::Matrix<state_policy::scalar_type, state_policy::state_vector_length, 6> Gamma_mat;

        float T = timestep;
        auto thetanb  = state_policy::thetanb(x);
        auto Rnb = state_policy::rotation(thetanb);
        auto Tnb = state_policy::transform(thetanb);

        Eigen::MatrixXd top_right = MatrixXd::Zero(Rnb.rows(), Tnb.cols());
        Eigen::MatrixXd bottom_left = MatrixXd::Zero(Tnb.rows(), Rnb.cols());

        J << Rnb        , top_right,
             bottom_left, Tnb;

        Gamma_mat Gamma = Gamma_mat::Zero();

        Gamma.block(0,0,J.rows(),J.cols()) = 0.5*T*T*J;
        Gamma.block<6,6>(J.rows(),0) = T*Matrix6s::Identity();

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
