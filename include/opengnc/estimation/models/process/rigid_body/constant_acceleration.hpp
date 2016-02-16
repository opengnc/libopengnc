#ifndef OPENGNC_ESTIMATION_MODELS_PROCESS_RIGID_BODY_CONSTANT_ACCELERATION_HPP
#define OPENGNC_ESTIMATION_MODELS_PROCESS_RIGID_BODY_CONSTANT_ACCELERATION_HPP

#include <opengnc/common/math.hpp>

namespace opengnc {
namespace estimation {
namespace models {
namespace process {
namespace rigid_body {

template<typename state_policy, typename covariance_policy>
class constant_acceleration
{
public:
    enum { input_length = state_policy::state_vector_length, output_length = state_policy::state_vector_length };

    typedef typename state_policy::scalar_type input_scalar_type;
    typedef typename state_policy::scalar_type output_scalar_type;

    typedef Eigen::Matrix<input_scalar_type, input_length, 1> x_vec;
    typedef Eigen::Matrix<output_scalar_type, output_length, 1> y_vec;
    typedef Eigen::Matrix<output_scalar_type, output_length, output_length> y_mat;

    typedef typename covariance_policy::param_type param_type;

    constant_acceleration()
        : _timestep(0.02)
    {}

    void set_timestep(float value) { _timestep = value; }

    void set_params(const param_type& value) { _params = value; }

    y_vec operator() (const x_vec& x)
    {
		using namespace opengnc::common;
		using namespace Eigen;

        auto thetanb  = state_policy::thetanb(x);
        auto vBNb = state_policy::vBNb(x);
		auto omegaBNb = state_policy::omegaBNb(x);

		auto Rnb = math::rotationQuaternion(thetanb);
		auto Tnb = math::transform(thetanb);

        // Position dynamics
        auto drBNn = vBNb;
        drBNn  = Rnb*vBNb;

        //Attitude dynamics
        auto dthetanb = thetanb;
		dthetanb = Tnb*omegaBNb;

        //Velocity dynamics
        auto dvBNb = vBNb;
        dvBNb.setConstant(0);

		auto domegaBNb = omegaBNb;
        domegaBNb.setConstant(0);

        x_vec dx;
		state_policy::pack_state(dx, drBNn, dthetanb, dvBNb, domegaBNb);

        y_vec x_hat = x + _timestep*dx;
        state_policy::apply_constraints(x_hat);

        return x_hat;
    }

    y_mat uncertainty(const x_vec& x)
    {
        return _covariance_policy.apply(x, _timestep, _params);
    }

private:
    float _timestep;
    param_type _params;
    covariance_policy _covariance_policy;
};

}
}
}
}
}


#endif // OPENGNC_ESTIMATION_MODELS_PROCESS_RIGID_BODY_CONSTANT_ACCELERATION_HPP
