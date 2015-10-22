#ifndef OPENGNC_ESTIMATION_MEASUREMENT_UPDATE_HPP
#define OPENGNC_ESTIMATION_MEASUREMENT_UPDATE_HPP

#include <Eigen/Core>
#include <Eigen/Cholesky>
#include "measurement_condition.hpp"
#include "measurement_model_traits.hpp"

namespace opengnc {
namespace estimation {

template<typename state_density_type,
         typename measurement_density_type,
         typename measurement_model_type,
         typename propagater_type,
         typename state_constraint_policy,
         typename measurement_exclusion_policy=void>
class measurement_update : public measurement_condition<measurement_model_traits<measurement_model_type>, measurement_exclusion_policy>
{
    enum { nx = measurement_model_type::input_length, ny = measurement_model_type::output_length };

	using condition_measurements_type = measurement_condition<measurement_model_traits<measurement_model_type>, measurement_exclusion_policy>;
	using traits = measurement_model_traits<measurement_model_type>;

	using input_scalar_type = typename traits::input_scalar_type;
	using output_scalar_type = typename traits::output_scalar_type;

	using y_vec = typename traits::y_vec;
	using y_mat = typename traits::y_mat;
	using xy_mat = typename traits::xy_mat;

public:
    measurement_update(measurement_model_type& model)
		: condition_measurements_type(_y_hat, _y_hat_used, _y_used)
		, _model(model)
    { }

    void operator () (state_density_type& state_density, const measurement_density_type& measurement_density)
    {
        const auto& x = state_density.mean();
        const auto& Px = state_density.covariance();

        y_mat Py_hat;
        xy_mat Pxy_hat;

        _propagater.init(x, Px);
        _propagater.propagate(_model, _y_hat, Py_hat, Pxy_hat);
        _y_hat_used = _y_hat;

        _y_used = measurement_density.mean();
        const y_mat& R = measurement_density.covariance();

		condition_measurements_type::apply(R, Py_hat, Pxy_hat);

        if(_y_used.rows() > 0 && _y_hat_used.rows() > 0)
        {
            //Conditioning: K = Pxy/Py
            xy_mat K = Py_hat.ldlt().solve(Pxy_hat.transpose()).transpose();

            state_density.mean() = x + (K*(_y_used - _y_hat_used)).template cast<typename state_density_type::scalar_type>();
            state_density.covariance() = Px - (K*Py_hat*K.transpose()).template cast<typename state_density_type::scalar_type>();

            //Apply Constaints
            state_constraint_policy::apply_constraints(state_density.mean());
        }
    }

    const y_vec& predicted_measurements() const { return _y_hat; }

    const y_vec& used_predicted_measurements() const { return _y_hat_used; }

    const y_vec& used_actual_measurements() const { return _y_used; }

private:
    measurement_model_type& _model;
    propagater_type _propagater;
    y_vec _y_hat;
    y_vec _y_hat_used;
    y_vec _y_used;
};

}
}

#endif // OPENGNC_ESTIMATION_MEASUREMENT_UPDATE_HPP
