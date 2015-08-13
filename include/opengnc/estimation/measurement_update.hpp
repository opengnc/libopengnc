#ifndef OPENGNC_ESTIMATION_MEASUREMENT_UPDATE_HPP
#define OPENGNC_ESTIMATION_MEASUREMENT_UPDATE_HPP

#include <Eigen/Core>
#include <Eigen/Cholesky>

namespace opengnc {
namespace estimation {

template<typename state_density_type,
         typename measurement_density_type,
         typename measurement_model_type,
         typename propagater_type,
         typename state_constraint_policy,
         typename measurement_exclusion_policy>
class measurement_update
{
    enum { nx = measurement_model_type::input_length, ny = measurement_model_type::output_length };

    typedef typename measurement_model_type::input_scalar_type input_scalar_type;
    typedef typename measurement_model_type::output_scalar_type output_scalar_type;

    typedef Eigen::Matrix<output_scalar_type, ny, 1> y_vec;
    typedef Eigen::Matrix<output_scalar_type, ny, ny> y_mat;
    typedef Eigen::Matrix<output_scalar_type, nx, ny> xy_mat;

public:
    measurement_update(measurement_model_type& model)
        : _model(model)
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

        if(_y_used.rows() > 0 && _y_hat.rows() > 0)
        {
            //Add measurement noise
            Py_hat +=  R;
            measurement_exclusion_policy::apply_exclusions(_y_used, _y_hat_used, Py_hat, Pxy_hat);
        }

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
