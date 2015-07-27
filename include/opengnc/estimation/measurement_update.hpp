#ifndef OPENGNC_ESTIMATION_MEASUREMENT_UPDATE_HPP
#define OPENGNC_ESTIMATION_MEASUREMENT_UPDATE_HPP

#include <Eigen/Core>
#include <Eigen/Cholesky>

namespace opengnc {
namespace estimation {

template<typename density_type,
         typename measurement_model_type,
         typename measurement_service_type,
         typename propagater_type,
         typename exclusion_policy,
         typename constraint_policy>
class measurement_update
        : exclusion_policy
        , constraint_policy
{
    enum { nx = measurement_model_type::input_length, ny = measurement_model_type::output_length };

    typedef typename measurement_model_type::input_scalar_type input_scalar_type;
    typedef typename measurement_model_type::output_scalar_type output_scalar_type;

    typedef Eigen::Matrix<output_scalar_type, ny, 1> y_vec;
    typedef Eigen::Matrix<output_scalar_type, ny, ny> y_mat;
    typedef Eigen::Matrix<output_scalar_type, nx, ny> xy_mat;

public:
    measurement_update(measurement_service_type& measurement_service, measurement_model_type& model)
        : _measurement_service(measurement_service)
        , _model(model)
    { }

    void operator () (density_type& density)
    {
        const auto& x = density.mean();
        const auto& Px = density.covariance();

        y_vec y_hat;
        y_mat Py_hat;
        xy_mat Pxy_hat;

        _propagater.init(x, Px);
        _propagater.propagate(_model, y_hat, Py_hat, Pxy_hat);

        y_vec y = _measurement_service.get_measurement();
        y_mat R = _measurement_service.get_measurement_uncertainty();

        if(y.rows() > 0 && y_hat.rows() > 0)
        {
            //Add measurement noise
            Py_hat +=  R;

            //Apply Exclusions
            exclusion_policy::apply(y, y_hat, Py_hat, Pxy_hat);
        }

        if(y.rows() > 0 && y_hat.rows() > 0)
        {
            //Conditioning: K = Pxy/Py
            xy_mat K = Py_hat.ldlt().solve(Pxy_hat.transpose()).transpose();

            density.mean() = x + (K*(y - y_hat)).template cast<typename density_type::scalar_type>();
            density.covariance() = Px - (K*Py_hat*K.transpose()).template cast<typename density_type::scalar_type>();

            //Apply Constaints
            constraint_policy::apply(density.mean());
        }
    }

private:
    measurement_service_type& _measurement_service;
    measurement_model_type& _model;
    propagater_type _propagater;
};

}
}

#endif // OPENGNC_ESTIMATION_MEASUREMENT_UPDATE_HPP
