#ifndef OPENGNC_ESTIMATION_TIME_UPDATE_HPP
#define OPENGNC_ESTIMATION_TIME_UPDATE_HPP

namespace opengnc {
namespace estimation {

template<typename density_type,
         typename process_model_type,
         typename propagater_type,
         typename constraint_policy>
class time_update : constraint_policy
{
public:
    time_update() { }

    void operator () (density_type& density)
    {
        _propagater.init(density.mean(), density.covariance());
        _propagater.propagate(_model, density.mean(), density.covariance());
        constraint_policy::apply(density.mean());
        density.covariance() = density.covariance() + _model.uncertainty(density.mean());
    }

protected:
    process_model_type _model;
    propagater_type _propagater;
};

}
}

#endif // OPENGNC_ESTIMATION_TIME_UPDATE_HPP
