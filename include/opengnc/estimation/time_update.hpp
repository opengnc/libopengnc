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
        propagater.init(density.mean(), density.covariance());
        propagater.propagate(model, density.mean(), density.covariance());
        constraint_policy::apply(density.mean());
        density.covariance() = density.covariance() + model.uncertainty(density.mean());
    }

protected:
    process_model_type model;
    propagater_type propagater;
};

}
}

#endif // OPENGNC_ESTIMATION_TIME_UPDATE_HPP
