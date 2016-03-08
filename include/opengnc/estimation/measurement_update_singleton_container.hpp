#ifndef MEASUREMENT_UPDATE_SINGLETON_CONTAINER_HPP
#define MEASUREMENT_UPDATE_SINGLETON_CONTAINER_HPP

#include <opengnc/estimation/measurement_update.hpp>

namespace opengnc {
namespace estimation {

struct measurement_update_singleton_container
{
    template<typename measurement_traits>
    typename measurement_traits::model_type& get_measurement_model()
    {
        static typename measurement_traits::model_type instance;
        return instance;
    }

    template<typename measurement_traits, typename state_density_type>
    typename opengnc::estimation::measurement_update
                    < state_density_type
                    , typename measurement_traits::density_type
                    , typename measurement_traits::model_type
                    , typename measurement_traits::propagator_type
                    , typename measurement_traits::constraint_policy
                    , typename measurement_traits::exclusion_policy>& get_measurement_update()
    {
        static typename opengnc::estimation::measurement_update
                < state_density_type
                , typename measurement_traits::density_type
                , typename measurement_traits::model_type
                , typename measurement_traits::propagator_type
                , typename measurement_traits::constraint_policy
                , typename measurement_traits::exclusion_policy> instance(get_measurement_model<measurement_traits>());

        return instance;
    }
};

}
}

#endif // MEASUREMENT_UPDATE_SINGLETON_CONTAINER_HPP
