#ifndef OPENGNC_ESTIMATION_UNSCENTED_KALMAN_FILTER_HPP
#define OPENGNC_ESTIMATION_UNSCENTED_KALMAN_FILTER_HPP

#include <cstdint>

#include <opengnc/estimation/time_update.hpp>
#include <opengnc/estimation/measurement_update.hpp>
#include <opengnc/estimation/unscented_transform.hpp>

namespace opengnc {
namespace estimation {

template<typename state_density_type,
         typename state_policy,
         typename process_model_type,
         typename container>
class unscented_kalman_filter : public container
{
    using process_propagater_type = opengnc::estimation::unscented_transform<process_model_type, state_policy>;
    using time_update_type = opengnc::estimation::time_update<state_density_type, process_model_type, process_propagater_type, state_policy>;

public:
    unscented_kalman_filter()
        : last_updated_ms(0)
        , process_update(process_model)
    { }

    void init(const state_density_type& state, uint64_t ms)
    {
        state_density = state;
        set_last_updated_ms(ms);
    }

    const state_density_type& get_state_density() const { return state_density; }

    uint64_t get_last_updated_ms() const { return last_updated_ms; }

    void set_last_updated_ms(uint64_t value) { last_updated_ms = value; }

    void time_update(uint64_t time_ms)
    {
        double time_step_sec = (time_ms - last_updated_ms)/1000.0;
        process_model.set_timestep(time_step_sec);
        process_update(state_density);
        last_updated_ms = time_ms;
    }

    process_model_type& get_process_model()
    {
        return process_model;
    }

    template<typename measurement_traits>
    void measurement_update(typename measurement_traits::density_type& measurement_density)
    {
        auto& model = container::template get_measurement_model<measurement_traits>();
        auto& update = container::template get_measurement_update<measurement_traits, state_density_type>();
        update(state_density, measurement_density);
    }


private:
    uint64_t last_updated_ms;
    state_density_type state_density;
    process_model_type process_model;
    time_update_type process_update;

};

}
}
#endif // OPENGNC_ESTIMATION_UNSCENTED_KALMAN_FILTER_HPP
