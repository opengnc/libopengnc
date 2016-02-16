#ifndef OPENGNC_ESTIMATION_STATE_POLICY_6DOF_WITH_GYRO_HPP
#define OPENGNC_ESTIMATION_STATE_POLICY_6DOF_WITH_GYRO_HPP

#include "state_policy_6dof.hpp"

namespace opengnc {
namespace estimation {

template <typename scalar>
class state_policy_6dof_with_gyro : public state_policy_6dof<scalar, 16>
{
private:
    using parent = state_policy_6dof<scalar, 16>;

public:
    enum { state_vector_length = 16 };
    enum { parameter_vector_length = 3 };

    static const typename parent::Vector3s gbBNb(const typename parent::XVector& x)
    {
        return x.template segment<3>(parent::state_vector_length);
    }

    static void pack_parameters(typename parent::XVector& x,
                                const typename parent::Vector3s& gbBNb)
    {
        x.template segment<3>(parent::state_vector_length) = gbBNb;
    }
};

}
}

#endif // OPENGNC_ESTIMATION_STATE_POLICY_6DOF_WITH_GYRO_HPP
