#ifndef OPENGNC_ESTIMATION_STATE_POLICY_6DOF_WITH_GYRO_HPP
#define OPENGNC_ESTIMATION_STATE_POLICY_6DOF_WITH_GYRO_HPP

#include "state_policy_6dof.hpp"

namespace opengnc {
namespace estimation {

template <typename scalar>
class state_policy_6dof_with_gyro : public state_policy_6dof<scalar, 16>
{
private:
    using base = state_policy_6dof<scalar, 16>;

public:
    enum { state_vector_length = 16 };
    enum { parameter_vector_length = 3 };

    static const typename base::Vector3s gbBNi(const typename base::XVector& x)
    {
        return x.template segment<3>(base::state_vector_length);
    }

    static void pack_parameters(typename base::XVector& x,
                                const typename base::Vector3s& gbBNi)
    {
        x.template segment<3>(base::state_vector_length) = gbBNi;
    }
};

}
}

#endif // OPENGNC_ESTIMATION_STATE_POLICY_6DOF_WITH_GYRO_HPP
