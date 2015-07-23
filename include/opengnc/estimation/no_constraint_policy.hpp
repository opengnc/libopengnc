#ifndef OPENGNC_ESTIMATION_NO_CONSTRAINT_POLICY_HPP
#define OPENGNC_ESTIMATION_NO_CONSTRAINT_POLICY_HPP

#include <Eigen/Core>

namespace opengnc {
namespace estimation {

template<typename scalar_type, int size>
struct no_constraint_policy
{
    typedef Eigen::Matrix<scalar_type, size, 1> x_vec;
    void apply(x_vec&) { }
};

}
}
#endif // OPENGNC_ESTIMATION_NO_CONSTRAINT_POLICY_HPP
