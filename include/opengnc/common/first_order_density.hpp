#ifndef OPENGNC_COMMON_FIRST_ORDER_DENSITY_HPP
#define OPENGNC_COMMON_FIRST_ORDER_DENSITY_HPP

#include <Eigen/Core>

namespace opengnc {
namespace common {

template<typename scalar_t, int nx=Eigen::Dynamic>
class first_order_density {
public:

    enum { size = nx };

    typedef scalar_t scalar_type;
    typedef Eigen::Matrix<double, size, 1> vec_type;
    typedef Eigen::Matrix<double, size, size> mat_type;

    first_order_density() { }

    first_order_density(vec_type mean, mat_type covariance)
        : _mean(mean)
        , _covariance(covariance)
    { }

    const vec_type& mean() const { return _mean; }
    vec_type& mean() { return _mean; }

    const mat_type& convariance() const { return _covariance; }
    mat_type& convariance() { return _covariance; }

private:
    vec_type _mean;
    mat_type _covariance;
};

}
}

#endif // OPENGNC_COMMON_FIRST_ORDER_DENSITY_HPP
