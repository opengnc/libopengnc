#ifndef NO_EXCLUSION_POLICY_HPP
#define NO_EXCLUSION_POLICY_HPP

#include <Eigen/Core>

namespace opengnc {
namespace estimation {

template<typename measurement_model_type>
class no_exclusion_policy
{
private:
	enum {	ny = measurement_model_type::output_length,
			nx = measurement_model_type::input_length };
	using scalar = typename measurement_model_type::output_scalar_type;

	typedef Eigen::Matrix<scalar, ny, 1> y_vec;
	typedef Eigen::Matrix<scalar, ny, ny> y_mat;
	typedef Eigen::Matrix<scalar, nx, ny> xy_mat;

public:
	static void apply_exclusions(const y_vec&, const y_vec&, const y_mat&, const xy_mat&) { }
};

}
}

#endif // NO_EXCLUSION_POLICY_HPP

