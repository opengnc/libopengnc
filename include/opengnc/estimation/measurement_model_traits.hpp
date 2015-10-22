#ifndef OPENGNC_ESTIMATION_MEASUREMENT_MODEL_TRAITS_HPP
#define OPENGNC_ESTIMATION_MEASUREMENT_MODEL_TRAITS_HPP

#include <Eigen/Core>

namespace opengnc {
namespace estimation {

template<typename measurement_model_type>
struct measurement_model_traits
{
	enum { nx = measurement_model_type::input_length, ny = measurement_model_type::output_length };

	typedef typename measurement_model_type::input_scalar_type input_scalar_type;
	typedef typename measurement_model_type::output_scalar_type output_scalar_type;

	typedef Eigen::Matrix<output_scalar_type, ny, 1> y_vec;
	typedef Eigen::Matrix<output_scalar_type, ny, ny> y_mat;
	typedef Eigen::Matrix<output_scalar_type, nx, ny> xy_mat;
};

}
}

#endif // OPENGNC_ESTIMATION_MEASUREMENT_MODEL_TRAITS_HPP

