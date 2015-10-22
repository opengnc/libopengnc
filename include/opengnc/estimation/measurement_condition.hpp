#ifndef OPENGNC_ESTIMATION_MEASUREMENT_CONDITION_HPP
#define OPENGNC_ESTIMATION_MEASUREMENT_CONDITION_HPP

#include "measurement_model_traits.hpp"

namespace opengnc {
namespace estimation {

template <typename traits, typename exclusion_policy>
class measurement_condition
{
private:
	typename traits::y_vec& _y_hat;
	typename traits::y_vec& _y_hat_used;
	typename traits::y_vec& _y_used;
	typename exclusion_policy::param_type _exclusion_param;

protected:
	void apply (const typename traits::y_mat& R,
				typename traits::y_mat& Py_hat,
				typename traits::xy_mat& Pxy_hat)
	{
		if(_y_used.rows() > 0 && _y_hat.rows() > 0)
		{
			//Add measurement noise
			Py_hat +=  R;
			exclusion_policy::apply_exclusions(_y_used, _y_hat_used, Py_hat, Pxy_hat, _exclusion_param);
		}
	}

public:
	measurement_condition(typename traits::y_vec& y_hat,
						  typename traits::y_vec& y_hat_used,
						  typename traits::y_vec& y_used)
		: _y_hat(y_hat)
		, _y_hat_used(y_hat_used)
		, _y_used(y_used)
	{}

	void set_exclusion_parameters(typename exclusion_policy::param_type param) { _exclusion_param = param; }
};

template <typename traits>
class measurement_condition<traits, void>
{
private:
	typename traits::y_vec& _y_hat;
	typename traits::y_vec& _y_hat_used;
	typename traits::y_vec& _y_used;

protected:
	void apply (const typename traits::y_mat& R,
				typename traits::y_mat& Py_hat,
				typename traits::xy_mat& Pxy_hat)
	{
		if(_y_used.rows() > 0 && _y_hat.rows() > 0)
		{
			//Add measurement noise
			Py_hat +=  R;
		}
	}

public:
	measurement_condition(typename traits::y_vec& y_hat,
						  typename traits::y_vec& y_hat_used,
						  typename traits::y_vec& y_used)
		: _y_hat(y_hat)
		, _y_hat_used(y_hat_used)
		, _y_used(y_used)
	{}
};

}
}

#endif // OPENGNC_ESTIMATION_MEASUREMENT_CONDITION_HPP

