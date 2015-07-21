#ifndef OPENGNC_ESTIMATION_MODELS_MEASUREMENT_GPS_HPP
#define OPENGNC_ESTIMATION_MODELS_MEASUREMENT_GPS_HPP

#include "opengnc/common/transforms/wgs84.hpp"

namespace opengnc {
namespace estimation {
namespace models {
namespace measurement {

template <typename state_policy>
class gps : state_policy
{
private:
	typedef Eigen::Matrix<typename state_policy::scalar_type, 3, 3> Matrix3s;
    bool _initialised;

protected:
	Eigen::Vector3d _rNOe;
	Matrix3s _Ren;

public:
    enum { input_length = state_policy::state_vector_length, output_length = 3 };
    typedef double output_scalar_type;
    typedef typename state_policy::scalar_type input_scalar_type;
    typedef Eigen::Vector3d YVector;
    typedef Eigen::Matrix<input_scalar_type, input_length, 1> XVector;

    gps()
        : _initialised(false)
    {}

    void initialise(Eigen::Vector3d& rNOe, Matrix3s& Ren)
    {
		_rNOe = rNOe;
		_Ren = Ren;
        _initialised = true;
    }

	bool isInitialised() const
    {
        return _initialised;
    }

    YVector operator() (const XVector& x) const
    {
		Eigen::Matrix<typename state_policy::scalar_type, 3, 1> rBNe = _Ren * state_policy::rBNn(x);
		Eigen::Matrix<scalar_type, 3, 1> rBOe;// = rBNe.cast<double>() + _rNOe;

        return opengnc::common::transforms::wgs84::rectangular_to_geodetic(rBOe);
    }

};

}
}
}
}

#endif // GPSMEASUREMENTMODEL_H
