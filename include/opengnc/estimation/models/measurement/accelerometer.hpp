#ifndef OPENGNC_ESTIMATION_MODELS_MEASUREMENT_ACCELEROMETER_HPP
#define OPENGNC_ESTIMATION_MODELS_MEASUREMENT_ACCELEROMETER_HPP

#include <opengnc/common/math.hpp>

namespace opengnc {
namespace estimation {
namespace models {
namespace measurement {

template <typename state_policy>
class accelerometer
{
private:
    using Vector3s = Eigen::Matrix<typename state_policy::scalar_type, 3, 1>;
    using Matrix3s = Eigen::Matrix<typename state_policy::scalar_type, 3, 3>;
    const Vector3s _gn;
    Matrix3s _Rib;
    Vector3s _rIBb;                                         //!< Distance from Body (B) to IMU (I) in body coordinates {b}
    bool _initialised;

public:
    enum { input_length = state_policy::state_vector_length, output_length = 3 };
    using output_scalar_type = float;
    using input_scalar_type = typename state_policy::scalar_type;
    using YVector = Vector3s;
    using XVector = Eigen::Matrix<input_scalar_type, input_length, 1>;

    accelerometer()
        : _gn(0, 0, -9.8)
		, _Rib(Matrix3s::Identity())
		, _rIBb(Vector3s::Zero())
        , _initialised(false)
    {}

    void initialise(const Matrix3s& Rib, const Vector3s& rIBb)
    {
        _Rib = Rib;
        _rIBb = rIBb;
        _initialised = true;
    }

    bool initialised() const { return _initialised; }

    YVector operator() (const XVector& x) const
    {
		using namespace opengnc::common;

		const auto Rnb = math::rotationQuaternion(state_policy::thetanb(x));
        const auto omegaBNb = state_policy::omegaBNb(x);

		YVector y = _Rib*(math::Skew(omegaBNb)*_rIBb - Rnb.transpose()*_gn);

        return y;
    }

};

}
}
}
}

#endif // OPENGNC_ESTIMATION_MODELS_MEASUREMENT_ACCELEROMETER_HPP
