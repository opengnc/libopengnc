#ifndef OPENGNC_ESTIMATION_MODELS_MEASUREMENT_GYROSCOPE_HPP
#define OPENGNC_ESTIMATION_MODELS_MEASUREMENT_GYROSCOPE_HPP

#include <Eigen/Core>

namespace opengnc {
namespace estimation {
namespace models {
namespace measurement {

template <typename state_policy>
class gyroscope
{
private:
    using Vector3s = Eigen::Matrix<typename state_policy::scalar_type, 3, 1>;
    using Matrix3s = Eigen::Matrix<typename state_policy::scalar_type, 3, 3>;
    Matrix3s _Rib;
    bool _initialised;

public:
    enum { input_length = state_policy::state_vector_length, output_length = 3 };
    using output_scalar_type = float;
    using input_scalar_type = typename state_policy::scalar_type;
    using YVector = Vector3s;
    using XVector = Eigen::Matrix<input_scalar_type, input_length, 1>;

	gyroscope()
		: _Rib(Matrix3s::Identity())
		, _initialised(false)
    {}

    void initialise(const Matrix3s& Rib)
    {
        _Rib = Rib;
        _initialised = true;
    }

    bool initialised() const { return _initialised; }

    YVector operator() (const XVector& x) const
    {
        const auto omegaBNb = state_policy::omegaBNb(x);
        const auto gbBNi = state_policy::gbBNi(x);

        YVector y = _Rib*(omegaBNb) + gbBNi;

        return y;
    }

};

}
}
}
}

#endif // OPENGNC_ESTIMATION_MODELS_MEASUREMENT_GYROSCOPE_HPP
