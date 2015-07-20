#ifndef OPENGNC_COMMON_TRANSFORMS_WGS84_HPP
#define OPENGNC_COMMON_TRANSFORMS_WGS84_HPP

#include <Eigen/Core>

namespace opengnc {
namespace common {
namespace transforms {

class wgs84
{
private:
	static constexpr double a = 6378137.0;
	static constexpr double f = 0.00335281;
	static constexpr double b = 6356752.3142;
	static constexpr double e = 0.08181919;
	static constexpr double edash = std::sqrt((std::pow(wgs84::a, 2) - std::pow(wgs84::b, 2)) / (std::pow(wgs84::b, 2)));

public:
	static Eigen::Vector3d geodetic_to_rectangular(const Eigen::Vector3d& gps)
	{
		double phi = gps[0] * M_PI / 180;
		double lambda = gps[1] * M_PI / 180;
		double h = gps[2];

		double Rn = a / std::sqrt(1 - e*e * (std::sin(phi))*std::sin(phi));     // Normal Radius

		Eigen::Vector3d rBOe;
		rBOe << (Rn+h)*std::cos(phi)*std::cos(lambda),
				(Rn+h)*std::cos(phi)*std::sin(lambda),
				(Rn*(1-e*e) + h)*std::sin(phi);

		return rBOe;
	}

	static Eigen::Vector3d rectangular_to_geodetic(const Eigen::Vector3d rBOe)
	{
		// Position vector of body to origin in ECEF rectangular
		const double& xe = rBOe[0];
		const double& ye = rBOe[1];
		const double& ze = rBOe[2];

		Eigen::Vector3d gps;
		double& lattitude = gps[0];
		double& longitude = gps[1];
		double& altitude = gps[2];

		// Taken from Fundamentals of Inertial Navigation, Satellite based Positioning and their Integration.
		// p51
		const double p = sqrt(xe*xe + ye*ye);
		const double theta = std::atan2((ze * a), (p * b));

		lattitude = std::atan2((ze + edash*edash*b*pow(sin(theta), 3)), (p - e*e*a*pow(cos(theta), 3)));
		longitude = 2 * std::atan2(ye, (xe + p));

		double N = pow(a, 2) / sqrt(pow(a, 2) * pow(cos(lattitude), 2) + pow(b, 2) * pow(sin(lattitude), 2));
		altitude = p / cos(lattitude) - N;

		// Convert to degrees for GPS measurement
		lattitude *= 180 / M_PI;
		longitude *= 180 / M_PI;

		return gps;
	}

        template <typename scalar=double>
	static Eigen::Matrix<scalar,3,1> Ren_from_geodetic(const Eigen::Vector3d& gps)
	{
		const double& phi = gps[0] * M_PI / 180;
		const double& lambda = gps[1] * M_PI / 180;

		const double theta_z = lambda;
		const double theta_y = -M_PI/2 - phi;

		Eigen::Matrix<scalar,3,1> Rz, Ry, Ren;
		Rz <<   std::cos(theta_z),  std::sin(theta_z),  0,
				-std::sin(theta_z), std::cos(theta_z),  0,
				0, 0, 1;

		Ry <<   std::cos(theta_y),  0,  -std::sin(theta_y),
				0, 1, 0,
				std::sin(theta_y),  0,  std::cos(theta_y);

		Ren = Rz.transpose() * Ry.transpose();

		return Ren;
	}
};

}
}
}

#endif
