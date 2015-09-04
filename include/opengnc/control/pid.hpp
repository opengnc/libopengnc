#ifndef PID_H
#define PID_H

#include <Eigen/Core>

namespace opengnc {
namespace control {

template <typename scalar_type>
class pid
{
protected:
	using Matrix2s = Eigen::Matrix<scalar_type, 2, 2>;
	using Vector2s = Eigen::Matrix<scalar_type, 2, 1>;

	scalar_type P;
	scalar_type I;
	scalar_type D;
	scalar_type N;

	Matrix2s Ad;
	Vector2s Bd;
	Vector2s Cd;
	scalar_type Cinf;

	scalar_type slewMax;
	scalar_type slewMin;
	scalar_type satMax;
	scalar_type satMin;

	scalar_type uPrev;
	Vector2s z;
	scalar_type T;

public:
	pid()            // default P = 1, dT = 0.01
		: P(1)
		, I(0)
		, D(0)
		, N(100)
		, Ad(Matrix2s::Zero())
		, Bd(Vector2s::Zero())
		, Cd(Vector2s::Zero())
		, Cinf(0)
		, slewMax(std::numeric_limits<scalar_type>::infinity())
		, slewMin(-std::numeric_limits<scalar_type>::infinity())
		, satMax(std::numeric_limits<scalar_type>::infinity())
		, satMin(-std::numeric_limits<scalar_type>::infinity())
		, uPrev(0)
		, z(Vector2s::Zero())
		, T(0.02)
	{
		calculateParameters();
	}

	pid(scalar_type _P,
		scalar_type _I,
		scalar_type _D,
		scalar_type _N,
		scalar_type _timestep,
		scalar_type _slewMax,
		scalar_type _slewMin,
		scalar_type _satMax,
		scalar_type _satMin)
	{
		initialise(_P, _I, _D, _N, _timestep, _slewMax, _slewMin, _satMax, _satMin);
	}

	pid(scalar_type _P,
		scalar_type _I,
		scalar_type _D,
		scalar_type _N,
		scalar_type _timestep)
	{
		initialise(_P, _I, _D, _N, _timestep);
	}

	pid(scalar_type _P,
		scalar_type _D,
		scalar_type _N,
		scalar_type _timestep)
	{
		initialise(_P, 0.0, _D, _N, _timestep);
	}

	pid(scalar_type _P,
		scalar_type _I,
		scalar_type _timestep)
	{
		initialise(_P, _I, 0.0, 0.0, _timestep);
	}

	pid(scalar_type _P,
		scalar_type _timestep)
	{
		initialise(_P, 0.0, 0.0, 0.0, _timestep);
	}

	void initialise(scalar_type _P,
					scalar_type _I,
					scalar_type _D,
					scalar_type _N,
					scalar_type _timestep,
					scalar_type _slewMax,
					scalar_type _slewMin,
					scalar_type _satMax,
					scalar_type _satMin)
	{
		P = _P;
		I = _I;
		D = _D;
		N = _N;
		slewMax = _slewMax;
		slewMin = _slewMin;
		satMax = _satMax;
		satMin = _satMin;
		T = _timestep;

		calculateParameters();
	}

	void initialise(scalar_type _P,
					scalar_type _I,
					scalar_type _D,
					scalar_type _N,
					scalar_type _timestep)
	{
		P = _P;
		I = _I;
		D = _D;
		N = _N;
		T = _timestep;

		calculateParameters();
	}

	void setSlewRates(scalar_type _slewMax, scalar_type _slewMin)
	{
		slewMax = _slewMax;
		slewMin = _slewMin;
	}

	void setSaturation(scalar_type _satMax, scalar_type _satMin)
	{
		satMax = _satMax;
		satMin = _satMin;
	}

	scalar_type operator ()(scalar_type _x, scalar_type _ref, scalar_type _dT)
	{
		if (_dT != T)
		{
			T = _dT;
			calculateParameters();
		}

		return operator ()(_x, _ref);
	}

	scalar_type operator ()(scalar_type _x, scalar_type _ref)
	{
		scalar_type error = _ref - _x;
		scalar_type y = control(uPrev);

		scalar_type Ct = (error - y);
		scalar_type uHat = Cinf*Ct;

		slew(uHat);
		saturate(uHat);
		uPrev = uHat;

		return uHat;
	}

	Vector2s& getState()  { return z; }
	scalar_type& getUPrev()  { return uPrev; }

protected:
    // (1/C - 1/Cinf) * u_prev
	scalar_type control(const scalar_type &_u)
	{
		float y = Cd.transpose()*z;
		z = Ad*z + Bd*_u;
		return y;
	}

	void slew(scalar_type& _u)
	{
		scalar_type rate = (_u-uPrev)/T;

		if (rate > slewMax) {
			rate = T*slewMax + uPrev;
		}
		else if (rate < slewMin) {
			_u = T*slewMin + uPrev;
		}
	}

	void saturate(scalar_type &_u)
	{
		if (_u > satMax) {
			_u = satMax;
		}
		else if (_u < satMin) {
			_u = satMin;
		}
	}

	void calculateParameters()
	{
		Cinf  = P + D*N;

		Ad << 1-(T*(P*N + I))/Cinf, -T*I*N/Cinf,
			  T, 1;
		Bd << T, 0;
		Cd << (D*std::pow(N,2) - I)/std::pow(Cinf, 2),
			  -I*N/std::pow(Cinf, 2);
	}

};

}
}

#endif // PID_H
