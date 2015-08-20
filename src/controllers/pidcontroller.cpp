#include "GNC/controller/pidcontroller.h"
#include <limits.h>

using namespace GNC;

PIDController::PIDController()
    : P(1)
    , I(0)
    , D(0)
    , N(100)
    , Ad(Matrix2s::Zero())
    , Bd(Vector2s::Zero())
    , Cd(Vector2s::Zero())
    , Cinf(0)
    , slewMax(std::numeric_limits<Scalar>::infinity())
    , slewMin(-std::numeric_limits<Scalar>::infinity())
    , satMax(std::numeric_limits<Scalar>::infinity())
    , satMin(-std::numeric_limits<Scalar>::infinity())
    , uPrev(0)
    , z(Vector2s::Zero())
    , T(0.02)
{
    calculateParameters();
}

PIDController::PIDController(Scalar _P,
                             Scalar _I,
                             Scalar _D,
                             Scalar _N,
                             Scalar _timestep,
                             Scalar _slewMax,
                             Scalar _slewMin,
                             Scalar _satMax,
                             Scalar _satMin)
{
    initialise(_P, _I, _D, _N, _timestep, _slewMax, _slewMin, _satMax, _satMin);
}
// PID controller
PIDController::PIDController(Scalar _P,
                Scalar _I,
                Scalar _D,
                Scalar _N,
                Scalar _timestep)
{
    initialise(_P, _I, _D, _N, _timestep);
}
// PD controller
PIDController::PIDController(Scalar _P,
                Scalar _D,
                Scalar _N,
                Scalar _timestep)
{
    initialise(_P, 0.0, _D, _N, _timestep);
}
// PI controller
PIDController::PIDController(Scalar _P,
                Scalar _I,
                Scalar _timestep)
{
    initialise(_P, _I, 0.0, 0.0, _timestep);
}
// P controller
PIDController::PIDController(Scalar _P,
                Scalar _timestep)
{
    initialise(_P, 0.0, 0.0, 0.0, _timestep);
}

void PIDController::initialise(Scalar _P,
                               Scalar _I,
                               Scalar _D,
                               Scalar _N,
                               Scalar _timestep,
                               Scalar _slewMax,
                               Scalar _slewMin,
                               Scalar _satMax,
                               Scalar _satMin)
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

void PIDController::initialise(Scalar _P,
                               Scalar _I,
                               Scalar _D,
                               Scalar _N,
                               Scalar _timestep)
{
    P = _P;
    I = _I;
    D = _D;
    N = _N;
    T = _timestep;

    calculateParameters();
}

void PIDController::setSlewRates(Scalar _slewMax, Scalar _slewMin)
{
    slewMax = _slewMax;
    slewMin = _slewMin;
}
void PIDController::setSaturation(Scalar _satMax, Scalar _satMin)
{
    satMax = _satMax;
    satMin = _satMin;
}

Scalar PIDController::operator ()(Scalar _x, Scalar _ref, Scalar _dT)
{
    if (_dT != T)
    {
        T = _dT;
        calculateParameters();
    }

    return operator ()(_x, _ref);
}

Scalar PIDController::operator ()(Scalar _x, Scalar _ref)
{
    Scalar error = _ref - _x;
    Scalar y = control(uPrev);

    Scalar Ct = (error - y);
    Scalar uHat = Cinf*Ct;

    slew(uHat);
    staturate(uHat);
    uPrev = uHat;

    return uHat;
}

Scalar PIDController::control(const Scalar &_u)
{
    float y = Cd.transpose()*z;
    z = Ad*z + Bd*_u;
    return y;
}

void PIDController::slew(Scalar &_u)
{
    Scalar rate = (_u-uPrev)/T;

    if (rate > slewMax) {
        rate = T*slewMax + uPrev;
    }
    else if (rate < slewMin) {
        _u = T*slewMin + uPrev;
    }
}

void PIDController::staturate(Scalar &_u)
{
    if (_u > satMax) {
        _u = satMax;
    }
    else if (_u < satMin) {
        _u = satMin;
    }
}

void PIDController::calculateParameters()
{
    Cinf  = P + D*N;

    Ad << 1-(T*(P*N + I))/Cinf, -T*I*N/Cinf,
          T, 1;
    Bd << T, 0;
    Cd << (D*std::pow(N,2) - I)/std::pow(Cinf, 2),
          -I*N/std::pow(Cinf, 2);
}
