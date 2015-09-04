#ifndef NUGNC_PIDCONTROLLER_H
#define NUGNC_PIDCONTROLLER_H

#include "GNC/util/types.h"

namespace GNC
{
class PIDController
{

protected:
    Scalar P;
    Scalar I;
    Scalar D;
    Scalar N;

    Matrix2s Ad;
    Vector2s Bd;
    Vector2s Cd;
    Scalar Cinf;

    Scalar slewMax;
    Scalar slewMin;
    Scalar satMax;
    Scalar satMin;

    Scalar uPrev;
    Vector2s z;
    Scalar T;

public:
    PIDController();            // default P = 1, dT = 0.01
    PIDController(Scalar _P,
                    Scalar _I,
                    Scalar _D,
                    Scalar _N,
                    Scalar _timestep,
                    Scalar _slewMax,
                    Scalar _slewMin,
                    Scalar _satMax,
                    Scalar _satMin);
    PIDController(Scalar _P,
                    Scalar _I,
                    Scalar _D,
                    Scalar _N,
                    Scalar _timestep);      // PID controller
    PIDController(Scalar _P,
                    Scalar _D,
                    Scalar _N,
                    Scalar _timestep);      // PD controller
    PIDController(Scalar _P,
                    Scalar _I,
                    Scalar _timestep);      // PI controller
    PIDController(Scalar _P,
                    Scalar _timestep);      // P controller

    void initialise(Scalar _P,
                    Scalar _I,
                    Scalar _D,
                    Scalar _N,
                    Scalar _timestep,
                    Scalar _slewMax,
                    Scalar _slewMin,
                    Scalar _satMax,
                    Scalar _satMin);
    void initialise(Scalar _P,
                    Scalar _I,
                    Scalar _D,
                    Scalar _N,
                    Scalar _timestep);

    void setSlewRates(Scalar _slewMax, Scalar _slewMin);
    void setSaturation(Scalar _satMax, Scalar _satMin);

    Scalar operator ()(Scalar _x, Scalar _ref);
    Scalar operator ()(Scalar _x, Scalar _ref, Scalar _dT);
    Vector2s& getState()  {return z; }
    Scalar& getUPrev()  {return uPrev; }

protected:
    // (1/C - 1/Cinf) * u_prev
    Scalar control(const Scalar &_u);

    void slew(Scalar& _u);
    void staturate(Scalar &_u);
    void calculateParameters();

};




}

#endif // PIDCONTROLLER_H

