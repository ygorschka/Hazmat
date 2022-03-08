#ifndef _TIME_SCHEME_H

#include "DataFile.h"
#include "Function.h"
#include "FiniteVolume.h"
#include <vector>

class TimeScheme {
private:
    // Pointeur de la classe FiniteVolume
    FiniteVolume* _fin_vol;
    // Pointeur de la classe DataFile
    DataFile* _df;
    // Pointeur de la classe function
    Function* _fct;

    int _Nx, _dim, _Ny;
    double _dt, _dx, _coeff_rusanov, _Lx, _cfl, _Ly, _dy;

public:
    // Constructeur de la classe
    TimeScheme(DataFile* data_file, Function* function, FiniteVolume* finite_volume);

    // Schema d'Euler explicite
    void Euler_Explicite(double* ErFr, const double& deltat);

    // Schema SSP RK2
    void RK2(double* ErFr, const double& deltat);

    // Calcul le pas de temps dt
    double Get_dt();
};

#define _TIME_SCHEME_H
#endif
