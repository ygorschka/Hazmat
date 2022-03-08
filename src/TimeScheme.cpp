#ifndef _TIME_SCHEME_CPP

#include "TimeScheme.h"
#include <math.h>

using namespace std;

// Constructeur de la classe
TimeScheme::TimeScheme(DataFile* data_file, Function* function, FiniteVolume* finite_volume) :
    _fin_vol(finite_volume), _df(data_file), _fct(function)
{
    _cfl = _df->Get_cfl();
    _Lx = _df->Get_Lx();
    _Ly = _df->Get_Ly();
    _dim = _df->Get_dim();
}

// Schema d'Euler explicite
void TimeScheme::Euler_Explicite(double* ErFr, const double& deltat)
{
    _Nx = _df->Get_Nx();
    _Ny = _df->Get_Ny();
    std::vector<double> vec, ErFr2;

    // Cas 1D
    if(_dim == 1) {
        int n = 2*_Nx;
        vec.resize(n);
        ErFr2.resize(n);
        vec = _fin_vol->Get_vec_flux();
        _fct->produit_cst(deltat, &vec[0], n);
        // Copie du vecteur ErFr
        for(int i=0; i<_Nx; i++) {
            ErFr2[2*i] = ErFr[2*i];
            ErFr2[2*i+1] = ErFr[2*i+1];
        }
        // Mise à jour du vecteur
        for(int i=0; i<_Nx; i++) {
            ErFr[2*i] = ErFr2[2*i] - vec[2*i];
            ErFr[2*i+1] = ErFr2[2*i+1] - vec[2*i+1];
        }
    }
    // Cas 2D
    else {
        int n = 3*_Nx*_Ny;
        vec.resize(n);
        //ErFr2.resize(n);
        vec = _fin_vol->Get_vec_flux();
        _fct->produit_cst(deltat, &vec[0], n);
        // Copie du vecteur ErFr
        /*for(int i=0; i<_Nx*_Ny; i++) {
            ErFr2[2*i] = ErFr[2*i];
            ErFr2[2*i+1] = ErFr[2*i+1];
            ErFr2[2*i+2] = ErFr[2*i+2];
        }*/
        // Mise à jour du vecteur
        for(int i=0; i<_Nx*_Ny; i++) {
            ErFr[3*i] = ErFr[3*i] - vec[3*i];
            ErFr[3*i+1] = ErFr[3*i+1] - vec[3*i+1];
            ErFr[3*i+2] = ErFr[3*i+2] - vec[3*i+2];
        }
    }
}

// Schema SSP RK2
void TimeScheme::RK2(double* ErFr, const double& deltat)
{
    _Nx = _df->Get_Nx();
    _Ny = _df->Get_Ny();
    std::vector<double> vec, ErFr2, k1, k2;

    // Cas 1D
    if(_dim == 1) {
        int n = 2*_Nx;
        vec.resize(n);
        k1.resize(n);
        k2.resize(n);
        vec = _fin_vol->Get_vec_flux();
        _fct->produit_cst(deltat, &vec[0], n);
        //Remplissage de k1
        for(int i=0; i<_Nx; i++) {
            k1[2*i] = ErFr[2*i] - vec[2*i];
            k1[2*i+1] = ErFr[2*i+1] - vec[2*i+1];
        }
        // Calcul de k2
        _fin_vol->Vector_Flux(&k1[0]);
        vec = _fin_vol->Get_vec_flux();
        _fct->produit_cst(deltat, &vec[0], n);
        //Remplissage de k2
        for(int i=0; i<_Nx; i++) {
            k2[2*i] = 0.75*ErFr[2*i] +0.25*k1[2*i] -0.25*vec[2*i];
            k2[2*i+1] = 0.75*ErFr[2*i+1] +0.25*k1[2*i+1] -0.25*vec[2*i+1];
        }
        _fin_vol->Vector_Flux(&k2[0]);
        vec = _fin_vol->Get_vec_flux();
        _fct->produit_cst(deltat, &vec[0], n);
        // Mise à jour du vecteur
        double tier = 1.0/3.0;
        for(int i=0; i<_Nx; i++) {
            ErFr[2*i] = tier*ErFr[2*i] + tier*k2[2*i] - deltat*tier*vec[2*i];
            ErFr[2*i+1] = tier*ErFr[2*i+1] + tier*k2[2*i+1] - deltat*tier*vec[2*i+1];
        }
    }

}

// Calcul du pas de temps dt
double TimeScheme::Get_dt()
{
    _Nx = _df->Get_Nx();
    _Ny = _df->Get_Ny();
    _dx = _Lx/(_Nx+1);
    _dy = _Ly/(_Ny+1);
    int n;
    std::vector<double> vec_coeff;
    if(_dim == 1){
        n = _Nx;
        vec_coeff.resize(n);
    }
    else {
        n = _Nx*_Ny;
        vec_coeff.resize(n);
    }
    // Vecteur mis à jour dans la boucle en temps dans le main
    vec_coeff = _fct->Get_vec_coeff();
    double coeff;
    coeff = _fct->max_vec(&vec_coeff[0], n);
    if (_dim == 1) {
        _dt = _cfl*_dx/(2.0*coeff);
    }
    else {
        _dt = _cfl*_dy*_dx/(2.0*coeff);
    }


    return _dt;
}

#define _TIME_SCHEME_CPP
#endif
