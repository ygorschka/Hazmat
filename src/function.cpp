#ifndef _FUNCTION_CPP

#include "function.h"
#include <math.h>
#include <cmath>
#include <vector>
#include <string>

// Constructeur de la classe function
Function::Function(DataFile* data_file) :
    _df(data_file)
{
    _dim = _df->Get_dim();
    _c = _df->Get_c();
    _Lx = _df->Get_Lx();
    _a = _df->Get_a();
    _case = _df->Get_case();
    _x0 = _df->Get_x0();
    _nb_simulation = _df->Get_nb_simulation();
    _vec_error.resize(3*_nb_simulation);
}

// Calcul le facteur d'eddington
double Function::eddington_factor(const double& x)
{
    double x2 = pow(x,2);
    _khi = (3.0 + 4.0*x2)/(5.0+2.0*sqrt(4.0-3.0*x2));

    return _khi;
}

// Calcul la norme 2 d'un vecteur
double Function::norme_2(const double* vec)
{
    double sum = 0.0;
    for(int i=0; i<_dim; i++) {
        sum  += pow(vec[i],2);
    }
    _norme = sqrt(sum);

    return _norme;
}

// Calcul le facteur d'anisotropie
double Function::anisotropy_factor(const double& Er, const double* Fr)
{
    double nor;
    nor = norme_2(&Fr[0]);
    _f = nor/(_c*Er);

    return _f;
}

// Calcul le coefficient de Rusanov
double Function::coeff_rusanov(const int& i, const double* ErFr)
{
    _Nx = _df->Get_Nx();
    _Ny = _df->Get_Ny();
    if(_dim == 1) {
        _vec_coeff.resize(_Nx);
        _vec_coeff[i] = _c;
        return _c;
    }
    else {
        _vec_coeff.resize(_Nx*_Ny);
        /*double coeffx, coeffy, coeff, energy, f, khi;
        std::vector<double> flux(2);
        flux[0] = ErFr[3*i+1];
        flux[1] = ErFr[3*i+2];
        energy = ErFr[3*i];
        f = anisotropy_factor(energy, &flux[0]);
        khi = eddington_factor(f);
        coeffx = _c*((3.0*khi-1.0)/2.0)*(flux[0]/norme_2(&flux[0]))*_c*energy;
        coeffy = _c*((3.0*khi-1.0)/2.0)*(flux[1]/norme_2(&flux[0]))*_c*energy;
        if(abs(coeffx) > abs(coeffy)) {
            _vec_coeff[i] = abs(coeffx);
            return abs(coeffx);
        }
        else {
            _vec_coeff[i] = abs(coeffy);
            return abs(coeffy);
        }*/
        /*double coeffx, coeffy;
        std::vector<double> flux(2);
        flux[0] = ErFr[3*i+1];
        flux[1] = ErFr[3*i+2];
        coeffx = _c*abs(flux[0]);
        coeffy = _c*abs(flux[1]);
        if(coeffx > coeffy) {
            _vec_coeff[i] = coeffx;
            return coeffx;
        }
        else {
            _vec_coeff[i] = coeffy;
            return coeffy;
        }*/
        _vec_coeff[i] = _c;
        return _c;
    }

}

// Calcul le produit d'un vecteur par une constante
void Function::produit_cst(const double a, double* b, const int& size)
{
    for (int i = 0; i < size; i++) {
        b[i]=a*b[i];
    }
};

// Renvoie la somme de 2 vecteurs
std::vector<double> Function::vector_sum(const double* vec1, const double* vec2, const int& size)
{
    std::vector<double> res(size, 0.0);

    for(int i=0; i<size; i++) {
        res[i] = vec1[i]+vec2[i];
    }

    return res;
}

// Renvoie le max d'un vecteur donné
double Function::max_vec(const double* vec, const int& size)
{
    double max= vec[0];

    for(int i=0; i<size; i++) {
        if (vec[i] > max) {
            max = vec[i];
        }
    }
    return max;
}

// Fonction minmod
double Function::minmod(const double& x, const double& y)
{
    double res = 0.0;
    if(x*y >= 0) {
        if(abs(x) < abs(y)) {
            res = x;
        }
        else {
            res = y;
        }
    }

    return res;
}

// Fonction Van Leer
double Function::van_leer(const double& x)
{
    double res;

    res = (x + abs(x))/(1.0 + abs(x));

    return res;
}

// Calcule l'erreur pour deux vecteurs donnés
void Function::Caculate_error(const double* vec_1, const double* vec_2, const int& nb_points_x,
                              const int& nb_points_y, const int& ind, const double& deltax,
                              const double& deltay)
{
    double err1 = 0.0;
    double err2 = 0.0;

    if(_dim == 1) {
        int ab = pow(2,ind+1);

        for(int i=0; i<nb_points_x; i++) {
            err1 = err1 + pow(vec_2[2*ab*i+2*ab-2]-vec_1[ab*i+ab-2],2);
            err2 = err2 + pow(vec_2[2*ab*i+2*ab-1]-vec_1[ab*i+ab-1],2);
        }
        err1 = deltax*err1;
        err1 = sqrt(err1);
        err2 = deltax*err2;
        err2 = sqrt(err2);

        _vec_error[3*ind] = deltax;
        _vec_error[3*ind+1] = err1;
        _vec_error[3*ind+2] = err2;
    }
    else {
        int ab = pow(2,ind);
        int ab2 = pow(2,ind+1);
        int ind_i_1, ind_j_1, ind_Nx_1, ind_g_1;
        int ind_i_2, ind_j_2, ind_Nx_2, ind_g_2;
        ind_Nx_1 = ab*nb_points_x+ab-1;
        ind_Nx_2 = ab2*nb_points_x+ab2-1;

        for(int i=0; i<nb_points_x; i++) {
            ind_i_1 = ab*i+ab-1;
            ind_i_2 = ab2*i+ab2-1;
            for(int j=0; j<nb_points_y; j++) {
                ind_j_1 = ab*j+ab-1;
                ind_j_2 = ab2*j+ab2-1;
                ind_g_1 = ind_j_1*ind_Nx_1+ind_i_1;
                ind_g_2 = ind_j_2*ind_Nx_2+ind_i_2;
                err1 = err1 + pow(vec_2[3*ind_g_2]-vec_1[3*ind_g_1],2);
                err2 = err2 + pow(vec_2[3*ind_g_2+1]-vec_1[3*ind_g_1+1],2)*deltax +
                       pow(vec_2[3*ind_g_2+2]-vec_1[3*ind_g_1+2],2)*deltay;
            }
        }

        err1 = deltax*err1;
        err1 = sqrt(err1);
        err2 = sqrt(err2);

        _vec_error[3*ind] = deltax;
        _vec_error[3*ind+1] = err1;
        _vec_error[3*ind+2] = err2;

    }

}

// Renvoie la donnée initiale
void Function::init(double* ErFr)
{
    _Nx = _df->Get_Nx();
    _Ny = _df->Get_Ny();
    _dx = _Lx/(_Nx+1);
    int ind;

    if (_dim == 1) {
        if(_case == "cas1") {
            for (int i = 0; i<_Nx; i++) {
                if(_x0 + (i+1)*_dx<0.0) {
                    ErFr[2*i] = 1.0;
                    ErFr[2*i+1] = 0.1;
                }
                else {
                    ErFr[2*i] = 10.0;
                    ErFr[2*i+1] = -1.0;
                }
            }
        }
        else if(_case == "cas2") {
            for (int i = 0; i<_Nx; i++) {
                if(_x0 + (i+1)*_dx<0.0) {
                    ErFr[2*i] = 2.0;
                    ErFr[2*i+1] = -0.5;
                }
                else {
                    ErFr[2*i] = 2.0;
                    ErFr[2*i+1] = 1.0;
                }
            }
        }
    }
    else if(_dim == 2) {
        if(_case == "cas1") {
            double T0 = 1000.0;
            double e0 = _a*pow(T0,4);
            double f0 = (1.0-10E-8)*_c*e0;
            for(int i=0; i<_Nx; i++) {
                for(int j=0; j<_Ny; j++) {
                    ind = j*_Nx+i;
                    ErFr[3*ind] = e0;
                    if(i<_Nx/2) {
                        if(j<_Ny/2) {
                            ErFr[3*ind+1] = f0;
                            ErFr[3*ind+2] = 0;
                        }
                        else {
                            ErFr[3*ind+1] = 0;
                            ErFr[3*ind+2] = f0;
                        }
                    }
                    else{
                        if(j<_Ny/2) {
                            ErFr[3*ind+1] = 0;
                            ErFr[3*ind+2] = -f0;
                        }
                        else {
                            ErFr[3*ind+1] = -f0;
                            ErFr[3*ind+2] = 0;
                        }
                    }
                }
            }
        }
        else if(_case == "cas2") {
            double e0 = 4.0;
            double f0 = 1.0;
            for(int i=0; i<_Nx; i++) {
                for(int j=0; j<_Ny; j++) {
                    ind = j*_Nx+i;
                    ErFr[3*ind] = e0;
                    if(i<_Nx/2) {
                        if(j<_Ny/2) {
                            ErFr[3*ind+1] = f0;
                            ErFr[3*ind+2] = 0;
                        }
                        else {
                            ErFr[3*ind+1] = 0;
                            ErFr[3*ind+2] = f0;
                        }
                    }
                    else{
                        if(j<_Ny/2) {
                            ErFr[3*ind+1] = 0;
                            ErFr[3*ind+2] = -f0;
                        }
                        else {
                            ErFr[3*ind+1] = -f0;
                            ErFr[3*ind+2] = 0;
                        }
                    }
                }
            }
        }
        else if (_case == "cas3") {
            double e0 = 4.0;
            double f0 = 0.0;
            for(int i=0; i<_Nx; i++) {
                for(int j=0; j<_Ny; j++) {
                    ind = j*_Nx+i;
                    ErFr[3*ind] = e0;
                    ErFr[3*ind+1] = f0;
                    ErFr[3*ind+2] = f0;
                }
            }
        }
    }
}

#define _FUNCTION_CPP
#endif
