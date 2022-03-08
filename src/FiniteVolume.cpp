#ifndef _FINITE_VOLUME_CPP

#include "FiniteVolume.h"
#include <vector>
#include <math.h>

using namespace std;

// Constructeur de la classe FiniteVolume
FiniteVolume::FiniteVolume(DataFile* data_file, Function* function) :
    _df(data_file), _fct(function)
{
    _dim = _df->Get_dim();
    _Lx = _df->Get_Lx();
    _Ly = _df->Get_Ly();
    _c = _df->Get_c();
    _space_ordre = _df->Get_order_space();
}

// Construit le vecteur Pr
void FiniteVolume::build_Pr(const double* ErFr)
{
    _Nx = _df->Get_Nx();
    _Ny = _df->Get_Ny();
    double khi, f, energy, norme, coeff1, coeff2;
    std::vector<double> flux(_dim);

    if(_dim == 1) {
        _Pr.resize(_Nx);
        for(int i=0; i<_Nx; i++) {
            flux[0] = ErFr[2*i+1];
            energy = ErFr[2*i];
            norme = _fct->norme_2(&flux[0]);
            // Comparaison avec l'erreur machine
            if(norme > 10E-16) {
                f = _fct->anisotropy_factor(energy, &flux[0]);
                khi = _fct->eddington_factor(f);
                // Test de différents cas de figures
                //_Pr[i] = energy*((1.0-khi)/2.0 + ((3.0*khi-1.0)/2.0)*(pow(flux[0],2)/(norme));
                _Pr[i] = energy*khi;

            }
            else {
                _Pr[i] = energy/3.0;
            }
        }
    }
    else {
        _Pr.resize(4*_Nx*_Ny);
        int ind;
        for(int i=0; i<_Nx; i++) {
            for(int j=0; j<_Ny; j++) {
                ind = j*_Nx+i;
                flux[0] = ErFr[3*ind+1];
                flux[1] = ErFr[3*ind+2];
                energy = ErFr[3*ind];
                norme = _fct->norme_2(&flux[0]);
                norme = pow(norme,2);
                // Comparaison avec l'erreur machine
                if(norme > 10E-16) {
                    f  = _fct->anisotropy_factor(energy, &flux[0]);
                    khi = _fct->eddington_factor(f);
                    coeff1 = (1.0-khi)/2.0;
                    coeff2 = (3.0*khi-1.0)/2.0;
                    _Pr[4*ind] = energy*(coeff1 + coeff2*(pow(flux[0],2)/norme));
                    _Pr[4*ind+1] = energy*(coeff2*flux[0]*flux[1]/norme);
                    _Pr[4*ind+2] = energy*(coeff2*flux[0]*flux[1]/norme);
                    _Pr[4*ind+3] = energy*(coeff1 + coeff2*(pow(flux[1],2)/norme));
                }
                else {
                    _Pr[4*ind] = energy/3.0;
                    _Pr[4*ind+1] = 0.0;
                    _Pr[4*ind+2] = 0.0;
                    _Pr[4*ind+3] = energy/3.0;
                }
            }
        }
    }
}

// Calcul du vecteur flux
void FiniteVolume::Vector_Flux(const double* ErFr)
{
    _Nx = _df->Get_Nx();
    _Ny = _df->Get_Ny();
    _dx = _Lx/(_Nx+1);
    _dy = _Ly/(_Ny+1);
    // Cas 1D
    if(_dim == 1) {

        _vec_flux.resize(2*_Nx);
        std::vector<double> Pr(_Nx, 0.0);
        std::vector<double> n(1);
        std::vector<double> flux1(2), flux2(2);
        double coeff, alpha;
        alpha = 1.0/_dx;

        // Recupère le vecteur Pr
        build_Pr(&ErFr[0]);
        Pr = Get_Pr();

        // Bord droit
        coeff = _fct->coeff_rusanov(0, &ErFr[0]);
        n[0] = -1.0;
        flux1 = Flux(&Pr[0], &ErFr[0], &n[0], 0, 0, coeff, 1);
        n[0] = 1.0;
        flux2 = Flux(&Pr[0], &ErFr[0], &n[0], 0, 1, coeff, 0);
        _vec_flux[0] = alpha*(flux1[0]+flux2[0]);
        _vec_flux[1] = alpha*(flux1[1]+flux2[1]);

        // Bord gauche
        coeff = _fct->coeff_rusanov(_Nx-1, &ErFr[0]);
        n[0] = -1.0;
        flux1 = Flux(&Pr[0], &ErFr[0], &n[0], _Nx-1, _Nx-2, coeff, 1);
        n[0] = 1.0;
        flux2 = Flux(&Pr[0], &ErFr[0], &n[0], _Nx-1, _Nx-1, coeff, 0);
        _vec_flux[2*_Nx-2] = alpha*(flux1[0]+flux2[0]);
        _vec_flux[2*_Nx-1] = alpha*(flux1[1]+flux2[1]);

        // Boucle pour faire le reste
        for(int i=1; i<_Nx-1; i++) {
            coeff = _fct->coeff_rusanov(i, &ErFr[0]);
            n[0] = -1.0;
            flux1 = Flux(&Pr[0], &ErFr[0], &n[0], i, i-1, coeff, 1);
            n[0] = 1.0;
            flux2 = Flux(&Pr[0], &ErFr[0], &n[0], i, i+1, coeff, 0);
            _vec_flux[2*i] = alpha*(flux1[0]+flux2[0]);
            _vec_flux[2*i+1] = alpha*(flux1[1]+flux2[1]);
        }

    }
    else {

        // Cas 2D
        _vec_flux.resize(3*_Nx*_Ny);
        int ind;

        // Remplie le vecteur flux
        if(_df->Get_order_space() == 1) {
            for(int i=0; i<_Nx; i++) {
                for(int j=0; j<_Ny; j++) {
                    ind = j*_Nx+i;
                    if((i == 0) && (j == 0)) {
                        Stencil1(ind, &ErFr[0], (j+1)*_Nx+i, j*_Nx+i, j*_Nx+i, j*_Nx+i+1);
                    }
                    else if((i == 0) && (j == _Ny-1)) {
                        Stencil1(ind, &ErFr[0], j*_Nx+i, (j-1)*_Nx+i, j*_Nx+i, j*_Nx+i+1);
                    }
                    else if((i == _Nx-1) && (j == _Ny-1)) {
                        Stencil1(ind, &ErFr[0], j*_Nx+i, (j-1)*_Nx+i, j*_Nx+i-1, j*_Nx+i);
                    }
                    else if((i == _Nx-1) && (j == 0)) {
                        Stencil1(ind, &ErFr[0], (j+1)*_Nx+i, j*_Nx+i, j*_Nx+i-1, j*_Nx+i);
                    }
                    else if(i == 0) {
                        Stencil1(ind, &ErFr[0], (j+1)*_Nx+i, (j-1)*_Nx+i, j*_Nx+i, j*_Nx+i+1);
                    }
                    else if(i == _Nx-1) {
                        Stencil1(ind, &ErFr[0], (j+1)*_Nx+i, (j-1)*_Nx+i, j*_Nx+i-1, j*_Nx+i);
                    }
                    else if(j == 0) {
                        Stencil1(ind, &ErFr[0], (j+1)*_Nx+i, j*_Nx+i, j*_Nx+i-1, j*_Nx+i+1);
                    }
                    else if(j == _Ny-1) {
                        Stencil1(ind, &ErFr[0], j*_Nx+i, (j-1)*_Nx+i, j*_Nx+i-1, j*_Nx+i+1);
                    }
                    else {
                        Stencil1(ind, &ErFr[0], (j+1)*_Nx+i, (j-1)*_Nx+i, j*_Nx+i-1, j*_Nx+i+1);
                    }
                }
            }
        }
    }
}

// Stencil d'ordre 1 dans le cas 2D
//Stencil1(ind, &ErFr[0], (j+1)*_Nx+i, (j-1)*_Nx+i, j*_Nx+i-1, j*_Nx+i+1);
void FiniteVolume::Stencil1(const int& ind, const double* ErFr, const int& indu, const int& indd,
                            const int& indl, const int& indr)
{
    _Nx = _df->Get_Nx();
    _Ny = _df->Get_Ny();
    _dx = _Lx/(_Nx+1);
    _dy = _Ly/(_Ny+1);
    std::vector<double> Pr(4*_Nx*_Ny, 0.0);
    std::vector<double> n(2);
    std::vector<double> flux1(3), flux2(3), flux3(3), flux4(3);
    double coeff, alpha;
    alpha = 1.0/(_dx*_dy);
    // Recupère le vecteur Pr
    build_Pr(&ErFr[0]);
    Pr = Get_Pr();
    coeff = _fct->coeff_rusanov(ind, &ErFr[0]);
    n[0] = 1.0, n[1] = 0.0;
    flux1 = Flux(&Pr[0], &ErFr[0], &n[0], ind, indr, coeff, 0);
    n[0] = -1.0, n[1] = 0.0;
    flux2 = Flux(&Pr[0], &ErFr[0], &n[0], ind, indl, coeff, 0);
    n[0] = 0.0, n[1] = 1.0;
    flux3 = Flux(&Pr[0], &ErFr[0], &n[0], ind, indu, coeff, 0);
    n[0] = 0.0, n[1] = -1.0;
    flux4 = Flux(&Pr[0], &ErFr[0], &n[0], ind, indd, coeff, 0);
    _vec_flux[3*ind] = alpha*(_dy*(flux1[0]+flux2[0])+_dx*(flux3[0]+flux4[0]));
    _vec_flux[3*ind+1] = alpha*(_dy*(flux1[1]+flux2[1])+_dx*(flux3[1]+flux4[1]));
    _vec_flux[3*ind+2] = alpha*(_dy*(flux1[2]+flux2[2])+_dx*(flux3[2]+flux4[2]));
}

// Calcul du flux-idée VF
// ind1 et ind2 : indices des mailles
std::vector<double> FiniteVolume::Flux(const double* Pr, const double* ErFr, const double* n,
                                       const int& ind1, const int& ind2, const double& coeff,
                                       const int& lr)
{
    std::vector<double> fl;

    if(_dim == 1) {
        if(_space_ordre == 1) {
            fl.resize(2);
            fl[0] = 0.5*(n[0]*(ErFr[2*ind1+1] + ErFr[2*ind2+1]) - coeff*(ErFr[2*ind2] - ErFr[2*ind1]));
            fl[1] = 0.5*(n[0]*pow(_c,2)*(Pr[ind1] + Pr[ind2]) - coeff*(ErFr[2*ind2+1] - ErFr[2*ind1+1]));
        }
    }
    else {
        fl.resize(3);
        fl[0] = 0.5*(n[0]*(ErFr[3*ind1+1]+ErFr[3*ind2+1]) + n[1]*(ErFr[3*ind1+2]+ErFr[3*ind2+2])) -
                0.5*coeff*(ErFr[3*ind2] - ErFr[3*ind1]);
        fl[1] = 0.5*pow(_c,2)*(n[0]*(Pr[4*ind2]+Pr[4*ind1]) + n[1]*(Pr[4*ind2+1]+Pr[4*ind1+1])) -
                0.5*coeff*(ErFr[3*ind2+1] - ErFr[3*ind1+1]);
        fl[2] = 0.5*pow(_c,2)*(n[0]*(Pr[4*ind2+2]+Pr[4*ind1+2]) + n[1]*(Pr[4*ind2+3]+Pr[4*ind1+3])) -
                0.5*coeff*(ErFr[3*ind2+2] - ErFr[3*ind1+2]);
    }
    return fl;

}

#define _FINITE_VOLUME_CPP
#endif
