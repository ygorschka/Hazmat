#ifndef _SAVE_CPP

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "save.h"

using namespace std;

// Constructeur de la classe
Save::Save(DataFile* data_file, Function* function) :
    _df(data_file), _fct(function)
{
    _dim = _df->Get_dim();
    _Lx = _df->Get_Lx();
    _Ly = _df->Get_Ly();
    _x0 = _df->Get_x0();
}

// Permet de sauvegarder la solution
void Save::Save_Data(const double* Vec, const std::string& file_name)
{
    _Nx = _df->Get_Nx();
    _Ny = _df->Get_Ny();
    _dx = _Lx/(_Nx+1);
    _dy = _Ly/(_Ny+1);
    double energy;
    double f;
    std::vector<double> flux(_dim);
    ofstream mon_flux; // Contruit un objet "ofstream"
    string name_file(file_name); // Le nom de mon fichier
    mon_flux.open(name_file, ios::out); // Ouvre un fichier appelé name_file
    if(_dim == 2) {
        if(mon_flux) {// Vérifie que le fichier est bien ouvert
            for (int i = 0; i < _Nx; ++i) {
                for (int j = 0; j < _Ny; ++j) {
                    energy = Vec[3*(j*_Nx+i)];
                    flux[0] = Vec[3*(j*_Nx+i)+1];
                    flux[1] = Vec[3*(j*_Nx+i)+2];
                    f = _fct->anisotropy_factor(energy, &flux[0]);
                    mon_flux << std::fixed << std::setprecision(12) << (i+1)*_dx << " " <<
                    (j+1)*_dy << " " << energy << " " << flux[0] << " " << flux[1] << " "
                    << f << " " << endl;
                }
            }
        }
        else { // Renvoie un message d’erreur si ce n’est pas le cas
            printf("ERREUR: Impossible d’ouvrir le fichier.\n");
        }
        mon_flux.close(); // Ferme le fichier
    }
    else {
        if(mon_flux) {// Vérifie que le fichier est bien ouvert
            for (int i = 0; i < _Nx; ++i) {
                energy = Vec[2*i];
                flux[0] = Vec[2*i+1];
                mon_flux << std::fixed << std::setprecision(12) <<_x0+(i+1)*_dx << " "
                << energy << " "<<flux[0] << " " << endl;
            }
        }
        else { // Renvoie un message d’erreur si ce n’est pas le cas
            printf("ERREUR: Impossible d’ouvrir le fichier.\n");
        }
        mon_flux.close(); // Ferme le fichier
    }
}

// Permet de sauvegarder la convergence
void Save::Save_convergence(const double* Vec, const int& nb_cv)
{
    ofstream mon_flux; // Contruit un objet "ofstream"
    string name_file("Convergence.txt"); // Le nom de mon fichier
    mon_flux.open(name_file, ios::out); // Ouvre un fichier appelé name_file
    if(mon_flux) {// Vérifie que le fichier est bien ouvert
        for (int i = 0; i < nb_cv; i++) {
            mon_flux << std::fixed << std::setprecision(12) << Vec[3*i] << " " << Vec[3*i+1]
            << " " << Vec[3*i+2] << " " << " " << endl;
        }
    }
    else { // Renvoie un message d’erreur si ce n’est pas le cas
        printf("ERREUR: Impossible d’ouvrir le fichier.\n");
    }
    mon_flux.close(); // Ferme le fichier
}


#define _SAVE_CPP
#endif
