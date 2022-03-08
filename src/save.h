#ifndef _SAVE_H

#include "DataFile.h"
#include "Function.h"
#include <string>

using namespace std;

class Save {
private :
    // Pointeur de la classe DataFile
    DataFile* _df;

    // Pointeur de la classe function
    Function* _fct;

    int _dim, _Nx, _Ny;
    double _Lx, _Ly, _dx, _dy, _x0;

public :
    // Constructeur de la classe
    Save(DataFile* data_file, Function* function);

    // Permet de sauvegarder la solution
    void Save_Data(const double* Vec, const std::string& file_name);

    // Sauvegarde la convergence
    void Save_convergence(const double* Vec, const int& nb_cv);
};

#define _SAVE_H
#endif
