#ifndef _DATA_FILE_CPP

#include "DataFile.h"
#include <iostream>
#include <fstream>

using namespace std;

// Constructeur de la classe DataFile
DataFile::DataFile(std::string file_name):
    _file_name(file_name)
{
    // Constantes physiques
    //_c = 299792458; // Vitesse de la lumière dans le vide
    //_c = 1.0;
    _a = 7.56573085e-16; // Constante de radiation
}

void DataFile::Read_data_file()
{
    ifstream data_file(_file_name.data());
    if (!data_file.is_open())
    {
        printf("Unable to open file\n");
        exit(0);
    }
    else
    {
        printf("-------------------------------------------------\n");
        printf("Reading data file\n");
    }

    data_file >> _Nx; // 1 Nombre de points selon x
    data_file >> _Ny; // 2 Nombre de points selon y
    data_file >> _Lx; // 3 Longueur selon x
    data_file >> _Ly; // 4 Longueur selon y
    data_file >> _x0; // 5 x0
    data_file >> _y0; // 6 y0
    data_file >> _order_space; // 7 Ordre du schéma numérique en espace
    data_file >> _order_time; // 8 Ordre du schéma numérique en temps
    data_file >> _dim; // 9 Dimension de l'espace
    data_file >> _Tf; // 10 Temps final de simulation
    data_file >> _cfl; // 11 Condition de CFL
    data_file >> _case; // 12 Cas de la simulation
    data_file >> _nb_simulation; // 13 Nombre de simulations pouir le calcul de l'erreur
    data_file >> _c; // 14 Vitesse de la lumière
}

// Change la valeur de Nx
void DataFile::change_Nx(const int& new_nx)
{
    _Nx = new_nx;
}

// Change la valeur de Ny
void DataFile::change_Ny(const int& new_ny)
{
    _Ny = new_ny;
}

#define _DATA_FILE_CPP
#endif
