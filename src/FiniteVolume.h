#ifndef _FINITE_VOLUME

#include "function.h"
#include "DataFile.h"
#include <vector>

using namespace std;

class FiniteVolume {
private:
    // Pointeur vers la classe fonction
    Function* _fct;
    // Pointeur vers la classe DataFile
    DataFile* _df;

    // Vecteurs
    std::vector<double> _Pr, _vec_flux;

    int _Nx, _Ny, _dim, _space_ordre;
    double _c, _dx, _dy, _Lx, _Ly;

public:
    // Constructeur de la classe
    FiniteVolume(DataFile* data_file, Function* function);

    // Construit du vecteur flux
    void Vector_Flux(const double* ErFr);
    void Stencil1(const int& ind, const double* ErFr, const int& indu, const int& indd,
                  const int& indl, const int& indr);
    // Calcul le flux
    std::vector<double> Flux(const double* Pr, const double* ErFr, const double* n, const int& ind1,
                             const int& ind2, const double& coeff, const int& lr);

    // Renvoie le vecteur flux
    const std::vector<double>& Get_vec_flux() const {return _vec_flux;};

    // Construit le vecteur Pr
    void build_Pr(const double* ErFr);
    // Retourne le vecteur Pr
    const std::vector<double>& Get_Pr() const {return _Pr;};

};

#define _FINITE_VOLUME
#endif
