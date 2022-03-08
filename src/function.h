#ifndef _FUNCTION_H

#include "DataFile.h"
#include <vector>

class Function {
private:
    // Pointeur de la classe DataFile
    DataFile* _df;

    std::vector<double> _vec_coeff, _vec_error;
    std::string _case;

    double _khi, _f, _norme, _c, _a, _x0, _dx, _Lx;
    int _dim, _Nx, _Ny, _nb_simulation;

public:
    // Contructeur de la classe
    Function(DataFile* data_file);

    double eddington_factor(const double& x);
    double norme_2(const double* vec);
    double anisotropy_factor(const double& Er, const double* Fr);
    double coeff_rusanov(const int& i, const double* ErFr);
    void produit_cst(const double a, double* b, const int& size);
    std::vector<double> vector_sum(const double* vec1, const double* vec2, const int& size);
    double max_vec(const double* vec, const int& size);

    // Renvoie la donnée initiale
    void init(double* ErFr);

    // Calcul de l'erreur pour deux résultats donnés
    void Caculate_error(const double* vec_1, const double* vec_2, const int& nb_points_x,
                        const int& nb_points_y, const int& ind, const double& deltax,
                        const double& deltay);
    // Renvoie le vecteur error
    const std::vector<double>& Get_vec_error() const {return _vec_error;};

    // Retourne le vecteur _vec_coeff pour pouvoir calculer le pas de temps
    const std::vector<double>& Get_vec_coeff() const {return _vec_coeff;};

    // Fonction minmod
    double minmod(const double& x, const double& y);

    // Fonction minmod
    double van_leer(const double& x);
};

#define _FUNCTION_H
#endif
