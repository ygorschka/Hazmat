#ifndef _DATA_FILE_H

#include <string>

class DataFile {
private:
    std::string _file_name, _case;
    double _a, _c, _Lx, _Ly, _coeff_rusanov, _Tf, _cfl, _x0, _y0;
    int _Nx, _Ny, _order_space, _order_time, _dim, _nb_simulation;

public:
    // Constructeur de la classe
    DataFile(std::string file_name);

    // Renvoie certaines variables
    const double Get_c() const {return _c;};
    const double Get_a() const {return _a;};
    const int Get_Nx() const {return _Nx;};
    const int Get_Ny() const {return _Ny;};
    const double Get_Lx() const {return _Lx;};
    const double Get_Ly() const {return _Ly;};
    const int Get_order_space() const {return _order_space;};
    const int Get_order_time() const {return _order_time;};
    const int Get_dim() const {return _dim;};
    const double Get_Tf() const {return _Tf;};
    const double Get_cfl() const {return _cfl;};
    const std::string Get_case() const {return _case;};
    const double Get_x0() const {return _x0;};
    const double Get_y0() const {return _y0;};
    const int Get_nb_simulation() const {return _nb_simulation;};

    // Lit le fichier d'entrée
    void Read_data_file();

    // Change la valeur de Nx
    void change_Nx(const int& new_nx);
    // Change la valeur de Ny
    void change_Ny(const int& new_ny);


};

#define _DATA_FILE_H
#endif
