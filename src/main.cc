#include <cstdio>
#include <string>
#include <cstdlib>
#include <unistd.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>
#include "DataFile.h"
#include "function.h"
#include "FiniteVolume.h"
#include "TimeScheme.h"
#include "save.h"

using namespace std;

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
      printf("Please, enter the name of your data file.\n");
      exit(0);
    }
    const string data_file_name = argv[1];

    // ----------------------- Fichier de données --------------------------------
    DataFile* p_df = new DataFile(data_file_name);
    p_df->Read_data_file();

    //---------------------------- Fonctions -------------------------------------
    Function* p_fct = new Function(p_df);

    //-------------------------- Volumes finis -----------------------------------
    FiniteVolume* p_fv = new FiniteVolume(p_df, p_fct);

    //-------------------------- Schémas en temps --------------------------------
    TimeScheme* p_ts = new TimeScheme(p_df, p_fct, p_fv);

    //-------------------------- Sauvegarde de la solution ------------------------
    Save* p_s = new Save(p_df, p_fct);

    // Initialisation de certaines variables
    int nb_simulation = p_df->Get_nb_simulation();
    double Tf = p_df->Get_Tf();
    std::string name;
    int nb_points_x, nb_points_y, Nx, Ny, dim;
    double dt, t, dx, Lx;
    std::vector<double> ErFr, ErFr2;
    Lx = p_df->Get_Lx();

    if (nb_simulation > 0) {
        if(dim == 1) {
            // Utile pour la convergence
            Nx = p_df->Get_Nx();
            nb_points_x = Nx;
            nb_points_x = 0;
        }
        else {
            Nx = p_df->Get_Nx();
            Ny = p_df->Get_Ny();
            nb_points_x = Nx;
            nb_points_y = Ny;
        }
        for(int k=0; k<nb_simulation; k++) {
            printf("k = %d\n", k);

            // Premier calcul
            // Variable qui change pour chaque pas de simulation
            t = 0.0;
            Nx = p_df->Get_Nx();
            Ny = p_df->Get_Ny();
            dx = Lx/(Nx+1);
            dim = p_df->Get_dim();

            if(dim == 1) {
                ErFr.resize(2*Nx);
            }
            else {
                ErFr.resize(3*Nx*Ny);
            }

            // Initialisation du cas
            p_fct->init(&ErFr[0]);

            // Boucle en temps
            while(t < Tf) {
                p_fv->Vector_Flux(&ErFr[0]);
                dt = p_ts->Get_dt();
                printf("dt = %e\n", dt);
                p_ts->Euler_Explicite(&ErFr[0], dt);
                t = t+dt;
            }

            name = "result_" + to_string(k) + ".txt";
            p_s->Save_Data(&ErFr[0], name);

            p_df -> change_Nx(2*Nx+1);
            p_df -> change_Ny(2*Ny+1);

            // Second calcul
            // Variable qui change pour chaque pas de simulation
            t = 0.0;
            Nx = p_df->Get_Nx();
            Ny = p_df->Get_Ny();
            dim = p_df->Get_dim();

            if(dim == 1) {
                ErFr2.resize(2*Nx);
            }
            else {
                ErFr2.resize(3*Nx*Ny);
            }

            // Initialisation du cas
            p_fct->init(&ErFr2[0]);

            // Boucle en temps
            while(t < Tf) {
                p_fv->Vector_Flux(&ErFr2[0]);
                dt = p_ts->Get_dt();
                printf("dt = %e\n", dt);
                p_ts->Euler_Explicite(&ErFr2[0], dt);
                t = t+dt;
            }

            p_fct->Caculate_error(&ErFr[0], &ErFr2[0], nb_points_x, nb_points_y, k, dx, 0.0);
        }

        std::vector<double> convergence(3*nb_simulation);
        convergence = p_fct->Get_vec_error();
        // Sauvegarde de la convergence
        p_s->Save_convergence(&convergence[0], nb_simulation);
    }

    else {
        // Sans convergence
        // Initialisation de certaines variables
        Tf = p_df->Get_Tf();
        t = 0.0;
        dt = 0.0;
        Nx = p_df->Get_Nx();
        Ny = p_df->Get_Ny();
        dim = p_df->Get_dim();
        if(dim == 1) {
            ErFr.resize(2*Nx);
        }
        else {
            ErFr.resize(3*Nx*Ny);
        }

        // Initialisation du cas
        p_fct->init(&ErFr[0]);

        // Boucle en temps
        while(t < Tf) {
            p_fv->Vector_Flux(&ErFr[0]);
            dt = p_ts->Get_dt();
            printf("dt = %e\n", dt);
            //p_ts->Euler_Explicite(&ErFr[0], dt);
            p_ts->RK2(&ErFr[0], dt);
            t = t+dt;
        }

        //Sauvegarde de la solution au temps Tf
        p_s->Save_Data(&ErFr[0], "result_Tf.txt");
    }

    return 0;
}
