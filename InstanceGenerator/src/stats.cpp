#include <iostream>
#include "gurobi_c++.h"
#include "InstanceGenerator.h"
#include <boost/filesystem.hpp>
#include "Reader.h"
using namespace std;
void stats() {
    ofstream myfile;
    double avg_radius = 0;
    double avg_bdy_dist = 0;
    double avg_ctr_dist = 0;
    myfile.open("/media/caigao/LENOVO/ROTK/RRARP_BD/InstanceGenerator/ret/instancesstat.txt");
    for(int nb_t = 6; nb_t <= 20; nb_t++){
        Instance pe; // --------easy
        for(int j = 1; j <= 10; j++){
          string file = "/media/caigao/LENOVO/ROTK/RRARP_BD/InstanceGenerator/ret/inst_n_"+ \
                         to_string(nb_t)+ "/n_"+to_string(nb_t)+"_e_"+to_string(j)+".txt";
          pe.read(file);
          pe.analyze();
          avg_radius += pe._avg_radius;
          avg_bdy_dist += pe._avg_bdy_dist;
          avg_ctr_dist += pe._avg_ctr_dist;
        }
        avg_radius /= 10;
        avg_bdy_dist /= 10;
        avg_ctr_dist /= 10;
        myfile << nb_t << '\t' << "e" << '\t' << avg_ctr_dist << '\t' << avg_bdy_dist << '\t' << avg_radius << endl;

        Instance pm; // -------medium
        avg_radius = 0;
        avg_bdy_dist = 0;
        avg_ctr_dist = 0;
          for(int j = 1; j <= 10; j++){
          string file = "/media/caigao/LENOVO/ROTK/RRARP_BD/InstanceGenerator/ret/inst_n_"+ \
                         to_string(nb_t)+ "/n_"+to_string(nb_t)+"_m_"+to_string(j)+".txt";
          pm.read(file);
          pm.analyze();
          avg_radius += pm._avg_radius;
          avg_bdy_dist += pm._avg_bdy_dist;
          avg_ctr_dist += pm._avg_ctr_dist;
        }
        avg_radius /= 10;
        avg_bdy_dist /= 10;
        avg_ctr_dist /= 10;
        myfile << nb_t << '\t' << "m" << '\t' << avg_ctr_dist << '\t' << avg_bdy_dist << '\t' << avg_radius << endl;

        Instance ph; // // -------hard
        avg_radius = 0;
        avg_bdy_dist = 0;
        avg_ctr_dist = 0;
          for(int j = 1; j <= 10; j++){
          string file = "/media/caigao/LENOVO/ROTK/RRARP_BD/InstanceGenerator/ret/inst_n_"+ \
                         to_string(nb_t)+ "/n_"+to_string(nb_t)+"_h_"+to_string(j)+".txt";
          ph.read(file);
          ph.analyze();
          avg_radius += ph._avg_radius;
          avg_bdy_dist += ph._avg_bdy_dist;
          avg_ctr_dist += ph._avg_ctr_dist;
        }
        avg_radius /= 10;
        avg_bdy_dist /= 10;
        avg_ctr_dist /= 10;
        myfile << nb_t << '\t' << "h" << '\t' << avg_ctr_dist << '\t' << avg_bdy_dist << '\t' << avg_radius << endl;
    }

    for(int nb_t = 21; nb_t <= 30; nb_t++){
      Instance pee;
      avg_radius = 0;
      avg_bdy_dist = 0;
      avg_ctr_dist = 0;
      for(int j = 1; j <= 10; j++){
        string file = "/media/caigao/LENOVO/ROTK/RRARP_BD/InstanceGenerator/ret/inst_n_"+ \
                       to_string(nb_t)+ "/n_"+to_string(nb_t)+"_e_"+to_string(j)+".txt";
        pee.read(file);
        pee.analyze();
        avg_radius += pee._avg_radius;
        avg_bdy_dist += pee._avg_bdy_dist;
        avg_ctr_dist += pee._avg_ctr_dist;
      }
      avg_radius /= 10;
      avg_bdy_dist /= 10;
      avg_ctr_dist /= 10;
      myfile << nb_t << '\t' << "e" << '\t' << avg_ctr_dist << '\t' << avg_bdy_dist << '\t' << avg_radius << endl;
    }

}
