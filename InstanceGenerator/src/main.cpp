#include <iostream>
#include "gurobi_c++.h"
#include "InstanceGenerator.h"
#include <boost/filesystem.hpp>


using namespace std;

int main(int argc, char* argv[]) {
//  int num_targets = atoi(argv[1]);
//  char* difflevel = argv[2];
  argc = argc;
  cout << "Instance Generator: " << endl;
  GRBEnv* evn = new GRBEnv();
  GRBModel model = GRBModel(*evn);
  for(int nb_t = 6; nb_t <= 20; nb_t++){
      cout << "--creating instace " << nb_t << endl;
      InstanceGenerator sample(nb_t, &model);
      boost::filesystem::path dir = "/media/caigao/LENOVO/ROTK/RRARP_BD/InstanceGenerator/ret/inst_n_"+to_string(nb_t);
      if(!boost::filesystem::exists(dir)){
          boost::filesystem::create_directories(dir);
      }
      for(int j = 1; j <= 10; j++){
        InstanceGenerator sample(nb_t, &model);
        sample.produce("e");
        string file = "/media/caigao/LENOVO/ROTK/RRARP_BD/InstanceGenerator/ret/inst_n_"+ \
                       to_string(nb_t)+ "/n_"+to_string(nb_t)+"_e_"+to_string(j)+".txt";
        sample.write_RRARP_instance(file);
        model.reset();
      }
      for(int j = 1; j <= 10; j++){
        InstanceGenerator sample(nb_t, &model);
        sample.produce("e");
        string file = "/media/caigao/LENOVO/ROTK/RRARP_BD/InstanceGenerator/ret/inst_n_"+ \
                      to_string(nb_t)+ "/n_"+to_string(nb_t)+"_m_"+to_string(j)+".txt";
        sample.write_RRARP_instance(file);
        model.reset();
      }
      for(int j = 1; j <= 10; j++){
        InstanceGenerator sample(nb_t, &model);
        sample.produce("e");
        string file = "/media/caigao/LENOVO/ROTK/RRARP_BD/InstanceGenerator/ret/inst_n_"+ \
                        to_string(nb_t)+ "/n_"+to_string(nb_t)+"_h_"+to_string(j)+".txt";
        sample.write_RRARP_instance(file);
        model.reset();
      }
  }

  for(int nb_t = 21; nb_t <= 30; nb_t++){
      cout << "--creating instace " << nb_t << endl;
      InstanceGenerator sample(nb_t, &model);
      boost::filesystem::path dir = "/media/caigao/LENOVO/ROTK/RRARP_BD/InstanceGenerator/ret/inst_n_"+to_string(nb_t);
      if(!boost::filesystem::exists(dir)){
          boost::filesystem::create_directories(dir);
      }
      for(int j = 1; j <= 10; j++){
        InstanceGenerator sample(nb_t, &model);
        sample.produce("e");
        string file = "/media/caigao/LENOVO/ROTK/RRARP_BD/InstanceGenerator/ret/inst_n_"+ \
                       to_string(nb_t)+ "/n_"+to_string(nb_t)+"_e_"+to_string(j)+".txt";
        sample.write_RRARP_instance(file);
        model.reset();
      }
  }

//  sample.write_RRARP_instance(dir);

  cout << "This is the end. Thank you!" << endl;
  return 0;

}
