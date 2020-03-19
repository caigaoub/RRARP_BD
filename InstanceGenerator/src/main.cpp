#include <iostream>
#include "gurobi_c++.h"
#include "InstanceGenerator.h"
#include <boost/filesystem.hpp>
#include "Reader.h"
#include "stats.cpp"

using namespace std;
void createInstance();
void createInstanceClustered();
int main() {
//  int num_targets = atoi(argv[1]);
//  char* difflevel = argv[2];

//  stats();
  createInstance();
  // createInstanceClustered();
  
  cout << "This is the end. Thank you!" << endl;
  return 0;

}


void createInstance(){
    // boost::filesystem::path cur_dir  = boost::filesystem::current_path();
    cout << "Instance Generator: " << endl;
    GRBEnv* evn = new GRBEnv();
    GRBModel model = GRBModel(*evn);
    string cur_dir = "/home/cai/Dropbox/Box_Research/Github/RRARP_BD/BendersDecomp/dat";
    for(int nb_t = 6; nb_t <= 35; nb_t++){
        cout << " ====>>> start creating instances with " << nb_t << " targets ~_~ " << endl;
        // boost::filesystem::path dir = cur_dir.string + "/dat/";
        // if(!boost::filesystem::exists(dir)){
        //     boost::filesystem::create_directories(dir);
        // }
        for(int j = 1; j <= 10; j++){
          cout << " ====>>> easy instance: " << j << endl;
          InstanceGenerator sample(nb_t, &model);
          sample.produce("e");
          string file = cur_dir + "/n_"+to_string(nb_t)+"_e_"+to_string(j)+".dat";
          sample.write_RRARP_instance(file);
          model.reset();
        }

        for(int j = 1; j <= 10; j++){
          cout << " ====>>> easy instance: " << j << endl;
          InstanceGenerator sample(nb_t, &model);
          sample.produce("m");
          string file = cur_dir + "/n_"+to_string(nb_t)+"_m_"+to_string(j)+".dat";
          sample.write_RRARP_instance(file);
          model.reset();
        }

        for(int j = 1; j <= 10; j++){
          cout << " ====>>> hard instance: " << j << endl;
          InstanceGenerator sample(nb_t, &model);
          sample.produce("h");
          string file = cur_dir + "/n_"+to_string(nb_t)+"_h_"+to_string(j)+".dat";
          sample.write_RRARP_instance(file);
          model.reset();
        }
    }
    delete  evn;
}

void createInstanceClustered(){
  cout << "--creating clustered instances " << endl;
  cout << "Instance Generator: " << endl;
  boost::filesystem::path cur_dir  = boost::filesystem::current_path();
  GRBEnv* evn = new GRBEnv();
  GRBModel model = GRBModel(*evn);

  for(int nb_t= 6; nb_t <=30; nb_t++){
    boost::filesystem::path dir = cur_dir.string() + "/ret/cluster_n_" + to_string(nb_t);
    if(!boost::filesystem::exists(dir)){
        boost::filesystem::create_directories(dir);
    }
    if (nb_t <= 12){
      int nb_cls = 2;
      for(int j = 1; j <= 10; j++){
        InstanceGenerator cluster(nb_t, &model);
        cluster.produce_clusters("m", nb_cls);
        string file = cur_dir.string() + "/ret/cluster_n_"+ to_string(nb_t)+ "/n_"+to_string(nb_t)+"_c_"+to_string(nb_cls) + "_" + to_string(j)+".txt";
        cluster.write_RRARP_cluster(file);
        model.reset();
      }
    }
    if (nb_t <= 16 && nb_t >12){
      int nb_cls = 3;
      for(int j = 1; j <= 10; j++){
        InstanceGenerator cluster(nb_t, &model);
        cluster.produce_clusters("m", nb_cls);
        string file = cur_dir.string() + "/ret/cluster_n_"+ to_string(nb_t)+ "/n_"+to_string(nb_t)+"_c_"+to_string(nb_cls) + "_" + to_string(j)+".txt";
        cluster.write_RRARP_cluster(file);
        model.reset();
      }
    }

    if (nb_t <= 20 && nb_t >16){
      int nb_cls = 4;
      for(int j = 1; j <= 10; j++){
        InstanceGenerator cluster(nb_t, &model);
        cluster.produce_clusters("m", nb_cls);
        string file = cur_dir.string() + "/ret/cluster_n_"+ to_string(nb_t)+ "/n_"+to_string(nb_t)+"_c_"+to_string(nb_cls) + "_" + to_string(j)+".txt";
        cluster.write_RRARP_cluster(file);
        model.reset();
      }
    }

    if (nb_t <= 30 && nb_t >20){
      int nb_cls = 5;
      for(int j = 1; j <= 10; j++){
        InstanceGenerator cluster(nb_t, &model);
        cluster.produce_clusters("m", nb_cls);
        string file = cur_dir.string() + "/ret/cluster_n_"+ to_string(nb_t)+ "/n_"+to_string(nb_t)+"_c_"+to_string(nb_cls) + "_" + to_string(j)+".txt";
        cluster.write_RRARP_cluster(file);
        model.reset();
      }
    }
  }
}
