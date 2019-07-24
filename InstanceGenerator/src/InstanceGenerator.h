#ifndef _INSTANCEGENERATOR_H_
#define _INSTANCEGENERATOR_H_
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <string.h>
#include <vector>
#include <math.h>
#include <cstdlib> // for std::rand() and std::srand()
#include <ctime> // for std::time()
#include <random>
#include "gurobi_c++.h"
using namespace std;
struct Vertex{
  int x;
  int y;
};

class InstanceGenerator{
private:
    int              num_targets;
    GRBModel *       model;
    GRBVar*           r;

    random_device     _rd;//obtain a random number from hardware
    mt19937           _eng;

    Vertex bot_left_corner;
    int scale;
    int panel_width;
    int panel_height;
    int nb_cls;


    Vertex depot1_loc;
    Vertex depot2_loc;
    vector<Vertex> targets_locs;
    vector<double> max_radii;
    vector<double> radii;
    vector<double> min_reward_pct;
    vector<double> max_risk_pct;

public:
    InstanceGenerator(int, GRBModel*);
    ~InstanceGenerator(){};
    /* --- Primary functions --- */
    void set_panel(Vertex, double, double);
    void set_locations();
    void set_locations(int nb_cls);
    void get_max_radii();
    void set_radii(const char*);
    void set_RR_threshold();
    void produce(const char*);
    void produce_clusters(const char*, int);
    /*--- Correctness Check ---*/
    bool is_same_loc();
    bool is_same_loc(Vertex, Vertex);
    bool is_point_in_panel(Vertex);
    bool is_intersected(int, int);
    /* --- Others --- */
    double eucl_distance(Vertex, Vertex);
    /*--- Output ---*/
    void print_instance();
    void write_RRARP_instance(string);
    void write_RRARP_cluster(string); // generate clusters
    inline static string itos(int i) {stringstream s; s << i; return s.str();};
    inline static string dtos(double i) {stringstream s; s << i; return s.str();};
};

#endif
