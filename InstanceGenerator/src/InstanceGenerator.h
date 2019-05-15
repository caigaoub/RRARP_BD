#ifndef _INSTANCEGENERATOR_H_
#define _INSTANCEGENERATOR_H_
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <math.h>
#include "gurobi_c++.h"
using namespace std;
struct Vertex{
  int x;
  int y;
};

class InstanceGenerator{
private:
    int num_targets;
    GRBModel * model;
    GRBVar* r;

    Vertex lower_left_corner;
    int scale;
    int panel_width;
    int panel_height;

    Vertex depot1_loc;
    Vertex depot2_loc;
    vector<Vertex> targets_locs;
    vector<double> max_radii;
    vector<double> radii;
    double aver_bry_dist; // average boundary distance
public:
    InstanceGenerator(int, GRBModel*);
    ~InstanceGenerator(){};
    void set_panel(Vertex);
    void set_locations();
    void set_radii(int);
    void get_max_radii();

    bool is_same_loc();
    bool is_same_loc(Vertex, Vertex);
    bool is_in_panel();
    bool is_disc_in_panel(int);
    bool is_point_in_panel(Vertex);
    bool is_intersected(int, int);

    void produce(int);
    double eucl_distance(Vertex, Vertex);
    double get_aver_bry_dist();


    void print_instance();
    void write_RRARP_instance();
    void write_TSP_instance();



    inline static string itos(int i) {stringstream s; s << i; return s.str();};
    inline static string dtos(double i) {stringstream s; s << i; return s.str();};
};

#endif
