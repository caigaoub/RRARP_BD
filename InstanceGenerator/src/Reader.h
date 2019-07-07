#ifndef _READER_H_
#define _READER_H_

#include <string>
#include <math.h>
#include <algorithm>
#include <map>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <ctype.h>
#include <float.h>
#include<string>
#include <stdarg.h>
#include <limits.h>
#include <sstream>
#include <locale.h>
#include <errno.h>
#include <typeinfo>
#include <fstream>
#include "InstanceGenerator.h"

class Instance{
public:
  string                    _name;
  int                       _nb_targets;
  Vertex                    _depot1;
  Vertex                    _depot2;
  vector<Vertex>            _centers;
  vector<double>            _radii;

  double                    _avg_radius;
  double                    _avg_bdy_dist; // average boundary distance
  double                    _avg_ctr_dist; // average centerwise distance
  void   read(string filename);
  void analyze();
};

/*
class InstanceSet{

public:
  int                       _nb_instances;
  vector<Instance>         _instances;


  void get_average();
};

*/
#endif
