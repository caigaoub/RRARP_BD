#include "Reader.h"


void Instance::read(string filename){
  _name = filename;
  auto pos = _name.find_last_of("/");
   _name = _name.substr(pos+1, _name.size());
  fstream file(filename);
  if (!file) {
   		cerr << "ERROR: could not open file '" << filename << "' for reading'" << endl;
   		throw(-1);
  }
  file >>  _nb_targets;
  file >> _depot1.x >> _depot1.y;
  file >> _depot2.x >> _depot2.y;
  _centers.resize(_nb_targets);
  _radii.resize(_nb_targets);
  for(int i=0; i < _nb_targets; i++){
    file >> _centers[i].x >> _centers[i].y >> _radii[i];
  }

}


void Instance::analyze(){
  // 1. average  radius length
  double total_dist = 0;
  for(int i=0; i < _nb_targets; i++){
    total_dist += _radii[i];
  }
  _avg_radius = total_dist/_nb_targets;
  // 2. average boundary distance
  total_dist = 0.0;
  int num_pairs = 0;
  total_dist += sqrt(pow(_depot1.x - _depot2.x, 2) + \
                pow(_depot1.y - _depot2.y, 2));
  num_pairs++;
  for(int i=0;i<_nb_targets;i++){
    total_dist += sqrt(pow(_centers[i].x - _depot1.x, 2) + \
                  pow(_centers[i].y - _depot1.y, 2)) - _radii[i];
    num_pairs++;
  }
  for(int i=0;i<_nb_targets;i++){
    total_dist += sqrt(pow(_centers[i].x - _depot2.x, 2) + \
                  pow(_centers[i].y - _depot2.y, 2)) - _radii[i];
    num_pairs++;
  }
  for(int i = 0; i < _nb_targets; i++){
    for(int j = 0; j < _nb_targets; j++){
      if (i != j){
        num_pairs++;
        total_dist += sqrt(pow(_centers[i].x - _centers[j].x, 2) + \
                      pow(_centers[i].y - _centers[j].y, 2)) - _radii[i] - _radii[j];
      }
    }
  }
  _avg_bdy_dist = total_dist / (double)num_pairs;
  // 3. average centerwise distance
  total_dist = 0.0;
  num_pairs = 0;
  total_dist += sqrt(pow(_depot1.x - _depot2.x, 2) + \
                pow(_depot1.y - _depot2.y, 2));
  num_pairs++;
  for(int i=0;i<_nb_targets;i++){
    total_dist += sqrt(pow(_centers[i].x - _depot1.x, 2) + \
                  pow(_centers[i].y - _depot1.y, 2));
    num_pairs++;
  }
  for(int i=0;i<_nb_targets;i++){
    total_dist += sqrt(pow(_centers[i].x - _depot2.x, 2) + \
                  pow(_centers[i].y - _depot2.y, 2));
    num_pairs++;
  }
  for(int i = 0; i < _nb_targets; i++){
    for(int j = 0; j < _nb_targets; j++){
      if (i != j){
        num_pairs++;
        total_dist += sqrt(pow(_centers[i].x - _centers[j].x, 2) + \
                      pow(_centers[i].y - _centers[j].y, 2));
      }
    }
  }
  _avg_ctr_dist = total_dist / (double)num_pairs;
}
