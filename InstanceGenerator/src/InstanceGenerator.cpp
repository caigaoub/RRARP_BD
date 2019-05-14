#include "InstanceGenerator.h"
#include <cstdlib> // for std::rand() and std::srand()
#include <ctime> // for std::time()
#include <fstream>
InstanceGenerator::InstanceGenerator(int var_num_targets){
  this->num_targets = var_num_targets;
  this->scale = 80;

  for(int i = 0;i < num_targets;i++){
    targets_locs.push_back({-1, -1});
    radii.push_back(-1);
  }
}

void InstanceGenerator::produce(int difflevel){
  set_panel({0, 0});
  int num_insts = 0;
  int total_num_insts = 1;
  while(true){
    set_locations();
    set_radii(difflevel);
    if(check_correctness()){
      /*write the instance*/

      write_RRARP_instance();
      write_TSP_instance();
      get_aver_bry_dist();
      print_instance();
      num_insts++;
    }
    if (num_insts >= total_num_insts){
      break;
    }
  }
}

void InstanceGenerator::set_panel(Vertex llc){
  /* a square panel with for points:
    [panel_height,0]  -------[panel_width,panel_height]
        |                          |
        |                          |
      [0 0] ------------------- [panel_width, 0]
  */
  this->lower_left_corner = {llc.x, llc.y};
  this->panel_width = scale * num_targets;
  this->panel_height = scale * num_targets;
}

void InstanceGenerator::set_locations(){
  unsigned int time_ui = static_cast<unsigned int>( time(NULL)%1000 );
  srand( time_ui );
  int crd_x, crd_y;
  for(int i = 0; i < num_targets + 2; i++){
    if (i == 0) {
      crd_x = rand()%(panel_width + 1) + lower_left_corner.x;
      crd_y = rand()%(panel_height + 1) + lower_left_corner.y;
      depot1_loc = {crd_x, crd_y};
    }
    if (i == 1) {
      crd_x = rand()%(panel_width + 1) + lower_left_corner.x;
      crd_y = rand()%(panel_height + 1) + lower_left_corner.y;
      depot2_loc = {crd_x, crd_y};
    }
    if (i > 1){
      crd_x = rand()%(panel_width - 1) + lower_left_corner.x + 1;
      crd_y = rand()%(panel_height - 1) + lower_left_corner.y + 1;
      targets_locs[i-2] = {crd_x, crd_y};
    }
  }
}

void InstanceGenerator::set_radii(int difflevel){
  unsigned int time_ui = static_cast<unsigned int>( time(NULL)%1000 );
  srand( time_ui );
  if (difflevel == 1){ // easy
    for(int i=0;i<num_targets;i++){
      radii[i] = ((double)rand()/(RAND_MAX) * 15) + 1;
    }
  }else if (difflevel == 2){ // medium
    for(int i=0;i<num_targets;i++){
      radii[i] = ((double)rand()/(RAND_MAX) * 15) + 15;
    }
  }else if (difflevel == 3){ // hard
    for(int i=0;i<num_targets;i++){
      radii[i] = ((double)rand()/(RAND_MAX) * 15) + 50;
    }
  }else{
    throw "Wrong diff-level input";
  }
}


bool InstanceGenerator::check_correctness(){
  /* check whether all discs stay within the panel */
  for(int i = 0; i < num_targets; i++){
    if(!is_disc_in_panel(i)){
//      cout << "A disc not in the panel" << endl;
      return false;
    }
  }
  /* check whether there exists two intersected discs */
  for(int i = 0; i < num_targets; i++){
    for(int j = 0; j < i; j++){
      if(is_intersected(i,j)){
//          cout << "two discs are intersected " << endl;
        return false;
      }
    }
  }
  return true;
}


bool InstanceGenerator::is_intersected(int i, int j){
  double len = sqrt(pow(targets_locs[i].x - targets_locs[j].x, 2) + \
                    pow(targets_locs[i].y - targets_locs[j].y, 2));
  if(len - radii[i] - radii[j] <= 0){
    return true;
  }else{
    return false;
  }
}

bool InstanceGenerator::is_disc_in_panel(int tar){
  Vertex p;
  p.x = targets_locs[tar].x + radii[tar];
  p.y = targets_locs[tar].y;
  if (!is_point_in_panel(p)){
    return false;
  }
  p.x = targets_locs[tar].x;
  p.y = targets_locs[tar].y + radii[tar];
  if (!is_point_in_panel(p)){
    return false;
  }
  p.x = targets_locs[tar].x - radii[tar];
  p.y = targets_locs[tar].y;
  if (!is_point_in_panel(p)){
    return false;
  }
  p.x = targets_locs[tar].x;
  p.y = targets_locs[tar].y  - radii[tar];
  if (!is_point_in_panel(p)){
    return false;
  }
  return true;
}

bool InstanceGenerator::is_point_in_panel(Vertex p){
  if (p.x >= lower_left_corner.x && \
      p.x <= lower_left_corner.x + panel_width && \
      p.y >= lower_left_corner.y && \
      p.y <= lower_left_corner.y + panel_height){
    return true;
  }else{
    return false;
  }
}

void InstanceGenerator::print_instance(){
  cout << "depot 1  " << depot1_loc.x << " " << depot1_loc.y << endl;
  cout << "depot 2  " << depot2_loc.x << " " << depot2_loc.y << endl;
  for(int i =0 ;i<num_targets;i++){
      cout << "target " << i+1 << " " << targets_locs[i].x << " " << \
          targets_locs[i].y << " radius " << radii[i] << endl;
  }
  cout << "average boundary distance: " << aver_bry_dist << endl;

}

void InstanceGenerator::write_RRARP_instance(){
  ofstream myfile;
  myfile.open("RRARP_instance.dat");
  myfile << num_targets << '\n';
  myfile <<  depot1_loc.x << '\t' << depot1_loc.y << '\n';
  myfile <<  depot2_loc.x << '\t' << depot2_loc.y << '\n';
  for(int i =0 ;i<num_targets;i++){
      myfile << radii[i]<< '\t';
  }
  myfile << '\n';
  for(int i =0 ;i<num_targets;i++){
      myfile << targets_locs[i].x << '\t' << targets_locs[i].y << '\n';
  }
  myfile.close();
}

void InstanceGenerator::write_TSP_instance(){
  ofstream of;
  of.open("TSP_instance.dat");
  of << num_targets << '\n';
  of <<  depot1_loc.x << '\t' << depot1_loc.y << '\n';
  for(int i =0 ;i<num_targets;i++){
      of << targets_locs[i].x << '\t' << targets_locs[i].y << '\n';
  }
  of <<  depot2_loc.x << '\t' << depot2_loc.y << '\n';
  of.close();
}

double InstanceGenerator::get_aver_bry_dist(){
  /*average boundary distance */
  double total_dist = 0.0;
  int num_pairs = 0;
  total_dist += sqrt(pow(depot1_loc.x - depot2_loc.x, 2) + \
                pow(depot1_loc.y - depot2_loc.y, 2));
  num_pairs++;
  for(int i=0;i<num_targets;i++){
    total_dist += sqrt(pow(targets_locs[i].x - depot1_loc.x, 2) + \
                  pow(targets_locs[i].y - depot1_loc.y, 2)) - radii[i];
    num_pairs++;
  }
  for(int i=0;i<num_targets;i++){
    total_dist += sqrt(pow(targets_locs[i].x - depot2_loc.x, 2) + \
                  pow(targets_locs[i].y - depot2_loc.y, 2)) - radii[i];
    num_pairs++;
  }
  for(int i = 0; i < num_targets; i++){
    for(int j = 0; j < num_targets; j++){
      if (i != j){
        num_pairs++;
        total_dist += sqrt(pow(targets_locs[i].x - targets_locs[j].x, 2) + \
                      pow(targets_locs[i].y - targets_locs[j].y, 2)) - radii[i] - radii[j];
      }
    }
  }

  aver_bry_dist = total_dist / (double)num_pairs;
  return aver_bry_dist;
}
