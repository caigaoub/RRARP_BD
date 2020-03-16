#include "InstanceGenerator.h"

InstanceGenerator::InstanceGenerator(int var_num_targets, GRBModel *model){
  this->num_targets = var_num_targets;
  this->model = model;
  this->scale = 1;

  for(int i = 0;i < num_targets;i++){
    targets_locs.push_back({-1, -1});
    radii.push_back(-1);
    max_radii.push_back(-1);
  } 

  model->set(GRB_IntAttr_ModelSense, 1);
  model->set(GRB_IntParam_OutputFlag, 0);
  r = new GRBVar[num_targets];
  double UB = 10000000;
  for(int i = 0; i< num_targets; i++){
    r[i] = model->addVar(0.0, UB, 1.0, GRB_CONTINUOUS, "r_" + itos(i));
  }
  min_reward_pct.resize(num_targets);
  max_risk_pct.resize(num_targets);
  model ->update();
}


void InstanceGenerator::produce(const char* difflevel){
  /*solve the mathemtical model */
  // set_panel({0, 0}, num_targets* scale, num_targets * scale);
  // set_locations();
  // get_max_radii();
  // set_radii(difflevel);
  // set_RR_threshold();
//  print_instance();
  /* strategy: choose the nearest target and take half length as its maximum radius */

  set_panel({0, 0}, 100.0, 100.0);
  set_locations();
  get_max_radii_strategyII();
  set_radii_strategyII(difflevel);
  set_RR_threshold();

}


// void InstanceGenerator::produce_clusters(const char* difflevel, int nb_cls){
//   this->nb_cls = nb_cls;
//   set_panel({0, 0}, num_targets* scale, num_targets * scale);
//   set_locations(nb_cls);
//   get_max_radii();
//   set_radii(difflevel);
//   set_RR_threshold();
// //  print_instance();
// }


void InstanceGenerator::set_panel(Vertex llc, double width, double height){
  /* a square panel with for points:
    [panel_height,0]  -------[panel_width,panel_height]
        |                          |
        |                          |
      [0 0] ------------------- [panel_width, 0]
  */
  this->bot_left_corner = {llc.x, llc.y};
  this->panel_width = width;
  this->panel_height = height;
}

void InstanceGenerator::set_locations(){
  int crd_x, crd_y;
  while(true){ // centers are required to be within the panel.
  //  unsigned int time_ui = static_cast<unsigned int>( time(NULL)%1000 );
  //    srand( time_ui );
    _eng = mt19937(_rd());// seed the random generator
    auto randx = uniform_int_distribution<>(bot_left_corner.x, bot_left_corner.x + panel_width);
    auto randy = uniform_int_distribution<>(bot_left_corner.y, bot_left_corner.y + panel_height);
    for(int i = 0; i < num_targets + 2; i++){
      if (i == 0) {
        crd_x = randx(_eng);
        crd_y = randy(_eng);
        depot1_loc = {crd_x, crd_y};
      }
      if (i == 1) {
        crd_x = randx(_eng);
        crd_y = randy(_eng);
        depot2_loc = {crd_x, crd_y};
      }
      if (i > 1){
        crd_x = randx(_eng);
        crd_y = randy(_eng);
        targets_locs[i-2] = {crd_x, crd_y};
      }
    }
    if (!is_same_loc()){
      break;
    }
  }

}

void InstanceGenerator::set_locations(int nb_cls){
  double cl_dist = scale * num_targets/2.0;
  while(true){
    vector<Vertex> clu(nb_cls);
    _eng = mt19937(_rd());// seed the random generator
    auto randx = uniform_int_distribution<>(bot_left_corner.x, bot_left_corner.x + panel_width);
    auto randy = uniform_int_distribution<>(bot_left_corner.y, bot_left_corner.y + panel_height);
    while (true){ // generate a dummy center for each cluster
        bool ctn_flag = true;
        for(unsigned i =0; i < clu.size(); i++){
          clu[i] = {randx(_eng), randy(_eng)};
        }
        for(unsigned i=0; i< clu.size(); i++){
          for(unsigned j=0; j<clu.size(); j++){
            if (i != j){
                if (eucl_distance(clu[i], clu[j])<cl_dist){
                  ctn_flag = false;
                  break;
                }
            }
          }
          if (ctn_flag == false)
            break;
        }
        if (ctn_flag == true){
          break;
        }
    }

    vector<int> nb_tars_cl(nb_cls, num_targets/nb_cls);
    for(int i = 0; i < num_targets % nb_cls; i++){ // number of targets in each cluster
      nb_tars_cl[i] += 1;
    }
    // generate centers
    depot1_loc = {randx(_eng), randy(_eng)};
    depot2_loc = {randx(_eng), randy(_eng)};
    int k = 0;
    for(int i=0; i< nb_cls; i++){
      auto randx = uniform_int_distribution<>(clu[i].x - scale/4, clu[i].x + scale/4);
      auto randy = uniform_int_distribution<>(clu[i].y - scale/4, clu[i].y + scale/4);
      for(int j=0; j< nb_tars_cl[i]; j++){
        targets_locs[k] = {randx(_eng), randy(_eng)};
        k++;
      }
    }
    if (!is_same_loc()){
      break;
    }
  }

}


void InstanceGenerator::set_radii(const char* difflevel){
  _eng = mt19937(_rd());// seed the random generator
  double ratio;
  char strE[] = "e";
  char strM[] = "m";
  char strH[] = "h";
  if (strcmp(difflevel, strE)==0){ // easy
    for(int i = 0;i < num_targets; i++){
      auto rand_real = uniform_real_distribution<>(0.1, 0.3);
      ratio = rand_real(_eng);
      radii[i] = max_radii[i] * ratio;
    }

  }else if (strcmp(difflevel, strM)==0){ // medium
    for(int i = 0;i < num_targets; i++){
      auto rand_real = uniform_real_distribution<>(0.3, 0.5);
      ratio = rand_real(_eng);
      radii[i] = max_radii[i] * ratio;
    }

  }else if (strcmp(difflevel, strH)==0){ // hard
    for(int i = 0;i < num_targets; i++){
      auto rand_real = uniform_real_distribution<>(0.5, 0.7);
      ratio = rand_real(_eng);
      radii[i] = max_radii[i] * ratio;
    }
  }else{
    throw "Wrong difficulty level input";
  }
}


void InstanceGenerator::set_radii_strategyII(const char* difflevel){
  _eng = mt19937(_rd());// seed the random generator
  double ratio;
  char strE[] = "e";
  char strM[] = "m";
  char strH[] = "h";
  if (strcmp(difflevel, strE)==0){ // easy
    for(int i = 0;i < num_targets; i++){
      auto rand_real = uniform_real_distribution<>(0.2, 0.4);
      ratio = rand_real(_eng);
      radii[i] = max_radii[i] * ratio;
    }

  }else if (strcmp(difflevel, strM)==0){ // medium
    for(int i = 0;i < num_targets; i++){
      auto rand_real = uniform_real_distribution<>(0.4, 0.6);
      ratio = rand_real(_eng);
      radii[i] = max_radii[i] * ratio;
    }

  }else if (strcmp(difflevel, strH)==0){ // hard
    for(int i = 0;i < num_targets; i++){
      auto rand_real = uniform_real_distribution<>(0.6, 0.8);
      ratio = rand_real(_eng);
      radii[i] = max_radii[i] * ratio;
    }
  }else{
    throw "Wrong difficulty level input";
  }
}


void InstanceGenerator::set_RR_threshold(){
  _eng = mt19937(_rd());// seed the random generator
  // double ratio = 0.0;
  for(int i = 0; i < num_targets; i++){
    auto rand_real_reward = uniform_real_distribution<>(0.3, 0.5);
    ratio = rand_real_reward(_eng);
    min_reward_pct[i] = 0.ratio;
    // auto rand_real_risk = uniform_real_distribution<>(0.4, 0.6);
    // ratio = rand_real_risk(_eng);
    // max_risk_pct[i] = 100000;
  }
}


void InstanceGenerator::get_max_radii(){
  /* Objective */
  GRBLinExpr expr_obj = 0;
  for (int i = 0; i < num_targets; i++){
    expr_obj += r[i];
  }
  model->setObjective(expr_obj, GRB_MAXIMIZE);
  /* constraints*/
  double smallest_dist = numeric_limits<int>::max();
  double len;
  for(int i =0; i< num_targets; i++){
    len = eucl_distance(depot1_loc, targets_locs[i]);
    cout << "d1" << "-" << i+1 << " : "<< len << endl;
    model->addConstr(r[i] <= len, "C_d1" + itos(i+1));
    if (len < smallest_dist){
      smallest_dist = len;
    }
  }
  for(int i =0; i< num_targets; i++){
    len = eucl_distance(depot2_loc, targets_locs[i]);
    cout << "d2" << "-" << i+1 << " : "<< len << endl;
    model->addConstr(r[i] <= len, "C_d2" + itos(i+1));
    if (len < smallest_dist){
      smallest_dist = len;
    }
  }
  for(int i = 0; i < num_targets; i++){
    for(int j = 0; j < i; j++){
      len = eucl_distance(targets_locs[i], targets_locs[j]);
     cout << i+1 << "-" << j+1 << " : "<< len << endl;
      model->addConstr(r[i] + r[j] <= len, "C_" + itos(i+1) + itos(j+1));
      if (len < smallest_dist){
        smallest_dist = len;
      }
    }
  }
  model->update();
//  double r_LB = min((double)num_targets, smallest_dist/2.0);
  double r_LB = smallest_dist/2.0;
  for(int i = 0; i < num_targets; i++){
    model->addConstr( r[i] >= r_LB, "Clb_" + itos(i));
  }
  model->update();
  try	{
    model->optimize();
		if (model->get(GRB_IntAttr_Status) == GRB_OPTIMAL){
      for(int i = 0; i < num_targets; i++){
        max_radii[i] = r[i].get(GRB_DoubleAttr_X);
        cout << max_radii[i] << "  ";
      }
    }
//    cout << endl;
	}
	catch (GRBException e) {
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during optimization" << endl;
	}
}



void InstanceGenerator::get_max_radii_strategyII(){
  double smallest_dist = numeric_limits<int>::max();
  double len = numeric_limits<int>::max();
  for(int i = 0; i < num_targets; i++){
    smallest_dist = numeric_limits<int>::max();
      len = eucl_distance(targets_locs[i], depot1_loc);
      if (smallest_dist > len){
        smallest_dist = len;
      }
      len = eucl_distance(targets_locs[i], depot2_loc);
      if (smallest_dist > len){
        smallest_dist = len;
      }

    for(int j = 0; j < num_targets; j++){
      if(i != j){
        len = eucl_distance(targets_locs[i], targets_locs[j]);
        if (smallest_dist > len){
          smallest_dist = len;
        }
      }
    }
    max_radii[i] = 0.5 * smallest_dist;
  }
}


double InstanceGenerator::eucl_distance(Vertex p, Vertex q){
  return sqrt(pow(p.x - q.x, 2) + pow(p.y - q.y, 2));
}

bool InstanceGenerator::is_same_loc(){
  bool flag= false;
  if(is_same_loc(depot1_loc, depot2_loc)){
    flag = true;
  }
  for(int i=0; i<num_targets;i++){
    if(is_same_loc(depot1_loc, targets_locs[i])){
      flag = true;
      break;
    }
  }
  for(int i=0; i<num_targets;i++){
    if(is_same_loc(depot2_loc, targets_locs[i])){
      flag = true;
      break;
    }
  }
  for(int i=0; i<num_targets;i++){
    for(int j=0;j<i;j++){
      if(is_same_loc(targets_locs[i], targets_locs[j])){
        flag = true;
        break;
      }
    }
  }
  return flag;
}


bool InstanceGenerator::is_same_loc(Vertex p, Vertex q){
  if (p.x == q.x && p.y == q.y){
    return true;
  }else{
    return false;
  }
}

bool InstanceGenerator::is_intersected(int i, int j){
  double len = eucl_distance(targets_locs[i], targets_locs[j]);
  if(len - radii[i] - radii[j] <= 0){
    return true;
  }else{
    return false;
  }
}

bool InstanceGenerator::is_point_in_panel(Vertex p){
  if (p.x >= bot_left_corner.x && \
      p.x <= bot_left_corner.x + panel_width && \
      p.y >= bot_left_corner.y && \
      p.y <= bot_left_corner.y + panel_height){
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

}

void InstanceGenerator::write_RRARP_instance(string path){
  ofstream myfile;
  myfile.open(path);
  // cout << path << endl;
  myfile << num_targets << '\n';
  myfile <<  depot1_loc.x << '\t' << depot1_loc.y << '\n';
  myfile <<  depot2_loc.x << '\t' << depot2_loc.y << '\n';
  for(int i =0 ;i<num_targets;i++){
      myfile << targets_locs[i].x << '\t' << targets_locs[i].y << '\t' << radii[i] << '\n';
  }
  for(int i = 0; i<num_targets;i++){
    if(i != num_targets-1){
      myfile << min_reward_pct[i] << '\t';
    }else{
      myfile << min_reward_pct[i] << '\n';
    }
  }
  for(int i = 0;i<num_targets;i++){
    if(i != num_targets-1){
      myfile << max_risk_pct[i] << '\t';
    }else{
      myfile << max_risk_pct[i] << '\n';
    }
  }
  myfile.close();
}

void InstanceGenerator::write_RRARP_cluster(string path){
  ofstream myfile;
  myfile.open(path);
  myfile << num_targets << '\t' << nb_cls << '\n';
  myfile <<  depot1_loc.x << '\t' << depot1_loc.y << '\n';
  myfile <<  depot2_loc.x << '\t' << depot2_loc.y << '\n';
  for(int i =0 ;i<num_targets;i++){
      myfile << targets_locs[i].x << '\t' << targets_locs[i].y << '\t' << radii[i] << '\n';
  }
  for(int i =0 ;i<num_targets;i++){
    if(i != num_targets - 1){
      myfile << min_reward_pct[i] << '\t';
    }else{
      myfile << min_reward_pct[i] << '\n';
    }
  }
  for(int i =0 ;i<num_targets;i++){
    if(i != num_targets - 1){
      myfile << max_risk_pct[i] << '\t';
    }else{
      myfile << max_risk_pct[i] << '\n';
    }
  }
  myfile.close();
}
