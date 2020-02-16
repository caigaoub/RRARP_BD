#ifndef  _DUALFORMULATION_H_
#define  _DUALFORMULATION_H_


#include "gurobi_c++.h"
#include "PartitionScheme.h"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

using namespace std;

class DualFormulation {
public:
	GRBModel * 												_model_dual; 
	vector<vector<pair<bool, double>>> * 					_G;
	int 													_nb_dstzn;
	int 													_nb_targets;

	GRBVar* 												_alpha;
	GRBVar** 												_beta;	
	int 													_size_alpha = 0;
	int 													_size_constrs = 0;
	int 													_status = -1;

	DualFormulation(){};
	~DualFormulation();
	void create_variables(GRBModel * model_dual, PartitionScheme * partition_);
	void set_objective(double ** );
	void set_constraints();
	void remove_constraints();
	double solve();
	void get_Benders_user_cut(GRBLinExpr &, GRBVar **);	
	string itos(int i) { stringstream s; s << i; return s.str(); }
};

#endif // ! _DUALFORMULATION_H_
