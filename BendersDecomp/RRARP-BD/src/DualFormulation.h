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


private:
	GRBModel * model_dual; 
	double ** G;
	int num_dstzn;
	int num_targets;

	GRBVar* beta;
	GRBVar** alpha;
	
	int size_alpha;
	int size_beta;

	int status_dual;
public:
	DualFormulation(GRBModel*, PartitionScheme *, int );
	~DualFormulation();
	void set_objective(double ** );
	void set_constraints();
	void remove_constraints();
	double solve();
	void get_Benders_user_cut(GRBLinExpr &, GRBVar **);	
	string itos(int i) { stringstream s; s << i; return s.str(); }
};

#endif // ! _DUALFORMULATION_H_
