#ifndef _SUPERCUTFORMULATION_H_
#define  _SUPERCUTFORMULATION_H_
#include <cmath>
#include <stdlib.h>
#include "gurobi_c++.h"
#include <sstream>
#include <tuple>
using namespace std;

class SuperCutFormulation{
public:
	int 						_size_var_x = -1;
	GRBModel* 					_model;
	GRBVar **				 	_var_x;
	GRBVar * 					_var_a;
	GRBVar * 					_var_z;
	int 						_idx = 1;
	~SuperCutFormulation();
	void create_variables(GRBModel*, int);
	void set_objective(double, vector<pair<pair<int, int>, double>> &);
	void set_objective();	
	void update_coefs(double, vector<pair<pair<int, int>, double>> &);
	void set_constraints();
	void fix_edge(int, int);
	void free_edge(int, int);
	double solve();
	void print_solution(double **);
	string itos(int i) { stringstream s; s << i; return s.str(); }
};


#endif // ! _SUPERCUTFORMULATION_H_
