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
	int 						_size_var_x;
	GRBModel* 					_model;
	GRBVar **				 	_var_x;

	void create_variables(GRBModel*, int);
	void set_objective(double, vector<tuple<int, int, double>> &);
	void set_constraints();
	void fix_edge(int, int);
	void free_edge(int, int);
	double solve();
	string itos(int i) { stringstream s; s << i; return s.str(); }
};


#endif // ! _SUPERCUTFORMULATION_H_
