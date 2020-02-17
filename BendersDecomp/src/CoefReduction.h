#ifndef _COEFREDUCTION_H_
#define  _COEFREDUCTION_H_
#include <cmath>
#include <stdlib.h>
#include "gurobi_c++.h"
#include <sstream>
#include <tuple>
using namespace std;

class CoefReduction{
private:
	int N;
	GRBModel* model;
	vector<vector<GRBVar>> w;
public:
	CoefReduction(GRBModel*, int);
	void set_objective(double, vector<tuple<int, int, double>> &);
	void set_constraints();
	void fix_edge(int, int);
	void free_edge(int, int);
	double solve();
	string itos(int i) { stringstream s; s << i; return s.str(); }
};


#endif // ! _COEFREDUCTION_H_
