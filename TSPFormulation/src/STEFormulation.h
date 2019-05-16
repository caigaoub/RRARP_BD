#ifndef _STEFORMULATION_H_
#define _STEFORMULATION_H_
#include "gurobi_c++.h"
#include "CostInfo.h"
#include "SubtourElimCuts.h"
using namespace std;
class STEFormulation
{
private:
	int N;
	GRBVar** y;
	CostInfo * cost;
	int numSubtourConst;
	unsigned int status;
public:
	STEFormulation(CostInfo * c, GRBModel * model);
	~STEFormulation();
	void solve(GRBModel * model, vector<int> &);
//	void writeSol(GRBModel * model);
	void printSol(GRBModel* model);
	string itos(int i) { stringstream s; s << i; return s.str();}
};
#endif // !_STEFORMULATION_H_
