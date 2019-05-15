#ifndef _BENDERSCUTS_H_
#define _BENDERSCUTS_H_

#include "GlobalMC.h"
#include <cmath>
#include <stdlib.h>
#include <ctime>
#include <chrono>
#include <queue>
#include "gurobi_c++.h"
#include "PartitionScheme.h"
#include "DualFormulation.h"
#include "CoefReduction.h"
using namespace std;

class BendersCuts : public GRBCallback {


private:

  int N;
  int num_dstzn;
	int num_targets;
	PartitionScheme* PS;
	DualFormulation* DL;
	double** G; // network G=(V,E)

	GRBVar ** y;
	GRBVar *v;
	GRBModel * model_tsp;

	vector<vector<double>> * SDS;
	vector<int> *fseq;
	vector<vector<int>> SeqPool;

	GRBLinExpr expr;
//	GRBEnv * evn_CoefRedc;
//	GRBModel * model_CoefRedc;
//	CoefReduction * CR;

	unsigned int num_Benders_cuts;
	unsigned int num_subtour_cuts;

	unsigned int flag;
	int count;

public:
	BendersCuts(GRBModel*, GRBVar**, GRBVar*, PartitionScheme*, DualFormulation *);
	static void findsubtour(int, double**, int*, int* );
	void check_cutting_point(int, double**, vector<int> &, int&, vector<int>&);
	void check_subcomponents(double**, vector<int>&, int&, vector<int>&);
	bool is_inSeqPool(vector<int> &);
	inline int get_num_Benders_cuts() { return num_Benders_cuts; }
	inline int get_num_subtour_cuts() { return num_subtour_cuts; }
	GRBLinExpr generate_Benderscut_SP(vector<int> *);
	string itos(int i) { stringstream s; s << i; return s.str(); }

	double improve_coef(int, int, double, vector<tuple<int, int, double>> &);
protected:
	void callback();
};


#endif // !_BENDERSCUTS_H_
