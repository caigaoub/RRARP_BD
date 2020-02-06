#ifndef _STEFORMULATION_H_
#define  _STEFORMULATION_H_

#include <sstream>
#include "gurobi_c++.h"
#include "PartitionScheme.h"
#include "BendersCuts.h"
//#include "CoefReduction.h"

using namespace std;

class STEFormulation {
public:
	PartitionScheme *					_partition = nullptr;
	// DualFormulation *					_formul_dual = nullptr;

	int									_size_var_y;
	GRBModel*							_model;
	GRBVar**							_var_y;
	GRBVar*								_var_v;

	int 								_total_nb_Benders_cuts = 0;
	int									_total_nb_subtour_cuts = 0;
	int 								_total_nb_user_cuts = 0;
	int 								_status;

	STEFormulation() {};
	~STEFormulation() {};
	void build_formul(GRBModel*, PartitionScheme*);
	
	pair<double, double> solve_LP_TSP();
	pair<double, double> solve_IP_TSP();
	void get_optimal_sol(double **);
	// double add_USER_cuts(double**);
	// bool add_SECs(double **);
	void check_cutting_point(int, double**, vector<int> &, int&, vector<int>&);
	// void check_subcomponents(double**, vector<int>&, int&, vector<int>&);
	string itos(int i) { stringstream s; s << i; return s.str(); }
	void set_model_MIP();
	void set_model_LP();
	void printSol(GRBModel*);
};





#endif
