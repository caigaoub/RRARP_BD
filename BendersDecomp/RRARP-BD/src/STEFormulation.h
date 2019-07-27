#ifndef _STEFORMULATION_H_
#define  _STEFORMULATION_H_

#include <sstream>
#include "gurobi_c++.h"
#include "PartitionScheme.h"
#include "BendersCuts.h"
//#include "CoefReduction.h"

using namespace std;

class STEFormulation {
private:
	PartitionScheme *					_partition = nullptr;
	DualFormulation *					_dual = nullptr;
	int									_size_var_y;
	GRBModel*							_model;
	GRBVar**							_var_y;
	GRBVar*								_var_v;

	int status;

	int num_Benders_cuts_const = 0;
	int num_subtour_cuts_const = 0;
	int num_user_cuts = 0;
public:
	STEFormulation() {};
	~STEFormulation() {};
	void build_model(GRBModel*, PartitionScheme*, DualFormulation *);
	
	pair<double, double> solve_LP_TSP();
	pair<double, double> solve_IP_TSP();
	void get_optimal_sol(double **);
	double add_USER_cuts(double**);
	bool add_SECs(double **);
	void check_cutting_point(int, double**, vector<int> &, int&, vector<int>&);
	void check_subcomponents(double**, vector<int>&, int&, vector<int>&);
	void print_num_Benders_cuts() { cout << "Benders' cut: " << num_Benders_cuts_const << endl; }
	void print_num_subtour_cuts() { cout << "Subtour' cut: " << num_subtour_cuts_const << endl; }
	int get_num_Benders_cuts() { return num_Benders_cuts_const; }
	int get_num_subtour_cuts() { return num_subtour_cuts_const; }
	string itos(int i) { stringstream s; s << i; return s.str(); }
	void set_model_MIP();
	void set_model_LP();
	void printSol(GRBModel*);
};





#endif
