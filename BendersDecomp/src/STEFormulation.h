#ifndef _STEFORMULATION_H_
#define  _STEFORMULATION_H_

#include <sstream>
#include "gurobi_c++.h"
#include "PartitionScheme.h"
#include "BendersCuts.h"
#include "SuperCutFormulation.h"
#include "DualFormulation.h"
#include "GlobalMC.h"
#include "ProgTime.h"
#include <boost/filesystem.hpp>
#include <sys/stat.h>


using namespace std;

class STEFormulation {
public:
	PartitionScheme *					_partition = nullptr;
	DualFormulation *					_formul_dual = nullptr;
	SuperCutFormulation *				_formul_supercut = nullptr;
	int									_size_var_y;
	GRBModel*							_model;
	GRBVar**							_var_y;
	GRBVar*								_var_v;

	int 								_total_nb_Benders_cuts = 0;
	int									_total_nb_subtour_cuts = 0;
	int 								_total_nb_user_cuts = 0;
	

	int 								_optimstatus;

	ProgTime * 							_time;


	STEFormulation() {};
	~STEFormulation();
	void build_formul(GRBModel*, PartitionScheme*);
	
	void add_dualformul(DualFormulation *);
	void add_SuperCutformul(SuperCutFormulation*);
	pair<double, double> solve_formul_wCB(int);// with callback
	pair<double, double> solve_formul_woCB(); // without callback
	void get_optimal_sol(double **);
	double add_USER_cuts(double**);
	pair<bool,int> add_SECs(double **);
	void check_cutting_point(int, double**, vector<int> &, int&, vector<int>&);
	void check_subcomponents(double**, vector<int>&, int&, vector<int>&);
	string itos(int i) { stringstream s; s << i; return s.str(); }
	void set_vars_integer();
	void set_vars_continuous();
	void print_solution();
	void write_solution(string, int);
};





#endif
