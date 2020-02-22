#ifndef _BENDERSCUTS_H_
#define _BENDERSCUTS_H_

#include "GlobalMC.h"
#include <cmath>
#include <stdlib.h>
#include <ctime>
#include <chrono>
#include <queue>
#include <sstream>
#include "gurobi_c++.h"
#include "PartitionScheme.h"
#include "DualFormulation.h"
#include "SuperCutFormulation.h"
using namespace std;

class BendersCuts : public GRBCallback {
private:
  	int 										_nb_dstzn;
	int 										_nb_targets;
	PartitionScheme* 							_partition = nullptr;
	DualFormulation* 							_formul_dual = nullptr;
	GRBModel *									_model_supercut = nullptr;
	SuperCutFormulation*						_formul_supercut = nullptr;

	vector<vector<pair<bool, double>>> * 		_G = nullptr; // network G=(V,E)
	int 										_which_Bcut;
	
	// GRBModel * 					_model_tsp;
	GRBVar ** 					_var_y;
	GRBVar *					_var_v;

	// vector<vector<double>>  	_SDS;
	vector<int> 				_fseq;
	vector<vector<int>> 		_SeqPool;
	int 						_max_supercuts = 10;
	int							_idx_supercut = 0;

	int 						_CB_nb_Benders_cuts=0;
	int 						_CB_nb_subtour_cuts=0;


public:
	~BendersCuts();
	// BendersCuts(GRBVar**, GRBVar*, PartitionScheme*);
	BendersCuts(GRBVar**, GRBVar*, PartitionScheme*, DualFormulation *, SuperCutFormulation*, int);
	static void findsubtour(int, double**, int*, int* );
	inline int get_nb_Benders_cuts() { return _CB_nb_Benders_cuts; }
	inline int get_nb_subtour_cuts() { return _CB_nb_subtour_cuts; }
	GRBLinExpr generate_Benderscut_SP(vector<int> *, vector<vector<double>> & );
  	GRBLinExpr generate_StrongBenderscut(vector<int> *, vector<vector<double>> & , bool);
	inline string itos(int i) { stringstream s; s << i; return s.str(); }
	void print_ySol(double**);
	// double improve_coef(int, int, double, vector<tuple<int, int, double>> &);
protected:
	void callback();
};


#endif // !_BENDERSCUTS_H_
