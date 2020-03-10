#ifndef _PARTITIONSCHEME_H_
#define _PARTITIONSCHEME_H_
#include "DataHandler.h"
#include <tuple>
#include <cmath>
#include <math.h>
#include <deque>
class myVector {
public:
	double _x;
	double _y;
	myVector(Vertex initialp, Vertex endp) {
		this->_x = endp._x - initialp._x;
		this->_y = endp._y - initialp._y;
	}
	~myVector() {};
	inline double get_vecLen() {return sqrt(_x*_x + _y*_y); };
};


class PartitionScheme {
public:
	DataHandler *								_dataset = nullptr;
	int											_nb_dstzn = -1; // 
	double										_subarc_angle = -1;
	vector<vector<Vertex>>						_points; // all turning points after boundary partitioning
	int											_size_G;
	vector<vector<pair<bool, double>>>			_G; // network constructed after partitioning pair<travelable edge, risk value>
	int 										_nb_adm_OutT = 0;
	int 										_nb_adm_InT = 0;
	vector<vector<double>>						_min_risk_tars;// minimum risk between targets
	vector<double>								_MAX_REWARD_LIN;
	vector<double> 								_MAX_REWARD_ROT;
	double 										_weight_onInnRisk = 5.0;
	vector<double>								_par_c;  // params in risk & reward functions
	double										_par_h;
	vector<double>								_par_optOBdist; // mean of optimal obervation distance
	vector<double>								_par_varOBdist; // variance of optimal obs. distance

	PartitionScheme() {};
	~PartitionScheme() {};

	void build(DataHandler& instance, int nb_dstzn);
	void build_nodes_crds();
	void calc_MAX_REWARD();

	void get_risk_reward_linearInnerTrajc();
	double get_risk_linearInnerTrajc(Vertex, Vertex, int);
	double get_reward_linearInnerTrajc(Vertex, Vertex, int);
	void get_risk_outerTrajc();
	double get_risk_outerTrajc(Vertex, Vertex);
	tuple<bool,double,double> is_intersected(Vertex v, Vertex u, int i);
	// risk&reward-related functions
	double inner_risk_function(double l1, double l2, int i);
	double outer_risk_function(double dist);
	// supporting funcitons
	double dot_product(myVector vec1, myVector vec2);
	double get_lineSeg_len(Vertex, Vertex);
	void solve_shortestpath(vector<int> &, vector<vector<double>> &);
	void solve_shortestpath_v2(vector<int> &, vector<vector<double>> &);
	double calc_sequence_distance(vector<int> &);
};


#endif // !_PAR
