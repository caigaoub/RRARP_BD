#ifndef _PARTITIONSCHEME_H_
#define _PARTITIONSCHEME_H_
#include "DataHandler.h"
#include "myNameClass.h"
#include <tuple>

class PartitionScheme {
private:
	// members for receiving geographical info
	int num_targets;
	Vertex  depot1_loc;
	Vertex  depot2_loc;
	vector<Vertex> target_locs;
	vector<double> radii;
	vector<double> bdg_rewards;
	vector<double> risk_thold;

	// discretization info
	int num_dstzn;
	double subarc_angle;

	// graph G = (V,E)
	vector<vector<Vertex>> TP;
	int num_V;
	double** G;	
	double** min_risk_mat;// minimum risk for every pair of boundaries

	vector<double> par_c1;
	double par_c_hat; 
	vector<double> par_optOBdist;
	vector<double> par_varOBdist;

public:
	PartitionScheme(int, DataHandler&);
	~PartitionScheme();

	inline int get_num_targets() { return num_targets; };
	inline double** get_min_risk_mat() { return min_risk_mat; };
	inline double** get_G() { return G; };
	inline int get_num_dstzn() { return num_dstzn; };

	vector<vector<Vertex>> get_nodes_crds();

	void get_risk_innerTrajc();
	double get_risk_innerTrajc(Vertex, Vertex, int);
	double get_reward_innerTrajc(Vertex, Vertex, int);
	void get_risk_outerTrajc();
	double get_risk_outerTrajc(Vertex, Vertex);
	tuple<bool,double,double> is_intersected(Vertex v, Vertex u, int i);

	// risk&reward-related functions
	double inner_risk_function(double l1, double l2, int i);
	double outer_risk_function(double dist);
	// supporting funcitons
	double dot_product(myVector vec1, myVector vec2);
	double get_lineSeg_len(Vertex, Vertex);
	void solve_shortestpath(vector<vector<double>> &, vector<int> &);
};



#endif // !_PAR
