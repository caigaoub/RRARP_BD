#ifndef _DATAHANDLER_H_
#define _DATAHANDLER_H_

#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include "myNameClass.h"
using namespace std;

class DataHandler {
private:
	int num_targets;
	Vertex depot1_loc;
	Vertex depot2_loc;
	double* radii;
	double* bdg_rewards; // budgetary rewards
	double* risk_thold;
	Vertex* target_locs;

public:
	DataHandler(const char* filename);
	~DataHandler();

	inline int get_num_targets() { return num_targets; }
	inline Vertex get_depot1_loc() { return depot1_loc; }
	inline Vertex get_depot2_loc() { return depot2_loc; }
	inline Vertex * get_target_locs() { return target_locs; }
	inline double * get_radii() { return radii; }
	inline double * get_risk_thold() { return	risk_thold; }
	inline double * get_bdg_rewards() { return bdg_rewards; };
	void print_data_info();

};


#endif // !_DATAHANDLER_H_
