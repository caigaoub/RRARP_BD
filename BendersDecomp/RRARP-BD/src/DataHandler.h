#ifndef _DATAHANDLER_H_
#define _DATAHANDLER_H_

#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "myNameClass.h"
using namespace std;
struct Vertex {
	double x;
	double y;

	Vertex& operator= (Vertex other) {
		swap(x, other.x);
		swap(y, other.y);
		return *this;
	}
};

class DataHandler {
public:
	int				_nb_targets;
	int				_nb_clusters = 0; // default no cluster
	Vertex			_depot1_loc;
	Vertex			_depot2_loc;
	Vertex*			_target_locs;
	double*			_radii;
	double*			_bdg_rewards_ratio; // budgetary rewards ratio
	double*			_risk_thold_ratio; // maximum risk taken at that inner trajectory	

	DataHandler() {};
	~DataHandler();
	void parse(const char* filename, bool is_cluster);
	void print();
};

#endif // !_DATAHANDLER_H_
