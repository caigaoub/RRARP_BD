#define _USE_MATH_DEFINES
#define INF numeric_limits<double>::infinity()
#define BigM  100000
#define SmallM 0.0000001
#include "PartitionScheme.h"

void exit_error(int line_num) {
    fprintf(stderr, "error occured at line %d\n", line_num);
    exit(1);
}

void PartitionScheme::build(DataHandler& dataset, int nb_dstzn, int type_trajc) {
	this->_dataset = &dataset;
	this->_nb_dstzn = nb_dstzn;
	this->_type_trajc = type_trajc;

	_size_G = 2 * _dataset->_nb_targets * _nb_dstzn + 2;
	_G.resize(_size_G);
	for (int i = 0; i < _size_G; i++) {
		_G[i].resize(_size_G, make_pair(false, 0));
	}
	_min_risk_tars.resize(_dataset->_nb_targets + 2);
	for (int i = 0; i < _dataset->_nb_targets + 2; i++) {
		_min_risk_tars[i].resize(_dataset->_nb_targets + 2, 0); // either 0 or positive for building TSP model
	}
	this->_subarc_angle = 2.0 * M_PI / _nb_dstzn;
	// parameters initialization for risk & reward functions
	for (int i = 0; i < _dataset->_nb_targets; i++){
		_par_c.push_back(_dataset->_radii[i]*_dataset->_radii[i]);
	}
	
	for (int i = 0; i < _dataset->_nb_targets; i++) //let the best observation distance be half radius
		_par_optOBdist.push_back(_dataset->_radii[i]/2.0);
	for (int i = 0; i < _dataset->_nb_targets; i++)
		_par_varOBdist.push_back(1.0);

	// double MINDIST_C2C = get_lineSeg_len(_dataset->_depot1_loc,_dataset->_depot2_loc);
	// double tempD = 0;
	// for(int i = 0; i < _dataset->_nb_targets; i++){
	// 	tempD = get_lineSeg_len(_dataset->_depot1_loc,_dataset->_target_locs[i]);
	// 	if(MINDIST_C2C > tempD){
	// 		MINDIST_C2C = tempD;
	// 	}
	// 	tempD = get_lineSeg_len(_dataset->_depot2_loc,_dataset->_target_locs[i]);
	// 	if(MINDIST_C2C > tempD){
	// 		MINDIST_C2C = tempD;
	// 	}
	// }
	// for(int i = 0; i < _dataset->_nb_targets; i++){
	// 	for(int j = i+1; j < _dataset->_nb_targets; j++){
	// 		tempD = get_lineSeg_len(_dataset->_target_locs[i],_dataset->_target_locs[j]);
	// 		if(MINDIST_C2C > tempD){
	// 			MINDIST_C2C = tempD;
	// 		}
	// 	}
	// }
	// _par_h = 1.0/((double)pow(_dataset->_nb_targets,2) * sqrt(2));
	// _par_h = 1.0/MINDIST_C2C;
	// cout << MINDIST_C2C << endl;
	// exit(0);

	_MAX_REWARD_LIN.resize(_dataset->_nb_targets, 0.0);
	_MAX_REWARD_ROT.resize(_dataset->_nb_targets, 0.0);

	// valid edges on graph G is generated by the following steps
	calc_MAX_REWARD();
	calc_avg_innerrisk();
	_par_h = 0.5* _AVG_RISK;
	build_nodes_crds();
	if(_type_trajc == 1){ // linear case
		get_risk_reward_linearInnerTrajc();
	}
	if(_type_trajc == 2){ // rotatory case
		get_risk_reward_rotatoryInnerTrajc();
	}
	if(_type_trajc != 1 && _type_trajc != 2){
		cerr << "ERROR: No such trajectory type(PartitionScheme.cpp::line 79)" << endl;
		exit(0);
	}
	get_risk_outerTrajc(); // **Important**:'_min_risk_tars' is generated in this function.
}

/**
 * discretize boundaries into sets of points
 * Out: _points(type: vector<vector<Vertex>>, rows: _nb_targets+2, cols: points on it)	
*/
void PartitionScheme::build_nodes_crds() {
	_points.resize(_dataset->_nb_targets + 2);
	_points[0].resize(1); //departure point
	_points[0][0] = _dataset->_depot1_loc;
	// _points[0][0].print();
	for (int pos = 0; pos < _dataset->_nb_targets; pos++) {
		_points[pos+1].resize(_nb_dstzn);
		for (int i = 0; i < _nb_dstzn; i++) {
			_points[pos+1][i]._x = _dataset->_radii[pos] * cos(i * _subarc_angle) + _dataset->_target_locs[pos]._x;
			_points[pos+1][i]._y = _dataset->_radii[pos] * sin(i * _subarc_angle) + _dataset->_target_locs[pos]._y;
			// _points[pos+1][i].print();
		}
	}
	_points[_dataset->_nb_targets + 1].resize(1); // landing depot
	_points[_dataset->_nb_targets + 1][0] = _dataset->_depot2_loc;
	// _points[_dataset->_nb_targets + 1][0].print();
}

/**
 * function: calculate the max reward on linear or rotatory trajectories that can be collected from each target
 * out: _MAX_REWARD_LIN (vector<double>, max reward on linear trajectory for target i)
 * out:	_MAX_REWARD_ROT (vector<double>, max reward on rotatory trajectory)
*/
void PartitionScheme::calc_MAX_REWARD(){
	if(_test_mod){
		for(int i = 0; i< _dataset->_nb_targets; i++){
			_MAX_REWARD_LIN[i] = 2.0*_dataset->_radii[i];
		}
		for(int i = 0; i< _dataset->_nb_targets; i++){
			_MAX_REWARD_ROT[i] = 2.0 * M_PI *_dataset->_radii[i];
		}
	}else{
		/* calculate max reward on linear trajectories*/
		double nb_chords = 100.0;
		double len_cur_chord = 0.0;
		double nb_SLSs = 200.0;//nb of sub line segments
		double len_SLS = 0.0;
		double total_reward = 0.0, point_reward = 0.0;
		for(int i = 0; i< _dataset->_nb_targets; i++){
			double len_chord_interval = 2.0*_dataset->_radii[i]/nb_chords;
			for(int k = 1; k < nb_chords; k++){
				len_cur_chord = len_chord_interval * k; 
				total_reward = 0.0;
				len_SLS = 0.5*len_cur_chord/nb_SLSs;
				double y_yoi_square = pow(_dataset->_radii[i],2) - pow(len_cur_chord/2.0, 2);
				for(int j = 0; j < nb_SLSs; j++){
					point_reward = len_SLS*exp(-pow(sqrt(pow(j*len_SLS,2) + y_yoi_square)- _par_optOBdist[i],2)/(2.0*_par_varOBdist[i]));
					// cout << point_reward << ", ";
					total_reward += point_reward;
				}
				// cout << endl;
				// cout << 2.0 * total_reward << '\t';
				_MAX_REWARD_LIN[i] = max(_MAX_REWARD_LIN[i], 2.0 * total_reward);
			}
			// cout << "target " << i << " (LINEAR): " << _MAX_REWARD_LIN[i] << endl;
		}
		double opt_z = 0.0;
		/* calculate max reward on rotation trajectories */
		for(int i = 0; i< _dataset->_nb_targets; i++){
			opt_z = 0.5 * (_par_optOBdist[i] + sqrt(_par_optOBdist[i]*_par_optOBdist[i] + 4.0 * _par_varOBdist[i] * _par_varOBdist[i]));
			opt_z = min(opt_z, _dataset->_radii[i]);
			_MAX_REWARD_ROT[i] = 2*M_PI*opt_z* exp(-pow(opt_z-_par_optOBdist[i],2)/(2.0*_par_varOBdist[i]*_par_varOBdist[i]));
			// cout << "target " << i << " (ROTATORY): " << _MAX_REWARD_ROT[i] << endl;
		}
	}
}

/**
 * function: calculate the averger inner, _AVG_RISK = \sum(average risk in target i)
 * out: _AVG_RISK (double)
*/
void PartitionScheme::calc_avg_innerrisk(){
	_AVG_RISK = 0.0;
	for(int i = 0; i < _dataset->_nb_targets; i++){
		_AVG_RISK += _par_c[i] -pow(_dataset->_radii[i],2)/3.0;
	}
	_AVG_RISK /= (double)_dataset->_nb_targets;
	// cout << "AVG RISK: " << _AVG_RISK << endl;
	// exit(0);
}

/* 
	calculate risk&reward on all inner circle trajectories.
*/
void PartitionScheme::get_risk_reward_linearInnerTrajc() {
	int idx_row, idx_col, flag;
	// double max_reward;
	// double max_risk;
	vector<double> risk_innerpath(_nb_dstzn, -1); // -1 means invalid risk or reward value
	vector<double> reward_innerpath(_nb_dstzn, -1);
	for (int s = 0; s < _dataset->_nb_targets; s++) { 
		// max_reward = 0.0;
		// max_risk = 0.0;
		for (int i = 0; i < _nb_dstzn; i++) {
			if (i == 0) { // entry point and exit point are the same
				risk_innerpath[i] = 0.0; // 0 means no risk or no reward
				reward_innerpath[i] = 0.0; //
			}
			else {
				if(_test_mod){
					/*  Test: the correctness by simply using euclidean distance */
					risk_innerpath[i] = get_lineSeg_len(_points[s+1][0], _points[s+1][i]);
					reward_innerpath[i] = get_lineSeg_len(_points[s+1][0], _points[s+1][i]);
				}else{
					/*  REAL: real risk evaluation over linear inner trajectory*/
					risk_innerpath[i] = get_risk_linearInnerTrajc(_points[s+1][0], _points[s+1][i], s);
					reward_innerpath[i] = get_reward_linearInnerTrajc(_points[s+1][0], _points[s+1][i], s);
				}
				// cout << "target: " << s+1 << ": " <<risk_innerpath[i] << " " << reward_innerpath[i] << endl; 
				// if (risk_innerpath[i] > max_risk && risk_innerpath[i] < BigM)
				// if (risk_innerpath[i] < BigM)
				// 	max_risk = risk_innerpath[i];
				// if (reward_innerpath[i] > max_reward && reward_innerpath[i] < BigM){
				// 	max_reward = reward_innerpath[i];
				// }
			}
		}
		// cout << "Target: " << s+1 << " max risk: " <<  max_risk << endl;
		// cout << "Target: " << s+1 << " max reward: " << max_reward << endl;
		/* write admissible trajectories (meet minimum reward requirement) to matrix G */ 
		idx_row = s * 2 * _nb_dstzn + 1; 
		idx_col = s * 2 * _nb_dstzn + _nb_dstzn + 1;
		for (int i = 0; i < _nb_dstzn; i++) {
			for (int j = 0; j < _nb_dstzn; j++) {
				flag = (_nb_dstzn - i + j) % _nb_dstzn;
				if (reward_innerpath[flag] >= _MAX_REWARD_LIN[s]*_dataset->_bdg_rewards_ratio[s]){
					// careful: weight on inner risk is 3 times as outer risk
					_G[idx_row + i][idx_col + j] = make_pair(true, risk_innerpath[flag]);
					_nb_adm_InT++;
				}else{
					;
					// _G[idx_row + i][idx_col + j] = make_pair(false, INF);
				}
			}
		}
	}
	/*print for debugging: */
	if(false){
		for (int i = 0; i < _size_G; i++) {
			for (int j = 0; j < _size_G; j++) {
				cout << _G[i][j].second << ' ';
			}
			cout << '\n';
		}
	}
	// exit(0);
}



double PartitionScheme::get_reward_linearInnerTrajc(Vertex entry, Vertex exit, int tar) { 
	double l1 = get_lineSeg_len(entry, exit) / 2.0;
	// cout << l1 << "   ";
	double l2 = sqrt(_dataset->_radii[tar] * _dataset->_radii[tar] - l1 *l1 + 0.0000000001);
	double num_intervals = 300.0; // partition the half chord into num_interval pieces
	double len_interval = l1 / num_intervals;
	double total_reward = 0.0;
	double d, point_reward;
	for (int j = 0; j < num_intervals; j++) {
		d = sqrt(l2 *l2 + pow(len_interval * j, 2));
		point_reward = len_interval* exp(-pow(d - _par_optOBdist[tar], 2) / (2.0 * pow(_par_varOBdist[tar], 2)));
		total_reward += point_reward;
	}
	total_reward *= 2.0;
	return total_reward;
}

double PartitionScheme::get_risk_linearInnerTrajc(Vertex entry, Vertex exit, int tar) {
	double l1 = get_lineSeg_len(entry, exit) / 2.0;
	double l2 = sqrt(_dataset->_radii[tar] * _dataset->_radii[tar] - l1* l1 + 0.00000000001);
	return inner_risk_function(l1, l2, tar);
}

/**
 * get the risk & reward on inner rotatory trajectories
*/
void PartitionScheme::get_risk_reward_rotatoryInnerTrajc() {
	int idx_row, idx_col, flag;
	vector<tuple<bool,double,double>> admissible_path(_nb_dstzn);
	// tuple<double,double> tmpRR;
	for (int s = 0; s < _dataset->_nb_targets; s++) { 
		for (int i = 0; i < _nb_dstzn; i++) {	
			if(_test_mod){
				admissible_path[i] = get_optEucl_rotatoryInnerTrajc(i*_subarc_angle, s);
			}else{
				/* search both ccw and cw  */
				admissible_path[i] = get_optimal_rotatoryInnerTrajc(i*_subarc_angle, s);
			}
		}
		/* write admissible trajectories (meet minimum reward requirement) to matrix G */ 
		idx_row = s * 2 * _nb_dstzn + 1; 
		idx_col = s * 2 * _nb_dstzn + _nb_dstzn + 1;
		for (int i = 0; i < _nb_dstzn; i++) {
			for (int j = 0; j < _nb_dstzn; j++) {
				flag = (_nb_dstzn - i + j) % _nb_dstzn;
				if (get<0>(admissible_path[i])){
					_G[idx_row + i][idx_col + j] = make_pair(true, get<1>(admissible_path[flag]));
					_nb_adm_InT++;
				}else{
					_G[idx_row + i][idx_col + j] = make_pair(false, INF);

				}
			}
		}
	}
	// exit(0);

	/*print for debugging: */
	if(false){
		for (int i = 0; i < _size_G; i++) {
			for (int j = 0; j < _size_G; j++) {
				cout << _G[i][j].second << ' ';
			}
			cout << '\n';
		}
	}
	// exit(0);
}
tuple<bool, double,double> PartitionScheme::get_optimal_rotatoryInnerTrajc(double phi, int tar){
	double nbp_onRadius = 100.0;
	double unit_rad = _dataset->_radii[tar]/nbp_onRadius;
	double tmpRad = 0.0, tmpReward = 0.0, tmpRisk = 0.0;
	double optRad = 0.0, optReward = 0.0, optRisk = INF;
	bool admissible = false;
	for(int i = 1; i <= nbp_onRadius; i++){
		tmpRad = i * unit_rad;
		/*rotate clockwise*/
		tmpReward = phi * tmpRad * (exp(-pow(tmpRad-_par_optOBdist[tar],2)/(2.0*pow(_par_varOBdist[tar],2))));
		if(tmpReward > _MAX_REWARD_ROT[tar]*_dataset->_bdg_rewards_ratio[tar]){
			tmpRisk = 0.0;
			tmpRisk += 2.0*(_par_c[tar]*(_dataset->_radii[tar]-tmpRad) - 1.0/3.0 *(pow(_dataset->_radii[tar],3)-pow(tmpRad,3)));
			tmpRisk += phi * tmpRad * (_par_c[tar] - pow(tmpRad,2));
			if(tmpRisk < optRisk){
				admissible = true;
				optRisk = tmpRisk;
				optRad = tmpRad;
				optReward = tmpReward;
			}
		}
		/*rotate counterclockwise*/
		tmpReward = (2.0*M_PI- phi) * tmpRad * (exp(-pow(tmpRad-_par_optOBdist[tar],2)/(2.0*pow(_par_varOBdist[tar],2))));
		if(tmpReward > _MAX_REWARD_ROT[tar]*_dataset->_bdg_rewards_ratio[tar]){
			tmpRisk = 0.0;
			tmpRisk += 2.0*(_par_c[tar]*(_dataset->_radii[tar]-tmpRad) - 1.0/3.0 *(pow(_dataset->_radii[tar],3)-pow(tmpRad,3)));
			tmpRisk += (2.0*M_PI- phi) * tmpRad * (_par_c[tar] - pow(tmpRad,2));
			if(tmpRisk < optRisk){
				admissible = true;
				optRisk = tmpRisk;
				optRad = tmpRad;
				optReward = tmpReward;
			}
		}
	}
	// cout << "Radius: " << _dataset->_radii[tar] << " opt depth: " << optRad << endl; 
	return make_tuple(admissible, optRisk, optReward);
}


tuple<bool, double,double> PartitionScheme::get_optEucl_rotatoryInnerTrajc(double phi, int tar){
	double nbp_onRadius = 100.0;
	double unit_rad = _dataset->_radii[tar]/nbp_onRadius;
	double tmpRad = 0.0, tmpReward = 0.0, tmpRisk = 0.0;
	double optReward = 0.0, optRisk = INF;
	double optRad = 0.0;
	bool admissible = false;
	for(int i = 1; i <= nbp_onRadius; i++){
		tmpRad = i * unit_rad;
		/*rotate clockwise*/
		tmpReward = phi * tmpRad;
		if(tmpReward > _MAX_REWARD_ROT[tar]*_dataset->_bdg_rewards_ratio[tar]){
			tmpRisk = 2.0*(_dataset->_radii[tar]-tmpRad) + phi * tmpRad;
			// cout << tmpRisk << endl;
			if(tmpRisk < optRisk){
				admissible = true;
				optRisk = tmpRisk;
				optRad = tmpRad;
				optReward = tmpReward;
			}
		}
		// cout << tmpRad << ", " << optRisk << ", " << optReward << endl;
		/*rotate counterclockwise*/
		tmpReward = (2.0*M_PI- phi) * tmpRad;
		if(tmpReward > _MAX_REWARD_ROT[tar]*_dataset->_bdg_rewards_ratio[tar]){
			tmpRisk = 2.0*(_dataset->_radii[tar]-tmpRad) + (2.0*M_PI- phi) * tmpRad;
			// cout << tmpRisk << endl;
			if(tmpRisk < optRisk){
				admissible = true;
				optRisk = tmpRisk;
				optRad = tmpRad;
				optReward = tmpReward;
			}
		}
		// cout << tmpRad << ", " << optRisk << ", " << optReward << endl;
	}
	// exit(0);
	// cout << "Radius: " << _dataset->_radii[tar] << " opt depth: " << optRad << " opt risk: " << optRisk << " opt reward: " << optReward << endl; 
	return make_tuple(admissible, optRisk, optReward);
}

void PartitionScheme::get_risk_outerTrajc() {
	/* Generate the accumulated risk on outer paths*/
	int idxmat_1, idxmat_2, idxmat_3, idxmat_4;
	double val_risk, val_min_risk;
	/* i) departure depot --> each entry turning point of a target */	
	for (int t = 1; t <= _dataset->_nb_targets; t++) {	
		// cout << 0 << " ---> " << t << endl;
		idxmat_1 = (t - 1) * 2 * _nb_dstzn + 1; 
		val_min_risk = INF;
		for (int i = 0; i < _nb_dstzn; i++) {
			if(_test_mod){
				val_risk = get_lineSeg_len(_points[t][i], _dataset->_depot1_loc);
			}else{
				val_risk = get_risk_outerTrajc(_points[t][i], _dataset->_depot1_loc);
			}
			_G[0][idxmat_1 + i] = make_pair(true, val_risk);
			_nb_adm_OutT += 1;

			// cout << _G[0][idxmat_1 + i].second << " ";
			if (val_risk < val_min_risk) {
				val_min_risk = val_risk;
			}
		}
		/*minimum risk matrix is symmetic*/
		_min_risk_tars[0][t] = val_min_risk;
		// _min_risk_tars[t][0] = val_min_risk;
		_min_risk_tars[t][0] = 0;

		// cout << " ====>> minimum risk from departure to target " << t << ": " << val_min_risk << endl;
	}
	// exit(0);
	/* each pair of entry and exit between targets*/
	for (int s = 1; s <= _dataset->_nb_targets; s++) { // ii) boundary s <-> boudary t
		idxmat_1 = (s - 1) * 2 * _nb_dstzn + _nb_dstzn + 1; // exit of target s
		idxmat_3 = (s - 1) * 2 * _nb_dstzn + 1; // entries  of target s
		for (int t = 1; t <= _dataset->_nb_targets; t++) {
			val_min_risk = INF;
			if(s != t){
				// cout << s << " ---> " << t << endl;

				idxmat_2 = (t - 1) * 2 * _nb_dstzn + 1; // vertices on boundary t acted as entries
				idxmat_4 = (t - 1) * 2 * _nb_dstzn + _nb_dstzn + 1; // vertices on boundary t acted as exits
				for (int i = 0; i < _nb_dstzn; i++) {
					for (int j = 0; j < _nb_dstzn; j++) {
						if(_test_mod){
							val_risk = get_lineSeg_len(_points[s][i], _points[t][j]); // for testing
						}else{
							val_risk = get_risk_outerTrajc(_points[s][i], _points[t][j]);
						}
						_G[idxmat_1 + i][idxmat_2 + j] = make_pair(true, val_risk);
						_nb_adm_OutT += 1;
						_G[idxmat_4 + j][idxmat_3 + i] = make_pair(true, val_risk);
						_nb_adm_OutT += 1;
						// cout<< s << " " << t << ": " << idxmat_1+i << "-->" << idxmat_2+j << " dist: "<< _G[idxmat_1 + i][idxmat_2 + j].second << endl;
						// cout << idxmat_1+i << "-->" << idxmat_2+j << " : " << _G[idxmat_1 + i][idxmat_2 + j].second << '\n';
						if (val_risk < val_min_risk)
							val_min_risk = val_risk;
					}
				}
				_min_risk_tars[s][t] = val_min_risk;
				_min_risk_tars[t][s] = val_min_risk;
				// cout << " ------ min-risk target s to target t: " << val_min_risk << endl;
			}
		}
	}
	/*all exits to arrival depot*/
	for (int s = 1; s <= _dataset->_nb_targets; s++) { 	
		// cout << s << " ---> " << _dataset->_nb_targets+1 << endl;
		idxmat_2 = (s - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
		val_min_risk = INF;
		for (int i = 0; i < _nb_dstzn; i++) {
			if(_test_mod){
				val_risk = get_lineSeg_len(_points[s][i], _dataset->_depot2_loc); // for testing 
			}else{
				val_risk = get_risk_outerTrajc(_points[s][i], _dataset->_depot2_loc);
			}
			_G[idxmat_2 + i][_size_G - 1] = make_pair(true, val_risk);
			_nb_adm_OutT += 1;

			if (val_risk < val_min_risk)
				val_min_risk = val_risk;
		}
		_min_risk_tars[s][_dataset->_nb_targets + 1] = val_min_risk;
		// _min_risk_tars[_dataset->_nb_targets + 1][s] = val_min_risk;
		_min_risk_tars[_dataset->_nb_targets + 1][s] = 0;
		// cout << " ------ depot 2 to target: " << val_min_risk << endl;
	}


	/* diagonal elements -- 
	for (s = 0; s <= _dataset->_nb_targets + 1; s++) {
		_min_risk_tars[s][s] = 0.0;
	}
	_min_risk_tars[0][_dataset->_nb_targets + 1] = 0.0;
	_min_risk_tars[_dataset->_nb_targets + 1][0] = 0.0;
	*/

	/*
	ofstream file_G("./matrixG.txt");
	for (int i = 0; i < _size_G; i++) {
		for (int j = 0; j < _size_G; j++) {
			file_G << _G[i][j].second << '\t';
		}
		file_G << '\n';
	}
	*/

	// ------------  Step 2: subtract the _min_risk_tars[s][t] from _G[i][j]   -------------
	if(true){
		for (int t = 1; t <= _dataset->_nb_targets; t++) {
		idxmat_1 = (t - 1) * 2 * _nb_dstzn + 1;
		// idxmat_2 = (t - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
		for (int i = 0; i < _nb_dstzn; i++) {
			if(_G[0][idxmat_1 + i].first == true){
				// cout << "D" << " " << t << ": " << "D" << "-->" << idxmat_1+i << " dist: ";
				// cout << _G[0][idxmat_1 + i].second << " "<< _min_risk_tars[0][t] << " ";
				_G[0][idxmat_1 + i].second -= _min_risk_tars[0][t];
				// cout << " = " << _G[0][idxmat_1 + i].second << endl;
			}

		}
		}
		for (int s = 1; s <= _dataset->_nb_targets; s++) {
			idxmat_1 = (s - 1) * 2 * _nb_dstzn + _nb_dstzn + 1; // vertices worked as exits in boundary s
			idxmat_3 = (s - 1) * 2 * _nb_dstzn + 1; // vertices worked as entries in boundary s
			for (int t = 1; t <= _dataset->_nb_targets; t++) {
				if(s != t && s < t){
					idxmat_2 = (t - 1) * 2 * _nb_dstzn + 1; // vertices worked as entries on boundary t
					idxmat_4 = (t - 1) * 2 * _nb_dstzn + _nb_dstzn + 1; // vertices worked as exits on boundary t
					for (int i = 0; i < _nb_dstzn; i++) {
						for (int j = 0; j < _nb_dstzn; j++) {
							if(_G[idxmat_1 + i][idxmat_2 + j].first == true){
								// cout<< s << " " << t << ": " << idxmat_1+i << "-->" << idxmat_2+j << " dist: ";
								// cout << _G[idxmat_1 + i][idxmat_2 + j].second << " "<< _min_risk_tars[s][t] << " ";
								_G[idxmat_1 + i][idxmat_2 + j].second -= _min_risk_tars[s][t];
								// cout << " = " << _G[idxmat_1 + i][idxmat_2 + j].second << endl;
							}
							if(_G[idxmat_4 + j][idxmat_3 + i].first == true){
								// cout<< t << " " << s << ": " << idxmat_4+j << "-->" << idxmat_3+i << " dist: ";
								// cout << _G[idxmat_4 + j][idxmat_3 + i].second << " "<< _min_risk_tars[t][s] << " ";
								_G[idxmat_4 + j][idxmat_3 + i].second -= _min_risk_tars[t][s];
						
								// cout << " = " << _G[idxmat_4 + j][idxmat_3 + i].second << endl;
							}
						}
					}	
				}			
			}
		}
		for (int s = 1; s <= _dataset->_nb_targets; s++) {
			idxmat_2 = (s - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
			for (int i = 0; i < _nb_dstzn; i++) {
				if(_G[idxmat_2 + i][_size_G - 1].first == true){
					// cout<< s << " " << "A" << ": " << idxmat_2+i << "-->" << _size_G - 1 << " dist: ";
					// cout << _G[idxmat_2 + i][_size_G - 1].second << " "<< _min_risk_tars[s][_dataset->_nb_targets + 1] << " ";
					_G[idxmat_2 + i][_size_G - 1].second -= _min_risk_tars[s][_dataset->_nb_targets + 1];
					// cout << " = " << _G[idxmat_2 + i][_size_G - 1].second << endl;
				}
			}
		}
	}
	

	/*print for debugging: */
	if(false){
		for (int i = 0; i < _size_G; i++) {
			for (int j = 0; j < _size_G; j++) {
				cout << _G[i][j].second << ' ';
			}
			cout << '\n';
		}
	}
	if(false){
		for (int s = 0; s < _dataset->_nb_targets+2; s++) {
			for (int t = 0; t < _dataset->_nb_targets+2; t++) {
				cout << _min_risk_tars[s][t] << ' ';
			}
			cout << '\n';
		}
	}
	// cout << "nb of admissible trajectories: " << _nb_adm_InT << ", " << _nb_adm_OutT << endl;
	// exit(0);

}

double PartitionScheme::get_risk_outerTrajc(Vertex v, Vertex u) {
	// double total_risk = get_lineSeg_len(v, u);
	// v.print();
	// u.print();
	double total_risk = outer_risk_function(get_lineSeg_len(v, u)); // total risk over the line segment
	// cout << total_risk << "_*_";
	for (int i = 0; i < _dataset->_nb_targets; i++) { // additional possible risk caused when passing through some region(s)
		tuple<bool, double, double> result = is_intersected(v, u, i);
		// cout << get<0>(result) << "**" << get<1>(result) << "**"<<get<2>(result) << endl;
		if (get<0>(result) == true){ // if current straight line intersects with region Ui, we add additional risk to this outer path
			// cout << "intersected target " << i+1 << endl;
			// total_risk = total_risk + get<1>(result) - get<2>(result);
			total_risk += _par_h/2.0 * get<1>(result); // add penalty
		}
	}
	// cout << total_risk << " ";
	// exit(0);
	return total_risk;
}

tuple<bool, double,double> PartitionScheme::is_intersected(Vertex v, Vertex u, int tar) {
	// current target center is denoted as o;
	Vertex o = { _dataset->_target_locs[tar]._x, _dataset->_target_locs[tar]._y };
	myVector v_2_o(v, o); // vector v->o
	myVector v_2_u(v, u); // vector v->u
	myVector o_2_u(o, u); // vector o->u
	double len_v_2_o = v_2_o.get_vecLen();
	double len_v_2_u = v_2_u.get_vecLen();
	double len_o_2_u = o_2_u.get_vecLen();	
	double angle_u = acos((len_o_2_u*len_o_2_u + len_v_2_u*len_v_2_u-len_v_2_o*len_v_2_o)/(2.0*len_o_2_u*len_v_2_u));
	double l2 = len_o_2_u * sin(angle_u); // vertical distance from center o to line segment vu
	double l1; // half length of the chord
	// double risk_vu_on_Ui; // risk of path vu when passing through target region Ui
	if (l2 + 0.0000001 < _dataset->_radii[tar] && angle_u < (M_PI / 2.0) && len_o_2_u * cos(angle_u)+ 0.0000001 < len_v_2_u) {
		l1 = sqrt(_dataset->_radii[tar] * _dataset->_radii[tar] - l2 * l2);
		// risk_vu_on_Ui = inner_risk_function(l1, l2, tar);
		// cout << "vo:" << len_v_2_o << " vu:" << len_v_2_u << " ou:" << len_o_2_u << endl;
		// cout << "v:" << v._x << "," << v._y << endl;
		// cout << "u:" << u._x << "," << u._y << endl;
		// cout << "o:" << o._x << "," << o._y << endl;
		// cout << len_v_2_u << " " << l1 << " " << l2 <<  " r=" << _dataset->_radii[tar] << endl;
		// return make_tuple(true, risk_vu_on_Ui, outer_risk_function(2 * l1));

		return make_tuple(true, 2.0*l1, outer_risk_function(2 * l1));

	}
	else {
		return make_tuple(false, NULL, NULL);
	}
}

double PartitionScheme::inner_risk_function(double l1, double l2, int tar) {
	double T = -2.0/3.0*pow(l1,3) + 2.0*(_par_c[tar] - pow(l2,2))*l1; // for function f_i = -(x-xoi)^2 - (y - yoi)^2 + c
	// T /= _par_c[tar];	 // normalize to [0, 1]
	// T *= 100.0;
	// cout << T << endl;

	return T;
}

double PartitionScheme::outer_risk_function(double dist) {
	double R = _par_h * dist;
	// cout << R << endl;
	return R;
	
}

double PartitionScheme::dot_product(myVector vec1, myVector vec2) {
	return vec1._x * vec2._x + vec1._y * vec2._y;
}


double PartitionScheme::get_lineSeg_len(Vertex v, Vertex u) {
	return sqrt(pow(u._x - v._x, 2) + pow(u._y - v._y, 2));
}

void PartitionScheme::solve_shortestpath(vector<int> & seq, vector<vector<double>> & SDS) {
	int idxmat_1, idxmat_2, idxmat_3;
	vector<double> entry_dist(_nb_dstzn, INF);
	vector<double> exit_dist(_nb_dstzn, INF);
	vector<double> INF_dist(_nb_dstzn, INF);
	double dist_endDepot = INF;
	double dist = 0.0;
	// i) SDS[start depot] = 0.0
	SDS[0].resize(1);
	SDS[0][0] = 0.0;
	// ii) SDS[*] of boundary points of the first visited target
	int pos = 1;
	int idx_circle = seq[pos];
	SDS[1].resize(_nb_dstzn);
	idxmat_1 = (idx_circle - 1) * 2 * _nb_dstzn + 1; //first entry of the first visited target
	for (int i = 0; i < _nb_dstzn; i++)	{
		if (_G[0][idxmat_1 + i].first == true) {
			entry_dist[i] = _G[0][idxmat_1 + i].second;
			SDS[1][i] = entry_dist[i];
		}
	}
	idxmat_2 = (idx_circle - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
	for (int i = 0; i < _nb_dstzn; i++) {
		for (int j = 0; j < _nb_dstzn; j++) {
			if (_G[idxmat_1 + i][idxmat_2 + j].first == true) {
				dist = entry_dist[i] + _G[idxmat_1 + i][idxmat_2 + j].second;
				if (dist < exit_dist[j])
					exit_dist[j] = dist;
			}
		}
	}
	entry_dist = INF_dist;

	// iii) nodes from seq[2] to seq[num_circles]
	int idx_fr_circle = 0;
	int idx_lat_circle = 0;
	for (pos = 2; pos <= _dataset->_nb_targets; pos++) {
		SDS[pos].resize(_nb_dstzn);
		idx_fr_circle = seq[pos - 1]; 
		idx_lat_circle = seq[pos]; 
		idxmat_1 = (idx_fr_circle - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
		idxmat_2 = (idx_lat_circle - 1) * 2 * _nb_dstzn + 1;
		for (int i = 0; i < _nb_dstzn; i++) {
			for (int j = 0; j < _nb_dstzn; j++) {
				if (_G[idxmat_1 + i][idxmat_2 + j].first == true) {
					dist = exit_dist[i] + _G[idxmat_1 + i][idxmat_2 + j].second;
					if (dist < entry_dist[j])
						entry_dist[j] = dist;
				}
			}
		}
		for (int i = 0; i < _nb_dstzn; i++) {
			SDS[pos][i] = entry_dist[i];
		}
		exit_dist = INF_dist;
		idxmat_3 = (idx_lat_circle - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
		// update exits
		for (int i = 0; i < _nb_dstzn; i++) {
			for (int j = 0; j < _nb_dstzn; j++) {
				if (_G[idxmat_2 + i][idxmat_3 + j].first == true) {
					dist = entry_dist[i] + _G[idxmat_2 + i][idxmat_3 + j].second;
					if (dist < exit_dist[j])
						exit_dist[j] = dist;
				}
			}
		}
		entry_dist = INF_dist;
	}
	// iv) the last target <-> the end depot
	idx_circle = seq[_dataset->_nb_targets];
	SDS[_dataset->_nb_targets + 1].resize(1);
	idxmat_1 = (idx_circle - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
	for (int i = 0; i < _nb_dstzn; i++) {
		if (_G[idxmat_1 + i][_size_G - 1].first == true) {
			dist = exit_dist[i] + _G[idxmat_1 + i][_size_G - 1].second;
			if (dist < dist_endDepot)
				dist_endDepot = dist;
		}
	}
	SDS[_dataset->_nb_targets + 1][0] = dist_endDepot;
}


void PartitionScheme::solve_shortestpath_v2(vector<int> & seq, vector<vector<double>> & SDS) {
	int idxmat_1, idxmat_2, idxmat_3;
	vector<double> entry_dist(_nb_dstzn, INF);
	vector<double> exit_dist(_nb_dstzn, INF);
	vector<double> INF_dist(_nb_dstzn, INF);
	double dist_endDepot = INF;
	double dist = 0.0;
	vector<vector<pair<int,int>>> prevnodes_tars(2 * _dataset->_nb_targets); //pair<target, turning point>
	for(int i = 0;i< 2* _dataset->_nb_targets; i++){
		prevnodes_tars[i].resize(_nb_dstzn);
	}

	pair<int,int> prevnode_arrival;
	int col_layer = 0;
	// i) SDS[start depot] = 0.0
	SDS[0].resize(1);
	SDS[0][0] = 0.0;
	// ii) SDS[*] of boundary points of the first visited target
	int pos = 1;
	int idx_circle = seq[pos];
	SDS[1].resize(_nb_dstzn);
	idxmat_1 = (idx_circle - 1) * 2 * _nb_dstzn + 1; //first entry of the first visited target
	for (int i = 0; i < _nb_dstzn; i++)	{
		if (_G[0][idxmat_1 + i].first == true) {
			entry_dist[i] = _G[0][idxmat_1 + i].second;
			SDS[1][i] = entry_dist[i];
			prevnodes_tars[col_layer][i] = make_pair(col_layer,0);
		}
	}

	col_layer++;
	idxmat_2 = (idx_circle - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
	for (int i = 0; i < _nb_dstzn; i++) {
		for (int j = 0; j < _nb_dstzn; j++) {
			if (_G[idxmat_1 + i][idxmat_2 + j].first == true) {
				dist = entry_dist[i] + _G[idxmat_1 + i][idxmat_2 + j].second;
				if (dist < exit_dist[j]){
					exit_dist[j] = dist;
					prevnodes_tars[col_layer][j] = make_pair(idx_circle,i);
				}
			}
		}
	}
	entry_dist = INF_dist;
	// iii) nodes from seq[2] to seq[num_circles]
	int idx_fr_circle = 0;
	int idx_lat_circle = 0;
	col_layer++;
	for (pos = 2; pos <= _dataset->_nb_targets; pos++) {
		SDS[pos].resize(_nb_dstzn);
		idx_fr_circle = seq[pos - 1]; 
		idx_lat_circle = seq[pos]; 
		idxmat_1 = (idx_fr_circle - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
		idxmat_2 = (idx_lat_circle - 1) * 2 * _nb_dstzn + 1;

		for (int i = 0; i < _nb_dstzn; i++) {
			for (int j = 0; j < _nb_dstzn; j++) {
				if (_G[idxmat_1 + i][idxmat_2 + j].first == true) {
					dist = exit_dist[i] + _G[idxmat_1 + i][idxmat_2 + j].second;
					if (dist < entry_dist[j]){
						entry_dist[j] = dist;
						prevnodes_tars[col_layer][j] = make_pair(idx_fr_circle, i);
					}
				}
			}
		}
		col_layer++;

		for (int i = 0; i < _nb_dstzn; i++) {
			SDS[pos][i] = entry_dist[i];
		}
		exit_dist = INF_dist;
		idxmat_3 = (idx_lat_circle - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
		// update exits
		for (int i = 0; i < _nb_dstzn; i++) {
			for (int j = 0; j < _nb_dstzn; j++) {
				if (_G[idxmat_2 + i][idxmat_3 + j].first == true) {
					dist = entry_dist[i] + _G[idxmat_2 + i][idxmat_3 + j].second;
					if (dist < exit_dist[j]){
						exit_dist[j] = dist;
						prevnodes_tars[col_layer][j] = make_pair(idx_lat_circle, i);
					}
				}
			}
		}
		col_layer++;
		entry_dist = INF_dist;
	}
	// iv) the last target <-> the end depot
	idx_circle = seq[_dataset->_nb_targets];
	SDS[_dataset->_nb_targets + 1].resize(1);
	idxmat_1 = (idx_circle - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
	for (int i = 0; i < _nb_dstzn; i++) {
		if (_G[idxmat_1 + i][_size_G - 1].first == true) {
			dist = exit_dist[i] + _G[idxmat_1 + i][_size_G - 1].second;
			if (dist < dist_endDepot){
				dist_endDepot = dist;
				prevnode_arrival = make_pair(idx_circle, i);
			}
		}
	}
	SDS[_dataset->_nb_targets + 1][0] = dist_endDepot;

	/* retrieve the shortest path */
	deque<pair<int,Vertex>> OptPath;
	OptPath.push_front(make_pair(_dataset->_nb_targets+1, _points[_dataset->_nb_targets+1][0]));
	OptPath.push_front(make_pair(prevnode_arrival.first, _points[prevnode_arrival.first][prevnode_arrival.second]));
	// int pre_tar = prevnode_arrival.first;
	int tpidx_pre_tar = prevnode_arrival.second;
	for(int i = 2*_dataset->_nb_targets-1; i >= 0; i--){
		OptPath.push_front(make_pair(prevnodes_tars[i][tpidx_pre_tar].first, _points[prevnodes_tars[i][tpidx_pre_tar].first][prevnodes_tars[i][tpidx_pre_tar].second]));
		tpidx_pre_tar = prevnodes_tars[i][tpidx_pre_tar].second;
	}

	cout << " =====>> optimal entries and exits selected: " << '\n';
	for(auto itr = OptPath.begin(); itr != OptPath.end(); itr++){
		cout << (*itr).first << ": (" << (*itr).second._x << ',' << (*itr).second._y << ")\n";
	}

	ofstream file_OptPath("/home/cai/Dropbox/Box_Research/Github/RRARP_BD/BendersDecomp/ret/OptimalPath.txt");
	for(auto itr = OptPath.begin(); itr != OptPath.end(); itr++){
		file_OptPath << (*itr).first << '\t' << (*itr).second._x << '\t' << (*itr).second._y << '\n';
	}
}


void PartitionScheme::solve_shortestpath_v3(vector<int> & seq, vector<vector<double>> & SDS) {
	int idxmat_1, idxmat_2, idxmat_3;
	vector<double> entry_dist(_nb_dstzn, INF);
	vector<double> exit_dist(_nb_dstzn, INF);
	vector<double> INF_dist(_nb_dstzn, INF);
	double dist_endDepot = INF;
	double dist = 0.0;
	vector<vector<pair<int,int>>> prevnodes_tars(2 * _dataset->_nb_targets); //pair<target, turning point>
	for(int i = 0; i< 2*_dataset->_nb_targets; i++){
		prevnodes_tars[i].resize(_nb_dstzn);
	}

	pair<int,int> prevnode_arrival;
	int col_layer = 0;
	// i) SDS[start depot] = 0.0
	SDS[0].resize(1);
	SDS[0][0] = 0.0;
	// ii) SDS[*] of boundary points of the first visited target
	int pos = 1;
	int idx_circle = seq[pos];
	SDS[1].resize(_nb_dstzn);
	idxmat_1 = (idx_circle - 1) * 2 * _nb_dstzn + 1; //first entry of the first visited target
	for (int i = 0; i < _nb_dstzn; i++)	{
		if (_G[0][idxmat_1 + i].first == true) {
			entry_dist[i] = _G[0][idxmat_1 + i].second;
			SDS[1][i] = entry_dist[i];
			prevnodes_tars[col_layer][i] = make_pair(col_layer,0);
		}
	}

	col_layer++;
	idxmat_2 = (idx_circle - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
	for (int i = 0; i < _nb_dstzn; i++) {
		for (int j = 0; j < _nb_dstzn; j++) {
			if (_G[idxmat_1 + i][idxmat_2 + j].first == true) {
				dist = entry_dist[i] + _G[idxmat_1 + i][idxmat_2 + j].second;
				if (dist < exit_dist[j]){
					exit_dist[j] = dist;
					prevnodes_tars[col_layer][j] = make_pair(idx_circle,i);
				}
			}
		}
	}
	entry_dist = INF_dist;
	// iii) nodes from seq[2] to seq[num_circles]
	int idx_fr_circle = 0;
	int idx_lat_circle = 0;
	col_layer++;
	for (pos = 2; pos <= _dataset->_nb_targets; pos++) {
		SDS[pos].resize(_nb_dstzn);
		idx_fr_circle = seq[pos - 1]; 
		idx_lat_circle = seq[pos]; 
		idxmat_1 = (idx_fr_circle - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
		idxmat_2 = (idx_lat_circle - 1) * 2 * _nb_dstzn + 1;

		for (int i = 0; i < _nb_dstzn; i++) {
			for (int j = 0; j < _nb_dstzn; j++) {
				if (_G[idxmat_1 + i][idxmat_2 + j].first == true) {
					dist = exit_dist[i] + _G[idxmat_1 + i][idxmat_2 + j].second;
					if (dist < entry_dist[j]){
						entry_dist[j] = dist;
						prevnodes_tars[col_layer][j] = make_pair(idx_fr_circle, i);
					}
				}
			}
		}
		col_layer++;

		for (int i = 0; i < _nb_dstzn; i++) {
			SDS[pos][i] = entry_dist[i];
		}
		exit_dist = INF_dist;
		idxmat_3 = (idx_lat_circle - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
		// update exits
		for (int i = 0; i < _nb_dstzn; i++) {
			for (int j = 0; j < _nb_dstzn; j++) {
				if (_G[idxmat_2 + i][idxmat_3 + j].first == true) {
					dist = entry_dist[i] + _G[idxmat_2 + i][idxmat_3 + j].second;
					if (dist < exit_dist[j]){
						exit_dist[j] = dist;
						prevnodes_tars[col_layer][j] = make_pair(idx_lat_circle, i);
					}
				}
			}
		}
		col_layer++;
		entry_dist = INF_dist;
	}
	// iv) the last target <-> the end depot
	idx_circle = seq[_dataset->_nb_targets];
	SDS[_dataset->_nb_targets + 1].resize(1);
	idxmat_1 = (idx_circle - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
	for (int i = 0; i < _nb_dstzn; i++) {
		if (_G[idxmat_1 + i][_size_G - 1].first == true) {
			dist = exit_dist[i] + _G[idxmat_1 + i][_size_G - 1].second;
			if (dist < dist_endDepot){
				dist_endDepot = dist;
				prevnode_arrival = make_pair(idx_circle, i);
			}
		}
	}
	SDS[_dataset->_nb_targets + 1][0] = dist_endDepot;

	/* retrieve the shortest path */
	deque<pair<int,Vertex>> OptPath;
	OptPath.push_front(make_pair(_dataset->_nb_targets+1, _points[_dataset->_nb_targets+1][0]));
	OptPath.push_front(make_pair(prevnode_arrival.first, _points[prevnode_arrival.first][prevnode_arrival.second]));
	// int pre_tar = prevnode_arrival.first;
	int tpidx_pre_tar = prevnode_arrival.second;
	for(int i = 2*_dataset->_nb_targets-1; i >= 0; i--){
		OptPath.push_front(make_pair(prevnodes_tars[i][tpidx_pre_tar].first, _points[prevnodes_tars[i][tpidx_pre_tar].first][prevnodes_tars[i][tpidx_pre_tar].second]));
		tpidx_pre_tar = prevnodes_tars[i][tpidx_pre_tar].second;
	}

	// cout << " =====>> optimal entries and exits selected: " << '\n';
	// for(auto itr = OptPath.begin(); itr != OptPath.end(); itr++){
	// 	cout << (*itr).first << ": (" << (*itr).second._x << ',' << (*itr).second._y << ")\n";
	// }

	ofstream file_OptPath("/home/cai/Dropbox/Box_Research/Github/RRARP_BD/BendersDecomp/ret/OptimalPath2.txt");
	for(auto itr = OptPath.begin(); itr != OptPath.end(); itr++){
		file_OptPath << (*itr).first << '\t' << (*itr).second._x << '\t' << (*itr).second._y << '\n';
	}
}


double PartitionScheme::calc_withdrawal_risk(vector<int> & fseq_){
	double val = 0.0;
	for(int i = 0; i <= _dataset->_nb_targets; i++){
		val += _min_risk_tars[fseq_[i]][fseq_[i+1]];
	}
	return val;
}


void PartitionScheme::calc_risk_C2C(){
	if(_test_mod){
		_risk_C2C.resize(_dataset->_nb_targets+2);
		for(int i = 0; i < _dataset->_nb_targets+2; i++)
			_risk_C2C[i].resize(_dataset->_nb_targets+2, 0.0);
		/*departure depot and arrival depot*/
		for(int i = 0; i < _dataset->_nb_targets; i++){
			_risk_C2C[0][i+1] = get_lineSeg_len(_dataset->_depot1_loc, _dataset->_target_locs[i]);
			_risk_C2C[i+1][_dataset->_nb_targets+1] = get_lineSeg_len(_dataset->_depot2_loc, _dataset->_target_locs[i]);
		}
		for(int i = 0; i < _dataset->_nb_targets; i++){
			for(int j = i+1 ; j < _dataset->_nb_targets; j++){
				_risk_C2C[i+1][j+1] = get_lineSeg_len(_dataset->_target_locs[i], _dataset->_target_locs[j]);
				_risk_C2C[j+1][i+1] = _risk_C2C[i+1][j+1];
			}
		}
	}else{
		_risk_C2C.resize(_dataset->_nb_targets+2);
		for(int i = 0; i < _dataset->_nb_targets+2; i++)
			_risk_C2C[i].resize(_dataset->_nb_targets+2);
		
		/* risk evaluated on radii*/
		// vector<double> risk_onRadii(_dataset->_nb_targets, 0);
		// for(int i = 0; i < _dataset->_nb_targets; i++)
		// 	risk_onRadii[i] = _par_c[i]*_dataset->_radii[i] - pow(_dataset->_radii[i],3)/3.0;
		
		/*departure depot and arrival depot*/
		for(int i = 0; i < _dataset->_nb_targets; i++){
			_risk_C2C[0][i+1] = _min_risk_tars[0][i+1];
			_risk_C2C[i+1][_dataset->_nb_targets+1] = _min_risk_tars[i+1][_dataset->_nb_targets+1];
		}
		for(int i = 0; i < _dataset->_nb_targets; i++){
			for(int j = i+1 ; j < _dataset->_nb_targets; j++){
				_risk_C2C[i+1][j+1] = _min_risk_tars[i+1][j+1];
				_risk_C2C[j+1][i+1] = _risk_C2C[i+1][j+1];
			}
		}
	}
	
	if(false){
		for (int s = 0; s < _dataset->_nb_targets+2; s++) {
			for (int t = 0; t < _dataset->_nb_targets+2; t++) {
				cout << _risk_C2C[s][t] << ' ';
			}
			cout << '\n';
		}
	}
}
