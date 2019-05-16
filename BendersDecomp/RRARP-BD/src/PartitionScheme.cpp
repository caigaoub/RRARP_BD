#include "DataHandler.h"
#include "PartitionScheme.h"
#include "myNameClass.h"
#include <cmath>
#include <math.h>

#define PI atan(1.0) * 4
#define INF numeric_limits<double>::infinity()

PartitionScheme::PartitionScheme(int k, DataHandler& instance) {
	/* retriving basic data info from current instance */
	this->num_dstzn = k;
	this->num_targets = instance.get_num_targets();
	this->depot1_loc = instance.get_depot1_loc();
	this->depot2_loc = instance.get_depot2_loc();

	/* for index ease, let all related data start counting from 1 to num_targets */
	target_locs.push_back({ INF, INF }); 
	for (int i = 0; i < num_targets; i++) {
		target_locs.push_back(instance.get_target_locs()[i]);
	}
	radii.push_back(INF);
	for (int i = 0; i < num_targets; i++)
		radii.push_back(instance.get_radii()[i]);
	bdg_rewards.push_back(INF); 
	for (int i = 0; i < num_targets; i++)
		bdg_rewards.push_back(instance.get_bdg_rewards()[i]);
	risk_thold.push_back(INF);
	for (int i = 0; i < num_targets; i++)
		risk_thold.push_back(instance.get_risk_thold()[i]);

	/* Graph: G = (V, E), where 
      - node set V: consisting of start depot, end depot, all entries and all exits
	  - edge set (p,q) in E: a) 0 means no connection existed or zero risk if connection is allowed;
	                         b) INF means inadmissible path from p to q;
							 c) positive value means path is admissible;
	*/
	num_V = 2 * num_targets * k + 2;
	G = new double*[num_V];
	for (int i = 0; i < num_V; i++) {
		G[i] = new double[num_V]();
	}
	/* every entity with positive value represents the minimum risk between boundary i and boundary j */
	min_risk_mat = new double*[num_targets + 2];
	for (int i = 0; i < num_targets + 2; i++) {
		min_risk_mat[i] = new double[num_targets + 2];
	}
	subarc_angle = 2.0 * PI / num_dstzn;

	// parameters initialization for risk&reward functions
	par_c1.push_back(0); 
	for (int i = 0; i < num_targets; i++)
		par_c1.push_back(1);
	par_optOBdist.push_back(0);
	for (int i = 0; i < num_targets; i++)
		par_optOBdist.push_back(radii[i+1]/2.0);
	par_varOBdist.push_back(0);
	for (int i = 0; i < num_targets; i++)
		par_varOBdist.push_back(1);

	par_c_hat = 0.01;
	// risk info on edges of graph G=(V,E) is generated by the following steps
	get_nodes_crds();
	get_risk_innerTrajc();
	get_risk_outerTrajc(); // 'min_risk_mat' is generated in this function.
}


PartitionScheme::~PartitionScheme() {
//	delete target_locations;
	for (int i = 0; i < num_V; i++) {
		delete[] G[i];
	}
	delete[] G;
}

vector<vector<Vertex>> PartitionScheme::get_nodes_crds() {
	// 0. all discretized boundary points are called turning point (TP)
	TP.resize(num_targets + 2);
	// 1. first turning point is the start depot
	TP[0].resize(1);
	TP[0][0] = depot1_loc;
	// 2. coordinates of all discretized points
	for (int pos = 1; pos <= num_targets; pos++) {
		TP[pos].resize(num_dstzn);
		for (int i = 0; i < num_dstzn; i++) {
			TP[pos][i].x = radii[pos] * cos(i * subarc_angle) + target_locs[pos].x;
			TP[pos][i].y = radii[pos] * sin(i * subarc_angle) + target_locs[pos].y;
		}
	}
	// 3. the last turning point is the end depot
	TP[num_targets + 1].resize(1);
	TP[num_targets + 1][0] = depot2_loc;

	return TP;
}

/*
As risk&reward functions have rotatory symmetry,  we select one turning point
and calculate its risk from itself to every other turning point in the each region.
Keep it in "risk_innerpath"
*/
void PartitionScheme::get_risk_innerTrajc() {
	int s, i, j, idx_row, idx_col, flag;
	vector<double> risk_innerpath(num_dstzn, -1);
	double temp_risk_lineseg, temp_reward_lineseg;

	for (s = 1; s <= num_targets; s++) {
		for (i = 0; i < num_dstzn; i++) {
			if (i == 0) {
			//	risk_innerpath[i] = 0.0;
				risk_innerpath[i] = INF;
			}
			else {
			//	temp_risk_lineseg = get_risk_innerTrajc(TP[s][0], TP[s][i], s);
			//	temp_reward_lineseg = get_reward_innerTrajc(TP[s][0], TP[s][i], s);

				// test the correctness by using euclidean distance
				temp_risk_lineseg = get_lineSeg_len(TP[s][0], TP[s][i]);  
				temp_reward_lineseg = get_lineSeg_len(TP[s][0], TP[s][i]);
				if (temp_risk_lineseg < risk_thold[s] && temp_reward_lineseg >= bdg_rewards[s]) {
					risk_innerpath[i] = temp_risk_lineseg; 	// find an admissible inner path
				}
				else {
					risk_innerpath[i] = INF;
				}
			}
		}
		// write in matrix G
		idx_row = (s - 1) * 2 * num_dstzn + 1;
		idx_col = (s - 1) * 2 * num_dstzn + num_dstzn + 1;
		for (i = 0; i < num_dstzn; i++) {
			for (j = 0; j < num_dstzn; j++) { 
				flag = (num_dstzn - i + j) % num_dstzn;
				G[idx_row + i][idx_col + j] = risk_innerpath[flag];
			}
		}
	}
}

double PartitionScheme::get_reward_innerTrajc(Vertex entry, Vertex exit, int i) {
	double l1 = get_lineSeg_len(entry, exit) / 2.0;
	double l2;
	double num_intervals = 200.0; // for half of the chord
	double len_interval = l1 / num_intervals;
	double total_reward = 0.0;

	double d;
	for (int j = 0; j < num_intervals; j++) {
		l2 = sqrt(radii[i] * radii[i] - l1 *l1);
		d = sqrt(l2 *l2 + pow(len_interval * j, 2));
		double temp = exp(-pow(d - par_optOBdist[i], 2) / (2 * pow(par_varOBdist[i], 2)));
		total_reward = total_reward + temp * len_interval;
	}
	total_reward = 2 * total_reward;
	return total_reward;
}

double PartitionScheme::get_risk_innerTrajc(Vertex entry, Vertex exit, int i) {
	double l1 = get_lineSeg_len(entry, exit) / 2.0;
	double l2 = sqrt(radii[i] * radii[i] - l1* l1);
	double val_risk = inner_risk_function(l1, l2, i);
	return val_risk;
}

/* Generate the accumulated risk on outer paths*/
void PartitionScheme::get_risk_outerTrajc() {
	int s, t, i, j, idxmat_1, idxmat_2, idxmat_3, idxmat_4;
	double val_risk, val_min_risk;
	// i) start depot <-> boundary t
	for (t = 1; t <= num_targets; t++) {
		idxmat_1 = (t - 1) * 2 * num_dstzn + 1; 
		val_min_risk = INF;
		for (i = 0; i < num_dstzn; i++) {
		//	val_risk = get_risk_outerTrajc(TP[t][i], depot1_loc);
			val_risk = get_lineSeg_len(TP[t][i], depot1_loc); // for testing 
			G[0][idxmat_1 + i] = val_risk;
			if (val_risk < val_min_risk) {
				val_min_risk = val_risk;
			}
		}
		min_risk_mat[0][t] = val_min_risk;
		min_risk_mat[t][0] = val_min_risk;
	}
	// ii) boundary s <-> boudary t
	for (s = 1; s <= num_targets; s++) {
		idxmat_1 = (s - 1) * 2 * num_dstzn + num_dstzn + 1; // vertices on boundary s acting as exits
		idxmat_3 = (s - 1) * 2 * num_dstzn + 1; // vertices on boundary s acting as entries
		for (t = 1; t <= num_targets && s != t; t++) {
			idxmat_2 = (t - 1) * 2 * num_dstzn + 1; // vertices on boundary t acted as entries
			idxmat_4 = (t - 1) * 2 * num_dstzn + num_dstzn + 1; // vertices on boundary t acted as exits
			val_min_risk = INF;
			for (i = 0; i < num_dstzn; i++) {
				for (j = 0; j < num_dstzn; j++) {
				//	val_risk = get_risk_outerTrajc(TP[s][i], TP[t][j]);
					val_risk = get_lineSeg_len(TP[s][i], TP[t][j]); // for testing
					G[idxmat_1 + i][idxmat_2 + j] = val_risk;
					G[idxmat_4 + j][idxmat_3 + i] = val_risk;
					if (val_risk < val_min_risk) {
						val_min_risk = val_risk;
					}
				}
			}
			min_risk_mat[s][t] = val_min_risk;
			min_risk_mat[t][s] = val_min_risk;
		}
	}
	// iii) boundary s <->  end depot
	for (s = 1; s <= num_targets; s++) {
		idxmat_2 = (s - 1) * 2 * num_dstzn + num_dstzn + 1;
		val_min_risk = INF;
		for (i = 0; i < num_dstzn; i++) {
		//	val_risk = get_risk_outerTrajc(TP[t][i], depot2_loc);
			val_risk = get_lineSeg_len(TP[s][i], depot2_loc); // for testing 
			G[idxmat_2 + i][num_V - 1] = val_risk;
			if (val_risk < val_min_risk) {
				val_min_risk = val_risk;
			}
		}
		min_risk_mat[s][num_targets + 1] = val_min_risk;
		min_risk_mat[num_targets + 1][s] = val_min_risk;
	}

	/*
	ofstream file_G("./matrixG.txt");
	for (int i = 0; i < num_V; i++) {
		for (int j = 0; j < num_V; j++) {
			file_G << G[i][j] << '\t';
		}
		file_G << '\n';
	}
	*/
	

	// ------------  Step 2: subtract the min_risk_mat[s][t] from G[i][j]   -------------
	for (t = 1; t <= num_targets; t++) {
		idxmat_1 = (t - 1) * 2 * num_dstzn + 1;
		idxmat_2 = (t - 1) * 2 * num_dstzn + num_dstzn + 1;
		for (i = 0; i < num_dstzn; i++) {
			G[0][idxmat_1 + i] = G[0][idxmat_1 + i] - min_risk_mat[0][t];
		}
	}
	for (s = 1; s <= num_targets; s++) {
		idxmat_1 = (s - 1) * 2 * num_dstzn + num_dstzn + 1; // vertices worked as exits in boundary s
		idxmat_3 = (s - 1) * 2 * num_dstzn + 1; // vertices worked as entries in boundary s
		for (t = 1; t <= num_targets && s != t; t++) {
			idxmat_2 = (t - 1) * 2 * num_dstzn + 1; // vertices worked as entries on boundary t
			idxmat_4 = (t - 1) * 2 * num_dstzn + num_dstzn + 1; // vertices worked as exits on boundary t
			for (i = 0; i < num_dstzn; i++) {
				for (j = 0; j < num_dstzn; j++) {
					G[idxmat_1 + i][idxmat_2 + j] = G[idxmat_1 + i][idxmat_2 + j] - min_risk_mat[s][t];
					G[idxmat_4 + j][idxmat_3 + i] = G[idxmat_4 + j][idxmat_3 + i] - min_risk_mat[t][s];
				}
			}
		}
	}
	for (s = 1; s <= num_targets; s++) {
		idxmat_2 = (s - 1) * 2 * num_dstzn + num_dstzn + 1;
		for (i = 0; i < num_dstzn; i++) {
			G[idxmat_2 + i][num_V - 1] = G[idxmat_2 + i][num_V - 1] - min_risk_mat[s][num_targets + 1];
		}
	}

	for (s = 0; s <= num_targets + 1; s++) {
		min_risk_mat[s][s] = 0.0;
	}
	min_risk_mat[0][num_targets + 1] = 0.0;
	min_risk_mat[num_targets + 1][0] = 0.0;

	/*
	cout << "\n";
	for (i = 0; i < num_V; i++) {
		for (j = 0; j < num_V; j++) {
			cout << G[i][j] << "  ";
		}
		cout << "\n";
	}

	cout << "\n";
	for (i = 0; i <= num_targets + 1; i++) {
		for (j = 0; j <= num_targets + 1; j++) {
			cout << min_risk_mat[i][j] << "  ";
		}
		cout << "\n";
	}
	*/
}


double PartitionScheme::get_risk_outerTrajc(Vertex v, Vertex u) {
	double total_risk = outer_risk_function(get_lineSeg_len(v, u)); // total risk over the line segment
	for (int i = 1; i <= num_targets; i++) { // additional possible risk caused when passing through some region(s)
		tuple<bool, double, double> result = is_intersected(v, u, i);
		if (get<0>(result) == true) {
			total_risk = total_risk + get<1>(result) - get<2>(result);
		}
	}
	return total_risk;
}


tuple<bool, double,double> PartitionScheme::is_intersected(Vertex v, Vertex u, int i) {

	// line segment v -> u, center is denoted as o;
	Vertex o = { target_locs[i].x, target_locs[i].y };
	myVector v_2_o(v, o); // vector v->o
	myVector v_2_u(v, u); // vector v->u
	double len_v_2_o = v_2_o.get_vecLen();
	double len_v_2_u = v_2_u.get_vecLen();
	double angle = acos(dot_product(v_2_o, v_2_u) / (len_v_2_o * len_v_2_u));
	double l2 = len_v_2_o * sin(angle); // vertical distance from center o to the chord
	double l1; // half length of the chord
	double risk_vu_on_Ui; // risk of path vu when passing through target region Ui
	if (angle < (PI / 2.0) && len_v_2_o * cos(angle)+ 0.0000001 < len_v_2_u && l2 + 0.0000001 < radii[i]) {
		l1 = sqrt(radii[i] * radii[i] - l2 * l2);
		risk_vu_on_Ui = inner_risk_function(l1, l2, i);
		return make_tuple(true, risk_vu_on_Ui, outer_risk_function(2 * l1));
	}
	else {
		return make_tuple(false, NULL, NULL);
	}
}

double PartitionScheme::inner_risk_function(double l1, double l2, int i) {
	double T = 2 * par_c1[i] * atan(l1 / l2) / l2;
	return T;
}

double PartitionScheme::outer_risk_function(double dist) {
	return par_c_hat * dist;
}

/*  supporting functions */
double PartitionScheme::dot_product(myVector vec1, myVector vec2) {
	return vec1.get_x() * vec2.get_x() + vec1.get_y() * vec2.get_y();
}

double PartitionScheme::get_lineSeg_len(Vertex v, Vertex u) {
	return sqrt(pow(u.x - v.x, 2) + pow(u.y - v.y, 2));
}


void PartitionScheme::solve_shortestpath(vector<vector<double>> & SDS, vector<int> & seq) {

	int idxmat_1, idxmat_2, idxmat_3;
	vector<double> entry_dist(num_dstzn, INF);
	vector<double> exit_dist(num_dstzn, INF);
	vector<double> INF_dist(num_dstzn, INF);
	double dist_endDepot = INF;

	double dist = 0.0;
	// i) SDS[start depot] = 0.0
	SDS[0].resize(1);
	SDS[0][0] = 0.0;

	// ii) SDS[*] of boundary points of the first visited target
	int pos = 1;
	int idx_circle = seq[pos];
	SDS[1].resize(num_dstzn);
	idxmat_1 = (idx_circle - 1) * 2 * num_dstzn + 1; // index of 1st entry at 1st target visited
	for (int i = 0; i < num_dstzn; i++)	{
		entry_dist[i] = G[0][idxmat_1 + i];
		SDS[1][i] = entry_dist[i];
	}
	idxmat_2 = (idx_circle - 1) * 2 * num_dstzn + num_dstzn + 1;  // index of 1st exit at 1st target visited
	for (int i = 0; i < num_dstzn; i++) {
		for (int j = 0; j < num_dstzn; j++) {
			dist = entry_dist[i] + G[idxmat_1 + i][idxmat_2 + j];
			if (dist < exit_dist[j])
				exit_dist[j] = dist;
		}
	}
	entry_dist = INF_dist;

	// iii) nodes from seq[2] to seq[num_circles]
	int idx_fr_circle = 0;
	int idx_lat_circle = 0;
	for (pos = 2; pos <= num_targets; pos++) {
		SDS[pos].resize(num_dstzn);
		idx_fr_circle = seq[pos - 1]; 
		idx_lat_circle = seq[pos]; 
		idxmat_1 = (idx_fr_circle - 1) * 2 * num_dstzn + num_dstzn + 1;
		idxmat_2 = (idx_lat_circle - 1) * 2 * num_dstzn + 1;
		for (int i = 0; i < num_dstzn; i++) {
			for (int j = 0; j < num_dstzn; j++) {
				dist = exit_dist[i] + G[idxmat_1 + i][idxmat_2 + j];
				if (dist < entry_dist[j])
					entry_dist[j] = dist;
			}
		}
		for (int i = 0; i < num_dstzn; i++) {
			SDS[pos][i] = entry_dist[i];
		}
		exit_dist = INF_dist;
		idxmat_3 = (idx_lat_circle - 1) * 2 * num_dstzn + num_dstzn + 1;
		// update exits
		for (int i = 0; i < num_dstzn; i++) {
			for (int j = 0; j < num_dstzn; j++) {
				dist = entry_dist[i] + G[idxmat_2 + i][idxmat_3 + j];
				if (dist < exit_dist[j])
					exit_dist[j] = dist;
			}
		}
		entry_dist = INF_dist;
	}
	// iv) the last target <-> the end depot
	idx_circle = seq[num_targets];
	SDS[num_targets + 1].resize(1);
	idxmat_1 = (idx_circle - 1) * 2 * num_dstzn + num_dstzn + 1;
	for (int i = 0; i < num_dstzn; i++) {
		dist = exit_dist[i] + G[idxmat_1 + i][num_V - 1];
		if (dist < dist_endDepot)
			dist_endDepot = dist;
	}
	SDS[num_targets + 1][0] = dist_endDepot;
}
