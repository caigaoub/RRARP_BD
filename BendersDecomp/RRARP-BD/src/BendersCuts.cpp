#include "BendersCuts.h"
#include "SubtourCuts.h"
#include <tuple>
//#include "CoefReduction.h"

void print_sequence(vector<int> * fseq) {
	cout << "sequence: ";
	for (vector<int>::iterator it = fseq->begin(); it != fseq->end(); ++it) {
		cout << *it << " ";
	}
	cout << endl;
}


BendersCuts::BendersCuts(GRBModel* m_tsp, GRBVar** yVars, GRBVar* vVar, PartitionScheme* PSVar, DualFormulation* dl) {

	this->model_tsp = m_tsp;
	this->DL = dl;
	y = yVars;
	v = vVar;
	PS = PSVar;
	num_targets = PSVar->get_num_targets();
	N = num_targets + 2;
	this->num_dstzn = PSVar->get_num_dstzn();
	G = PS->get_G();

	for (int i = 0; i <= num_targets; i++)
		fseq.push_back(-1);

	num_Benders_cuts = 0;
	num_subtour_cuts = 0;


//	evn_CoefRedc = new GRBEnv();
//	model_CoefRedc = new GRBModel(*evn_CoefRedc);
//	model_CoefRedc->getEnv().set(GRB_IntParam_OutputFlag, 0);
}
/*
double BendersCuts::improve_coef(int s, int t, double beta_sink, vector<tuple<int, int, double>> & CoefSet) {

	GRBEnv * evn_CoefRedc = new GRBEnv();
	GRBModel model_CoefRedc = GRBModel(*evn_CoefRedc);
	model_CoefRedc.getEnv().set(GRB_IntParam_OutputFlag, 0);
	CoefReduction CR(&model_CoefRedc, N);
	CR.set_constraints();
	CR.set_objective(beta_sink, CoefSet);
	CR.fix_edge(s, t);
	double gain = CR.solve();
	CR.free_edge(s, t);
	return gain;

}
*/

void BendersCuts::callback() {
	try {
		if (where == GRB_CB_MIPSOL) {
			double **y_sol = new double*[N];
			for (int i = 0; i < N; i++) {
				y_sol[i] = new double[N];
				y_sol[i] = getSolution(y[i], N);
			}
			int *tour = new int[N];
			int len;
			findsubtour(N, y_sol, &len, tour);
			if (len < N) {
				GRBLinExpr expr = 0;
				for (int i = 0; i < len; i++) {
					for (int j = i + 1; j < len; j++) {
						expr += y[tour[i]][tour[j]] + y[tour[j]][tour[i]];
					}
				}
				addLazy(expr <= len - 1);
				num_subtour_cuts++;
			}
			else {
				// add Benders optimality cuts by solving shortest path problem
				for (int i = 0; i < N - 1; i++) {
					fseq.at(i) = tour[i];
				}
				SDS = new vector<vector<double>>(num_targets + 2);
				PS->solve_shortestpath(*SDS, fseq);
				expr = 0;
		//		expr = generate_Benderscut_SP(&fseq);
				expr = generate_StrongBenderscut(&fseq);
				addLazy(expr >= 0);
				num_Benders_cuts++;
				vector<int> fseq2(N);

				/*
				// Test the correctness of adding Benders cuts by solving the dual model
				DL->set_objective(y_sol);
				double dist = DL->solve();
				expr = 0;
				DL->get_Benders_user_cut(expr, y);
				addLazy(expr <= *v);
				num_Benders_cuts++;
				*/
			}

			for (int i = 0; i < N; i++)
				delete[] y_sol[i];
			delete[] y_sol;
			delete[] tour;
		}
	}
	catch (GRBException e) {
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during callback" << endl;
	}
}

/* Function:  generate Benders optimality cuts by solving shortest path problem */
GRBLinExpr BendersCuts::generate_Benderscut_SP(vector<int> * fseq) {
	int i, j, idx_circle, idxmat_1, idxmat_2;
	double coef, dist;
	double sd = (*SDS)[num_targets + 1][0]; // sink node (depot)
	vector<tuple<int, int, double>>  CoefSet;
	double smallest_coef = INFINITY;
	// (1) node weight from 0 to all nodes in each circle
	for (int to = 2; to <= num_targets; to++) {
		coef = 0.0;
		idx_circle = fseq->at(to);
		idxmat_1 = (idx_circle - 1) * 2 * num_dstzn + 1;
		for (i = 0; i < num_dstzn; i++) {
			dist = G[0][idxmat_1 + i];
			coef += max(0.0, (*SDS)[to][i] - dist);
		//	coef = max(coef, max(0.0, (*SDS)[to][i] - dist));
			if (coef >= sd) {
				coef = sd;
				break;
			}
		}
		if (coef < smallest_coef) {
			smallest_coef = coef;
		}
		CoefSet.push_back(make_tuple(0, idx_circle, coef));
	}

	// (2)  node weights from circle to circle
	int circle_from, circle_to, from, to;

	for (from = 1; from <= num_targets - 2; from++) {
		circle_from = fseq->at(from);
		idxmat_1 = (circle_from - 1) * 2 * num_dstzn + num_dstzn + 1;
		for (to = from + 2; to <= num_targets; to++) {
			circle_to = fseq->at(to);
			idxmat_2 = (circle_to - 1) * 2 * num_dstzn + 1;
			coef = 0.0;
			for (i = 0; i < num_dstzn; i++) {
				for (j = 0; j < num_dstzn; j++) {
					dist = G[idxmat_1 + i][idxmat_2 + j];
					coef += max(0.0, (*SDS)[to][j] - (*SDS)[from][i] - dist);
			//		coef = max(coef, max(0.0, (*SDS)[to][j] - (*SDS)[from][i] - dist));
				}
				if (coef >= sd) {
					coef = sd;
					break;
				}
			}
			if (coef < smallest_coef) {
				smallest_coef = coef;
			}
			if (coef > 0.000001)
				CoefSet.push_back(make_tuple(circle_from, circle_to, coef));
		}
	}
	// (3) node weights from circle to sink (or we can call source)
	for (to = 1; to <= num_targets - 1; to++) {
		coef = 0.0;
		idx_circle = fseq->at(to);
		idxmat_1 = (idx_circle - 1) * 2 * num_dstzn + num_dstzn + 1;
		for (i = 0; i < num_dstzn; i++) {
			dist = G[idxmat_1 + i][0];
			coef += max(0.0, (*SDS)[num_targets + 1][0] - (*SDS)[to][i] - dist);
		//	coef = max(coef, max(0.0, (*SDS)[num_targets + 1][0] - (*SDS)[to][i] - dist));
			if (coef >= sd) {
				coef = sd;
				break;
			}
		}
		if (coef < smallest_coef) {
			smallest_coef = coef;
		}
		if (coef > 0.000001)
			CoefSet.push_back(make_tuple(idx_circle, num_targets + 1, coef));
	}
	for (unsigned int i = 0; i < CoefSet.size(); i++) {
		expr += get<2>(CoefSet[i]) * y[get<0>(CoefSet[i])][get<1>(CoefSet[i])];
	//	cout << get<2>(CoefSet[i]) << "*" << "y_" << get<0>(CoefSet[i]) << get<1>(CoefSet[i]) << " + ";
	}
//	cout << endl;
	expr = expr + *v - sd;
	delete SDS;
	return expr;
}

GRBLinExpr BendersCuts::generate_StrongBenderscut(vector<int> * fseq) {
	int i, j, idx_circle, idxmat_1, idxmat_2;
	double coef, dist;
	double sd = (*SDS)[num_targets + 1][0]; // sink node (depot)
	vector<tuple<int, int, double>>  CoefSet;
	double smallest_coef = INFINITY;
	// (1) node weight from 0 to all nodes in each circle
	for (int to = 2; to <= num_targets; to++) {
		coef = 0.0;
		idx_circle = fseq->at(to);
		idxmat_1 = (idx_circle - 1) * 2 * num_dstzn + 1;
		for (i = 0; i < num_dstzn; i++) {
			dist = G[0][idxmat_1 + i];
			coef = max(coef, max(0.0, (*SDS)[to][i] - dist));
			if (coef >= sd) {
				coef = sd;
				break;
			}
		}
		if (coef < smallest_coef) {
			smallest_coef = coef;
		}
		CoefSet.push_back(make_tuple(0, idx_circle, coef));
	}

	// (2)  node weights from circle to circle
	int circle_from, circle_to, from, to;
	for (from = 1; from <= num_targets - 2; from++) {
		circle_from = fseq->at(from);
		idxmat_1 = (circle_from - 1) * 2 * num_dstzn + num_dstzn + 1;
		for (to = from + 2; to <= num_targets; to++) {
			circle_to = fseq->at(to);
			idxmat_2 = (circle_to - 1) * 2 * num_dstzn + 1;
			coef = 0.0;
			for (i = 0; i < num_dstzn; i++) {
				for (j = 0; j < num_dstzn; j++) {
					dist = G[idxmat_1 + i][idxmat_2 + j];
					coef = max(coef, max(0.0, (*SDS)[to][j] - (*SDS)[from][i] - dist));
				}
				if (coef >= sd) {
					coef = sd;
					break;
				}
			}
			if (coef < smallest_coef) {
				smallest_coef = coef;
			}
			if (coef > 0.000001)
				CoefSet.push_back(make_tuple(circle_from, circle_to, coef));
		}
	}
	// (3) node weights from circle to sink (or we can call source)
	for (to = 1; to <= num_targets - 1; to++) {
		coef = 0.0;
		idx_circle = fseq->at(to);
		idxmat_1 = (idx_circle - 1) * 2 * num_dstzn + num_dstzn + 1;
		for (i = 0; i < num_dstzn; i++) {
			dist = G[idxmat_1 + i][0];
			coef = max(coef, max(0.0, (*SDS)[num_targets + 1][0] - (*SDS)[to][i] - dist));
			if (coef >= sd) {
				coef = sd;
				break;
			}
		}
		if (coef < smallest_coef) {
			smallest_coef = coef;
		}
		if (coef > 0.000001)
			CoefSet.push_back(make_tuple(idx_circle, num_targets + 1, coef));
	}

	// generate the cut
	double delta = sd - smallest_coef;
	if (smallest_coef < sd * 0.5) {
		for (unsigned int i = 0; i < CoefSet.size(); i++) {
			if (get<2>(CoefSet[i]) >= delta) {
		//		expr += delta * y[get<0>(CoefSet[i])][get<1>(CoefSet[i])];
				CoefSet.at(i) = make_tuple(get<0>(CoefSet[i]), get<1>(CoefSet[i]), delta);
		//		cout << delta << "*" << "y_" << get<0>(CoefSet[i]) << get<1>(CoefSet[i]) << " + ";
			}
			else {
		//		expr += get<2>(CoefSet[i]) * y[get<0>(CoefSet[i])][get<1>(CoefSet[i])];

	//			cout << get<2>(CoefSet[i]) << "*" << "y_" << get<0>(CoefSet[i]) << get<1>(CoefSet[i]) << " + ";
			}
		}
//		expr = expr + (*v) - sd;
	//	cout << endl;
	}
	else {
		for (unsigned int i = 0; i < CoefSet.size(); i++) {
		//	expr += sd * 0.5 * y[get<0>(CoefSet[i])][get<1>(CoefSet[i])];
			CoefSet.at(i) = make_tuple(get<0>(CoefSet[i]), get<1>(CoefSet[i]), sd * 0.5);
	//		cout << sd * 0.5 << "*" << "y_" << get<0>(CoefSet[i]) << get<1>(CoefSet[i]) << " + ";
		}
	//	expr = expr + *v - sd;
	//	cout << endl;
	}

	/*
	if (count <= 1) {
		int s, t;
		double val;
		for (int i = 0; i < CoefSet.size(); i++) {
			s = get<0>(CoefSet[i]);
			t = get<1>(CoefSet[i]);
			val = improve_coef(s, t, sd, CoefSet);
			if (val < 0) {
				cout << "............." << endl;
				val = max(0.0, get<2>(CoefSet[i]) + val);
				CoefSet.at(i) = make_tuple(s, t, val);
			}
		}
		count++;
	}
	cout << endl;
	*/
	for (unsigned int i = 0; i < CoefSet.size(); i++) {
		expr += get<2>(CoefSet[i]) * y[get<0>(CoefSet[i])][get<1>(CoefSet[i])];
	//	cout << get<2>(CoefSet[i]) << "*" << "y_" << get<0>(CoefSet[i]) << get<1>(CoefSet[i]) << " + ";
	}
//	cout << endl;
	expr = expr + *v - sd;
	delete SDS;
	return expr;
}


void BendersCuts::findsubtour(int  n, double** sol, int*  tourlenP, int*  tour) {
	bool* seen = new bool[n];
	int bestind, bestlen;
	int i, node, len, start;
	for (i = 0; i < n; i++)
		seen[i] = false;
	start = 0;
	bestlen = n + 1;
	bestind = -1;
	node = 0;
	while (start < n) {
		for (node = 0; node < n; node++)
			if (!seen[node])
				break;
		if (node == n)
			break;
		for (len = 0; len < n; len++) {
			tour[start + len] = node;
			seen[node] = true;
			for (i = 0; i < n; i++) {
				if (sol[node][i] > 0.5 && !seen[i]) {
					node = i;
					break;
				}
			}
			if (i == n) { // all adj(node) are visited
				len++;
				if (len < bestlen) {
					bestlen = len;
					bestind = start;
				}
				start += len;
				break;
			}
		}
	}

	for (i = 0; i < bestlen; i++)
		tour[i] = tour[bestind + i];
	*tourlenP = bestlen;

	delete[] seen;
}

void BendersCuts::print_ySol(double ** y_sol) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			cout << y_sol[i][j] << "   ";
		}
		cout << endl;
	}

}
