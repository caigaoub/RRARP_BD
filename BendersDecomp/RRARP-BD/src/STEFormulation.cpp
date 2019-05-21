#include "gurobi_c++.h"
#include "STEFormulation.h"
#include "DualFormulation.h"
#include "GlobalMC.h"
#define INF numeric_limits<double>::infinity()

STEFormulation::STEFormulation(GRBModel* model_TSP, PartitionScheme* PS, DualFormulation * dl) {
	this->PS = PS;
	this->N = PS->get_num_targets() + 2;
	min_risk_mat = PS->get_min_risk_mat();
	num_dstzn = PS->get_num_dstzn();
	this->model = model_TSP;
	DL = dl;
	num_user_cuts = 0;
	// Step 0: set up the model
	model->getEnv().set(GRB_IntParam_LazyConstraints, 1);
	model->getEnv().set(GRB_IntParam_PreCrush, 1);
	model->getEnv().set(GRB_IntParam_NumericFocus, 1);
	// Step 1: create master problem variables
	int i, j;
	y = new GRBVar*[N];
	for (i = 0; i < N; i++) {
		y[i] = new GRBVar[N];
	}
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			y[i][j] = model->addVar(0.0, 1.0, min_risk_mat[i][j], GRB_BINARY, "y_" + itos(i) + "_" + itos(j));
		}
	}
	model->update();
	v = new GRBVar; // artifical variable v
	*v = model->addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "v");
	model->update();
	// Step 2: add constraints: degree constraints and reverse-sequence-avoiding constraints
	GRBLinExpr expr1, expr2;
	GRBLinExpr exprA = 0;
	GRBLinExpr exprB = 0;
	for (i = 0; i < N; i++) {
		expr1 = 0;
		expr2 = 0;
		for (j = 0; j < N; j++) {
			expr1 += y[i][j];
			expr2 += y[j][i];
		}
		model->addConstr(expr1 == 1, "Deg1_Row" + itos(i));
		model->addConstr(expr2 == 1, "Deg1_Col" + itos(i));
		model->update();
		if (i > 0 && i < N - 1) {
			exprA += i * (y[0][i]);
			exprB += i * (y[i][N - 1]);
		}
	}
	model->addConstr(exprA <= exprB, "C_Inv_constr");
	model->update();
	// Step 3: transfer the shortest Hamiltonian path problem to TSP where the start depot and end depot are connected 'virtually'
	y[N - 1][0].set(GRB_DoubleAttr_LB, 1);
	for (int i = 0; i < N; i++) {
		y[i][i].set(GRB_DoubleAttr_UB, 0);
	}
}

pair<double, double> STEFormulation::solve_IP_TSP() {
	try {
		BendersCuts * cb = new BendersCuts(model, y, v, PS, DL);
		model->setCallback(cb);
		model->optimize();
		if (model->get(GRB_IntAttr_Status) == GRB_OPTIMAL|| model->get(GRB_IntAttr_Status) == GRB_TIME_LIMIT) {
			status = 0;
			double obj_val = model->get(GRB_DoubleAttr_ObjVal);
			double v_val = (*v).get(GRB_DoubleAttr_X);
			num_Benders_cuts_const = cb->get_num_Benders_cuts() + num_user_cuts;
			num_subtour_cuts_const = cb->get_num_subtour_cuts();
			delete cb;
			return make_pair(obj_val, v_val);
		}
		else {
			delete cb;
		}
	}
	catch (GRBException e) {
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during optimization" << endl;
	}
	return make_pair(-INF, -INF);
}

/*
  - solve the linear relaxation of the master problem TSP
  - return (objective value of LP-TSP, value of variable v of LP-TSP);
*/
pair<double, double> STEFormulation::solve_LP_TSP() {
	try {
		model->optimize();
		if (model->get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
			status = 0;
			double	obj_val = model->get(GRB_DoubleAttr_ObjVal);
			double	v_val = (*v).get(GRB_DoubleAttr_X);
			return make_pair(obj_val, v_val);
		}
	}
	catch (GRBException e) {
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during optimization" << endl;
	}
	return make_pair(-INF, -INF);
}

void STEFormulation::get_optimal_sol(double ** sol) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			sol[i][j] = y[i][j].get(GRB_DoubleAttr_X);
		}
	}
}

double STEFormulation::add_USER_cuts(double** y_sol) {
	DL->set_objective(y_sol);
	double obj_dual = DL->solve();
	GRBLinExpr expr = 0;
	DL->get_Benders_user_cut(expr, y);
	model->addConstr(expr <= *v, "Fischeti-cut");
	model->update();
	num_user_cuts++;
	return obj_dual;
}


/* add subtour elimination constraints (SEC) at the first root */
bool STEFormulation::add_SECs(double** y_sol) {
	bool is_disconnected = false; // one whole graph
	vector<int> tour;
	int num_compts;
	vector<int> size_subcompts;
	check_subcomponents(y_sol, tour, num_compts, size_subcompts);
	if (num_compts > 1) { // i) if multiple subtours exist in the solution
		GRBLinExpr expr = 0;
		int size;
		int num_subtourelim_constraints = ((num_compts == 2) ? 1 : num_compts);
		int start = 0;
		for (int i = 0; i < num_subtourelim_constraints; i++) {
			expr = 0;
			size = size_subcompts[i];
			for (int j = 0; j < size; j++) {
				for (int k = j + 1; k < size; k++) {
					expr += y[tour[start + j]][tour[start + k]] + y[tour[start + k]][tour[start + j]];
				}
			}
			start += size;
			model->addConstr(expr <= size - 1, "SEC1_Cut");
			model->update();
			num_subtour_cuts_const++;
		}
		is_disconnected = true;
	}
	else {
		for (int node = 0; node < N; node++) { // ii) if cutting points exist in the solution
			vector<int> tour2;
			int num_subcompts2;
			vector<int> size_subcompts2;
			check_cutting_point(node, y_sol, tour2, num_subcompts2, size_subcompts2);
			if (num_subcompts2 > 1) {
				int start = 0;
				int size;
				int num_subtour_constrs = (num_subcompts2 == 2 ? 1 : num_subcompts2);
				for (int i = 0; i < num_subtour_constrs; i++) {
					GRBLinExpr expr = 0;
					size = size_subcompts2[i];
					for (int j = 0; j < size; j++) {
						for (int k = j + 1; k < size; k++) {
							expr += y[tour2[start + j]][tour2[start + k]] + y[tour2[start + k]][tour[start + j]];
						}
					}
					start += size;
					model->addConstr(expr <= size - 1, "SEC2_Cut");
					model->update();
					num_subtour_cuts_const++;
				}
			}
			else { // iii)  check global min cut
				 /* y_sym is same solution as y, serving as input for finding the global min-cut */
				int m = 0; // number of weighted edges in the solution
				double ** y_sym = new double*[N];
				for (int i = 0; i < N; i++) {
					y_sym[i] = new double[N];
				}
				for (int i = 0; i < N; i++) {
					for (int j = 0; j < N; j++) {
						if (y_sol[i][j] <= 0 && y_sol[j][i] <= 0) {
							y_sym[i][j] = -1;
							y_sym[j][i] = -1;
							m += 1;
						}
						else {
							if (y_sol[i][j] > 0) {
								y_sym[j][i] = y_sol[i][j];
								y_sym[i][j] = y_sol[i][j];
							}
							if (y_sol[j][i] > 0) {
								y_sym[j][i] = y_sol[j][i];
								y_sym[i][j] = y_sol[j][i];
							}
						}
					}
				}
				m = N*N - m;
				GlobalMC gmc(N, m, y_sym);
				tuple<double, unordered_set<Edge*>, vector<int>> res = gmc.GlobalMinCut();
				if (get<0>(res) < 2.0) {
					vector<int> S = get<2>(res);
					int lens = S.size();
					GRBLinExpr expr = 0;
					for (int j = 0; j < lens; j++) {
						for (int k = j + 1; k < lens; k++) {
							expr += y[tour[j]][tour[k]] + y[tour[k]][tour[j]];
						}
					}
					model->addConstr(expr <= lens - 1, "SEC3_Cut");
					model->update();
					num_subtour_cuts_const++;
				}

				for (int i = 0; i < N; i++)
					delete[] y_sym[i];
				delete[] y_sym;
			}
		}
	}
	return is_disconnected;
}


/* Function: list all subtours by breadth first search */
void STEFormulation::check_subcomponents(double** sol, vector<int>& tour, int& num_subcompts, vector<int>& size_subcompts) {
	int i, node, cur_node;
	bool* seen = new bool[N];
	for (i = 0; i < N; i++) {
		seen[i] = false; // all targets are unvisited at first
	}
	queue<int> que;
	que.push(0); // always start with node 0.
	seen[0] = true; // node 0 is visited.
	tour.push_back(0);
	int total_num_visited = 0;
	num_subcompts = 1; // at least one component exists
	int num_visited_cursubtour; // number of visited nodes in the current subtour
	while (true) {
		num_visited_cursubtour = 1; // when this line gets excuted, there must be a node in the queue
		while (!que.empty()) {
			cur_node = que.front();
			que.pop();
			for (i = 0; i < N; i++) {
				if (sol[cur_node][i] > 0.000001 && !seen[i]) {
					que.push(i);
					seen[i] = true;
					tour.push_back(i);
					num_visited_cursubtour++;
				}
			}
		}
		size_subcompts.push_back(num_visited_cursubtour);
		total_num_visited += num_visited_cursubtour;
		if (num_visited_cursubtour == N) { // only one component exists and all nodes are visited
			break;
		}
		else if (total_num_visited == N) { // multiple components exist and all nodes are visited
			break;
		}
		for (node = 0; node < N; node++) { // search one unvisited target and put it in the queue
			if (!seen[node]) {
				num_subcompts++;
				que.push(node);
				seen[node] = true;
				tour.push_back(node);
				break;
			}
		}
	}
	delete[] seen;
}


/*Function: search one */
void STEFormulation::check_cutting_point(int cuttingPoint, double** sol, vector<int> & tour, int& num_subcompts, vector<int>& size_subcompts) {
	int i, node, cur_node;
	vector<bool> seen(N, false);
	seen[cuttingPoint] = true;
	queue<int> que;
	if (cuttingPoint == 0) {
		que.push(1);
		tour.push_back(1);
		seen[1] = true;
	}
	else {
		que.push(0);
		tour.push_back(0);
		seen[0] = true;
	}

	int total_num_visited = 0;
	num_subcompts = 1;
	int cur_subtour_visited; // number of visited nodes

	while (true) {
		cur_subtour_visited = 1;
		while (!que.empty()) {
			// pop out the first node of the queue
			cur_node = que.front();
			que.pop();
			for (i = 0; i < N; i++) {
				if (sol[cur_node][i] > 0.001 && !seen[i]) {
					que.push(i);
					seen[i] = true;
					tour.push_back(i);
					cur_subtour_visited++;
				}
			}
		}
		size_subcompts.push_back(cur_subtour_visited);
		total_num_visited += cur_subtour_visited;
		if (cur_subtour_visited == N - 1) {
			break;
		}
		else if (total_num_visited == N - 1) {
			break;
		}
		for (node = 0; node < N; node++) {
			if (!seen[node]) {
				num_subcompts++;
				que.push(node);
				tour.push_back(node);
				seen[node] = true;
				break;
			}
		}

	}

}

void STEFormulation::set_model_MIP() {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			y[i][j].set(GRB_CharAttr_VType, GRB_INTEGER);
		}
	}
	model->update();
}

void STEFormulation::set_model_LP() {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			y[i][j].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
		}
	}
	model->update();
}

void STEFormulation::printSol(GRBModel *model) {
	try {
		if (model->get(GRB_IntAttr_SolCount) > 0) {
			int i, j;
			double **sol = new double*[N];
			for (i = 0; i < N; i++) {
				sol[i] = new double[N];
			}
			for (i = 0; i < N; i++) {
				for (j = 0; j < N; j++) {
					sol[i][j] = y[i][j].get(GRB_DoubleAttr_X);
				}
			}
			for (i = 0; i < N; i++) {
				for (j = 0; j < N; j++) {
					cout << sol[i][j] << '\t';
				}
				cout << endl;
			}
			int len;
			int *tour2 = new int[N];
			BendersCuts::findsubtour(N, sol, &len, tour2);
			cout << '\n';
			for (i = 0; i < len; i++) {
				cout << tour2[i] << '\t';
			}
			for (i = 0; i < N; i++)
				delete[] sol[i];
			delete[] sol;
		}
		for (int i = 0; i < N; i++)
			delete[] y[i];
		delete[] y;
	}
	catch (const GRBException& ex) {
		cout << "Error number: " << ex.getErrorCode() << endl;
		cout << ex.getMessage() << endl;
	}
	catch (...) {
		cerr << "Error" << endl;
	}
}
