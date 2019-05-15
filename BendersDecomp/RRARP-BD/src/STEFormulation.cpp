#include "gurobi_c++.h"
#include "STEFormulation.h"
#include "DualFormulation.h"
#include "GlobalMC.h"
#define INF numeric_limits<double>::infinity()

STEFormulation::STEFormulation(GRBModel* model_TSP, PartitionScheme* PS, DualFormulation * dl) {
	this->PS = PS;
	this->N = PS->get_num_targets() + 2; // In TSP model, both the first and last point represent the depot at the same time
	min_risk_mat = PS->get_min_risk_mat();
	num_dstzn = PS->get_num_dstzn();
	this->model = model_TSP;
	DL = dl;
	// set up the model
	model->getEnv().set(GRB_IntParam_LazyConstraints, 1);
	model->getEnv().set(GRB_IntParam_PreCrush, 1);
	model->getEnv().set(GRB_IntParam_NumericFocus, 1);

	/*----------------------------------------------------------------------------------------------*/
	// create variables
	int i, j;
	y = new GRBVar*[N];
	for (i = 0; i < N; i++) {
		y[i] = new GRBVar[N];
	}

	// min_risk_mat: is N-1 * N-1 matrix.
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			if (i == N - 1 && j!= N-1) {
				y[i][j] = model->addVar(0.0, 1.0, min_risk_mat[0][j], GRB_BINARY, "y_" + itos(i) + "_" + itos(j));
			}
			else if (j == N - 1 &&  i!= N-1) {
				y[i][j] = model->addVar(0.0, 1.0, min_risk_mat[i][0], GRB_BINARY, "y_" + itos(i) + "_" + itos(j));
			}
			else if (j == N - 1 && i == N - 1) {
				y[i][j] = model->addVar(0.0, 1.0, min_risk_mat[0][0], GRB_BINARY, "y_" + itos(i) + "_" + itos(j));
			}
			else {
				y[i][j] = model->addVar(0.0, 1.0, min_risk_mat[i][j], GRB_BINARY, "y_" + itos(i) + "_" + itos(j));
			}
		}
	}
	model->update();


	v = new GRBVar;
	*v = model->addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "v");
	model->update();
	GRBLinExpr expr1, expr2;

	/*----------------------------------------------------------------------------------------------*/

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
		if (i > 0 && i < N - 1 ) {
			exprA += i * (y[0][i]);
			exprB += i * (y[i][N - 1]);
		}

	}
	model->addConstr(exprA <= exprB, "C_Inv_constr");
	model->update();

	y[N - 1][0].set(GRB_DoubleAttr_LB, 1);
	for (int i = 0; i < N; i++) {
		y[i][i].set(GRB_DoubleAttr_UB, 0);
	}
//	model->write("tsp.lp");
}

pair<double, double> STEFormulation::solve_IP_TSP() {

	try {
	//	model->write("tsp.lp");
		BendersCuts * cb = new BendersCuts(model, y, v, PS, DL);
		model->setCallback(cb);

		model->optimize();
		if (model->get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
			status = 0;
			double obj_val = model->get(GRB_DoubleAttr_ObjVal);
			double v_val = (*v).get(GRB_DoubleAttr_X);
			delete cb;
			return make_pair(obj_val, v_val);
		}else{
			delete cb;
		}
		num_Benders_cuts_const = cb->get_num_Benders_cuts();
		num_subtour_cuts_const = cb->get_num_subtour_cuts();
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


pair<double,double> STEFormulation::solve_LP_TSP(){
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
		return make_pair(-INF,-INF);
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
	model->addConstr(expr <= *v, "User_cut");
	model->update();
	return obj_dual;
}



bool STEFormulation::add_SEP_cuts(double** y_sol) {

	bool is_disconnected = false; // one whole graph
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
	/*
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			cout << y_sym[i][j] << '\t';
		}
		cout << '\n';
	}
	cout << "--------------------------------------------" << endl;
	*/
	m = N*N - m;

	// 1) - multiple components
	vector<int> tour;
	int num_compts;
	vector<int> size_subcompts;
	check_subcomponents(y_sym, tour, num_compts, size_subcompts);
	if (num_compts > 1) {
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
			model->addConstr(expr <= size - 1, "SEP1_Cut");
			cout << "SEP_1_Cut is added!!" << endl;
			model->update();
		}
		is_disconnected = true;
	}
	else {
		for (int node = 0; node < N; node++) {
			vector<int> tour2;
			int num_subcompts2;
			vector<int> size_subcompts2;
			check_cutting_point(node, y_sym, tour2, num_subcompts2, size_subcompts2);
			if (num_subcompts2 > 1) {
				int start = 0;
				int size;
				int num_subtour_constrs = (num_subcompts2 == 2 ? 1 : num_subcompts2);
				//	int tmp = 1;
				for (int i = 0; i < num_subtour_constrs; i++) {
					GRBLinExpr expr = 0;
					size = size_subcompts2[i];
					for (int j = 0; j < size; j++) {
						for (int k = j + 1; k < size; k++) {
							expr += y[tour2[start + j]][tour2[start + k]] + y[tour2[start + k]][tour[start + j]];
						}
					}
					start += size;
					model->addConstr(expr <= size - 1, "SEP2_Cut");
					cout << "SEP_2_Cut is added!!" << endl;
					model->update();
				}
			}
			else { //3.  check global min cut
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
					model->addConstr(expr <= lens - 1, "SEP3_Cut");
					cout << "SEP_3_Cut is added!!" << endl;
					model->update();
				}
			}
		}

	}


	for (int i = 0; i < N; i++)
		delete[] y_sym[i];
	delete[] y_sym;

	return is_disconnected;
}



void STEFormulation::check_subcomponents(double** sol, vector<int>& tour, int& num_subcompts, vector<int>& size_subcompts) {
	int i, node, cur_node;

	bool* seen = new bool[N];
	for (i = 0; i < N; i++) {
		seen[i] = false; // all targets are unvisited at first
	}

	queue<int> que;
	que.push(0); // always begin by node 0.
	seen[0] = true; // node 0 is labeled as 'seen'.


	tour.push_back(0);
	int total_num_visited = 0;
	num_subcompts = 1; // at least one component exists
	int num_visited_cursubtour; // number of visited nodes in the current subtour

	while (true) {
		num_visited_cursubtour = 1;
		while (!que.empty()) {
			// extract the first node out of queue
			cur_node = que.front();
			que.pop();
			for (i = 0; i < N; i++) {
				if (sol[cur_node][i] > 0.001 && !seen[i]) {
					que.push(i); // push those unvisited nodes into the queue
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

		for (node = 0; node < N; node++) { // if there exists unvisited target(s), start with visiting one of its node
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
	for (int i = 0; i<N; i++) {
		for (int j = 0; j < N; j++) {
			y[i][j].set(GRB_CharAttr_VType, GRB_INTEGER);
		}
	}
	model->update();
}

void STEFormulation::set_model_LP() {
	for (int i = 0; i<N; i++) {
		for (int j = 0; j < N; j++) {
			y[i][j].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
		}
	}
	model->update();
}



void STEFormulation::printSol(GRBModel *model)
// Extract solution
{
	// ************Optimal Solution **********************
	if (model->get(GRB_IntAttr_SolCount) > 0) {
		int i, j;
		double **sol = new double*[N];
		for (i = 0; i < N; i++) {
			sol[i] = new double[N];
		}

//		double tmp;
		for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++) {
				//sol[i][j] = *model->get(GRB_DoubleAttr_X, var_y[i][j], 1);
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
