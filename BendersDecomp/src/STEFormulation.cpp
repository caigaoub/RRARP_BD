#include "STEFormulation.h"
#define INF numeric_limits<double>::infinity()

void STEFormulation::add_dualformul(DualFormulation* dual_){
	this->_formul_dual = dual_;
}

void STEFormulation::add_SuperCutformul(SuperCutFormulation* supercut_){
	this->_formul_supercut = supercut_;
}


/**
 	serve the master formulation with Gurobi model and the underlying network 
 	@	
*/
void STEFormulation::build_formul(GRBModel* model_MP_, PartitionScheme* network_) {
	this->_partition = network_;
	this->_size_var_y = _partition->_dataset->_nb_targets + 2;
	this->_model = model_MP_;

	/* set up model parameters */
	_model->getEnv().set(GRB_IntParam_LazyConstraints, 1);
	_model->getEnv().set(GRB_IntParam_PreCrush, 1);
	_model->getEnv().set(GRB_IntParam_NumericFocus, 1);
	_model->getEnv().set(GRB_DoubleParam_TimeLimit, 7200);

	// Step 1: create variables for master problem 
	_var_y = new GRBVar*[_size_var_y];
	for(int i = 0; i < _size_var_y; i++) {
		_var_y[i] = new GRBVar[_size_var_y];
	}
	for(int i = 0; i <_size_var_y; i++) {
		for(int j = 0; j < _size_var_y; j++) {
			_var_y[i][j] = _model->addVar(0.0, 1.0, _partition->_min_risk_tars[i][j], GRB_BINARY, "y_" + itos(i) + "," + itos(j));
		}
	}
	_model->update();
	_var_v = new GRBVar; // artifical variable v
	*_var_v = _model->addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "v");
	_model->update();

	// Step 2: add constraints: degree constraints and reverse-sequence-avoiding constraints
	GRBLinExpr expr1, expr2, exprA = 0, exprB = 0;
	for (int i = 0; i < _size_var_y; i++) {
		expr1 = 0;
		expr2 = 0;
		for (int j = 0; j < _size_var_y; j++) {
			expr1 += _var_y[i][j];
			expr2 += _var_y[j][i];
		}
		_model->addConstr(expr1 == 1, "deg1_row" + itos(i));
		_model->addConstr(expr2 == 1, "deg1_col" + itos(i));
		_model->update();
		if (i > 0 && i < _size_var_y - 1) {
			exprA += i * (_var_y[0][i]);
			exprB += i * (_var_y[i][_size_var_y - 1]);
		}
	}
	// _model->addConstr(exprA <= exprB, "rm_inv_seq");
	_model->update();

	_var_y[_size_var_y - 1][0].set(GRB_DoubleAttr_LB, 1);
	for (int i = 0; i < _size_var_y; i++) {
		_var_y[i][i].set(GRB_DoubleAttr_UB, 0);
	}
	// _model->write("./ret/model.lp");
}
/**
	solve master problem(TSP)
*/
pair<double, double> STEFormulation::solve_formul_wCB(int which_cut) {
	try {
		BendersCuts * cb = new BendersCuts(_var_y, _var_v, _partition, _formul_dual, _formul_supercut, which_cut);
		_model->setCallback(cb);
		_model->optimize();
		if (_model->get(GRB_IntAttr_Status) == GRB_OPTIMAL|| _model->get(GRB_IntAttr_Status) == GRB_TIME_LIMIT) {
			_status = 0;
			double obj_val = _model->get(GRB_DoubleAttr_ObjVal);
			double v_val = (*_var_v).get(GRB_DoubleAttr_X);
			_total_nb_Benders_cuts = cb->get_nb_Benders_cuts();
			_total_nb_subtour_cuts = cb->get_nb_subtour_cuts();
			cout << "**************** ***** ****************"<< endl;
			cout << "**************** ***** ****************"<< endl;
			cout << "**************** ***** ****************"<< endl;
			cout << "====>> objective value: " << obj_val << endl;
			cout << "====>> v value: " << v_val  << endl;
			cout << "====>> nb of Benders cuts(not including user cut): " << _total_nb_Benders_cuts << endl;
			cout << "====>> nb of Benders cuts(user cuts): " <<  _total_nb_user_cuts << endl;
			cout << "====>> nb of subtour cuts: " << _total_nb_subtour_cuts << endl;
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
pair<double, double> STEFormulation::solve_formul_woCB() {
	try {
		_model->optimize();
		if (_model->get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
			_status = 0;
			double	obj_val = _model->get(GRB_DoubleAttr_ObjVal);
			double	v_val = (*_var_v).get(GRB_DoubleAttr_X);
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
	// cout << " flag ***" << endl;
	for (int i = 0; i < _size_var_y; i++) {
		for (int j = 0; j < _size_var_y; j++) {
			sol[i][j] = _var_y[i][j].get(GRB_DoubleAttr_X);
		}
	}
	// cout << "end flag======" << endl;
}

double STEFormulation::add_USER_cuts(double** y_sol) {
	_formul_dual->set_objective(y_sol);
	double obj_dual = _formul_dual->solve();
	GRBLinExpr expr = 0;
	_formul_dual->get_Benders_user_cut(expr, _var_y);

	_model->addConstr(expr <= *_var_v, "Fischeti-cut");
	_model->update();
	_total_nb_user_cuts++;

	return obj_dual;
}


/* add subtour elimination constraints (SEC) at the root */
pair<bool,int> STEFormulation::add_SECs(double** y_sol) {
	bool is_disconnected = false; // one whole graph in default
	int new_subtour_cuts = 0;
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
					expr += _var_y[tour[start + j]][tour[start + k]] + _var_y[tour[start + k]][tour[start + j]];
				}
			}
			start += size;
			_model->addConstr(expr <= size - 1, "SEC1_Cut");
			new_subtour_cuts++;
			_model->update();
			_total_nb_subtour_cuts++;
		}
		is_disconnected = true;
	}
	else {
		for (int node = 0; node < _size_var_y; node++) { // ii) if cutting points exist in the solution
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
							expr += _var_y[tour2[start + j]][tour2[start + k]] + _var_y[tour2[start + k]][tour[start + j]];
						}
					}
					start += size;
					_model->addConstr(expr <= size - 1, "SEC2_Cut");
					new_subtour_cuts++;
					_model->update();
					_total_nb_subtour_cuts++;
				}
			}
			else { // iii)  check global min cut
				 /* y_sym is same solution as y, serving as input for finding the global min-cut */
				int m = 0; // number of weighted edges in the solution
				double ** y_sym = new double*[_size_var_y];
				for (int i = 0; i < _size_var_y; i++) {
					y_sym[i] = new double[_size_var_y];
				}
				for (int i = 0; i < _size_var_y; i++) {
					for (int j = 0; j < _size_var_y; j++) {
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
				m = _size_var_y * _size_var_y - m;
				GlobalMC gmc(_size_var_y, m, y_sym);
				tuple<double, unordered_set<Edge*>, vector<int>> res = gmc.GlobalMinCut();
				if (get<0>(res) < 2.0) {
					vector<int> S = get<2>(res);
					int lens = S.size();
					GRBLinExpr expr = 0;
					for (int j = 0; j < lens; j++) {
						for (int k = j + 1; k < lens; k++) {
							expr += _var_y[tour[j]][tour[k]] + _var_y[tour[k]][tour[j]];
						}
					}
					_model->addConstr(expr <= lens - 1, "SEC3_Cut");
					new_subtour_cuts++;
					_model->update();
					_total_nb_subtour_cuts++;
				}

				for (int i = 0; i < _size_var_y; i++)
					delete[] y_sym[i];
				delete[] y_sym;
			}
		}
	}
	return make_pair(is_disconnected,new_subtour_cuts);
}


/* Function: list all subtours by breadth first search */
void STEFormulation::check_subcomponents(double** sol, vector<int>& tour, int& num_subcompts, vector<int>& size_subcompts) {
	int i, node, cur_node;
	bool* seen = new bool[_size_var_y];
	for (i = 0; i < _size_var_y; i++) {
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
			for (i = 0; i < _size_var_y; i++) {
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
		if (num_visited_cursubtour == _size_var_y) { // only one component exists and all nodes are visited
			break;
		}
		else if (total_num_visited == _size_var_y) { // multiple components exist and all nodes are visited
			break;
		}
		for (node = 0; node < _size_var_y; node++) { // search one unvisited target and put it in the queue
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
	vector<bool> seen(_size_var_y, false);
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
			for (i = 0; i < _size_var_y; i++) {
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
		if (cur_subtour_visited == _size_var_y - 1) {
			break;
		}
		else if (total_num_visited == _size_var_y - 1) {
			break;
		}
		for (node = 0; node < _size_var_y; node++) {
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

void STEFormulation::set_vars_integer() {
	for (int i = 0; i < _size_var_y; i++) {
		for (int j = 0; j < _size_var_y; j++) {
			_var_y[i][j].set(GRB_CharAttr_VType, GRB_INTEGER);
		}
	}
	_model->update();
}

void STEFormulation::set_vars_continuous() {
	for (int i = 0; i < _size_var_y; i++) {
		for (int j = 0; j < _size_var_y; j++) {
			_var_y[i][j].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
		}
	}
	_model->update();
}

void STEFormulation::print_solution(GRBModel *model) {
	try {
		if (model->get(GRB_IntAttr_SolCount) > 0) {
			int i, j;
			double **sol = new double*[_size_var_y];
			for (i = 0; i < _size_var_y; i++) {
				sol[i] = new double[_size_var_y];
			}
			for (i = 0; i < _size_var_y; i++) {
				for (j = 0; j < _size_var_y; j++) {
					sol[i][j] = _var_y[i][j].get(GRB_DoubleAttr_X);
				}
			}
			cout << "====>> optimal TSP Solution matrix: " << '\n';
			for (i = 0; i < _size_var_y; i++) {
				for (j = 0; j < _size_var_y; j++) {
					if(abs(sol[i][j]) < 1e-6){
							cout << 0 << "  ";
					}
					if(abs(sol[i][j]-1) < 1e-6){
						cout << 1 << "  ";
					}
					// cout << sol[i][j] << "   ";	
				}
				cout << endl;
			}
			int len;
			int *tour2 = new int[_size_var_y];
			BendersCuts::findsubtour(_size_var_y, sol, &len, tour2);
			cout << "====>> current best visiting sequence: ";
			for (i = 0; i < len; i++) {
				cout << tour2[i] << "  ";
			}
			cout << '\n';
			for (i = 0; i < _size_var_y; i++)
				delete[] sol[i];
			delete[] sol;
		}
	}
	catch (const GRBException& ex) {
		cout << "Error number: " << ex.getErrorCode() << endl;
		cout << ex.getMessage() << endl;
	}
	catch (...) {
		cerr << "Error" << endl;
	}
}

STEFormulation::~STEFormulation(){
	for(int i = 0; i < _size_var_y; i++)
		delete[] _var_y[i];
	delete[] _var_y;
	delete _var_v;

}