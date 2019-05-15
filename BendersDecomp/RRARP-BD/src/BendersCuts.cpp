#include "BendersCuts.h"
#include "SubtourCuts.h"
#include <tuple>
#include "CoefReduction.h"

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

	this->fseq = new vector<int>(N, 0);
	flag = 1;

	num_Benders_cuts = 0;
	num_subtour_cuts = 0;
	count = 0;
//	evn_CoefRedc = new GRBEnv();
//	model_CoefRedc = new GRBModel(*evn_CoefRedc);
//	model_CoefRedc->getEnv().set(GRB_IntParam_OutputFlag, 0);
}

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


void BendersCuts::callback() {
	try {
		if (where == GRB_CB_MIPSOL) {
			int i;

			double **y_sol = new double*[N];
			for (i = 0; i < N; i++) {
				y_sol[i] = new double[N];
				y_sol[i] = getSolution(y[i], N);
			}

			vector<int> tour;
			int num_subcompts = 0;
			vector<int> size_subcompts;
			check_subcomponents(y_sol, tour, num_subcompts, size_subcompts);

			if (num_subcompts > 1) {
				GRBLinExpr expr = 0;
				int size;
				// if the graph has only two components, add just one constraints
				int num_subtourelim_constraints = ((num_subcompts == 2) ? 1 : num_subcompts);
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
					addLazy(expr <= size - 1);
					num_subtour_cuts++;
				}
			}
			else {

				int *tour2 = new int[N];

				int i, len;
				findsubtour(N, y_sol, &len, tour2);

				// add Benders optimality cuts by solving shortest path problem
				vector<int> fseq2(N);
				for (i = 1; i < len; i++) {
					fseq2.at(i) = tour2[i];
				}

				if (!is_inSeqPool(fseq2)) {
					SeqPool.push_back(fseq2);
			//		print_sequence(&fseq2);

					SDS = new vector<vector<double>>(num_targets + 2);
					PS->solve_shortestpath(*SDS, fseq2);

					expr = 0;
					expr = generate_Benderscut_SP(&fseq2);

					addLazy(expr >= 0);
					num_Benders_cuts++;
				}

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

			for (i = 0; i < N; i++)
				delete[] y_sol[i];
			delete[] y_sol;
		}

	//	if (where == GRB_CB_MIPNODE && getDoubleInfo(GRB_CB_MIPNODE_SOLCNT) < 2 )

	}
	catch (GRBException e) {
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during callback" << endl;
	}


}


// generate Benders optimality cuts by solving shortest path problem
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
		//	coef += max(0.0, (*SDS)[to][i] - dist);
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
	//	expr += coef * y[0][idx_circle];
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
				//	coef += max(0.0, (*SDS)[to][j] - (*SDS)[from][i] - dist);
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

		//	expr += coef * y[circle_from][circle_to];

		}
	}
	// (3) node weights from circle to sink (or we can call source)
	for (to = 1; to <= num_targets - 1; to++) {
		coef = 0.0;
		idx_circle = fseq->at(to);
		idxmat_1 = (idx_circle - 1) * 2 * num_dstzn + num_dstzn + 1;
		for (i = 0; i < num_dstzn; i++) {
			dist = G[idxmat_1 + i][0];
		//	coef += max(0.0, (*SDS)[num_targets + 1][0] - (*SDS)[to][i] - dist);
			coef = max(coef, max(0.0, (*SDS)[num_targets + 1][0] - (*SDS)[to][i] - dist));
			if (coef >= sd) {
				coef = sd;
				break;
			}
		}

	//	expr += coef * y[idx_circle][num_targets + 1];
		if (coef < smallest_coef) {
			smallest_coef = coef;
		}
		if (coef > 0.000001)
			CoefSet.push_back(make_tuple(idx_circle, num_targets + 1, coef));
	}
	/*
	// generate the cut
	double delta = sd - smallest_coef;
	if (smallest_coef < sd * 0.5) {
		for (int i = 0; i < CoefSet.size(); i++) {
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
		for (int i = 0; i < CoefSet.size(); i++) {
		//	expr += sd * 0.5 * y[get<0>(CoefSet[i])][get<1>(CoefSet[i])];
			CoefSet.at(i) = make_tuple(get<0>(CoefSet[i]), get<1>(CoefSet[i]), sd * 0.5);
	//		cout << sd * 0.5 << "*" << "y_" << get<0>(CoefSet[i]) << get<1>(CoefSet[i]) << " + ";
		}
	//	expr = expr + *v - sd;
	//	cout << endl;
	}
	*/
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
			if (i == n) {
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


bool BendersCuts::is_inSeqPool(vector<int> & seqVar) {

	bool isIn = false;
	bool isSame;
	for (unsigned int i = 0; i < SeqPool.size(); i++) {
		isSame = true;
		for (int j = 0; j < N; j++) {
			if (SeqPool[i][j] != seqVar.at(j)) {
				isSame = false;
				break;
			}
		}
		if (isSame == true) {
			isIn = true;
		}
	}

	return isIn;

}


void BendersCuts::check_subcomponents(double** sol, vector<int>& tour, int& num_subcompts, vector<int>& size_subcompts) {
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

void BendersCuts::check_cutting_point(int cuttingPoint, double** sol, vector<int> & tour, int& num_subcompts, vector<int>& size_subcompts) {
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
