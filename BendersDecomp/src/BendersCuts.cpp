#include "BendersCuts.h"
#include "SubtourCuts.h"
#include <tuple>
// #include "SuperCutFormulation.h"

void print_sequence(vector<int> * fseq) {
	// cout << " sequence: ";
	for (vector<int>::iterator it = fseq->begin(); it != fseq->end(); ++it) {
		cout << *it << " ";
	}
	cout << endl;
}
void print_sequence(vector<int> fseq) {
	// cout << " sequence: ";
	for (unsigned int i=0; i < fseq.size(); i++) {
		cout << fseq[i] << " ";
	}
	cout << endl;
}

BendersCuts::BendersCuts(GRBVar** y_, GRBVar* v_, PartitionScheme* partition_, DualFormulation* dual_, SuperCutFormulation* supercut_formul_, int which_cut){
// BendersCuts::BendersCuts(GRBVar** y_, GRBVar* v_, PartitionScheme* partition_) {
	this->_var_y = y_;
	this->_var_v = v_;
	this->_partition = partition_;
	this->_nb_targets = _partition->_dataset->_nb_targets;
	// this->_N = _nb_targets + 2;
	this->_nb_dstzn = _partition->_nb_dstzn;
	this->_G = &(_partition->_G);
	this->_formul_dual = dual_;
	this->_formul_supercut = supercut_formul_;
	this->_which_Bcut = which_cut;

	for (int i = 0; i <= _nb_targets; i++)
		_fseq.push_back(-1);

}


void BendersCuts::callback() {
	try {
		if (where == GRB_CB_MIPSOL) {
			double **y_sol = new double*[_nb_targets+2];
			for (int i = 0; i < _nb_targets+2; i++) {
				y_sol[i] = new double[_nb_targets+2];
				y_sol[i] = getSolution(_var_y[i], _nb_targets+2);
			}
			int *tour = new int[_nb_targets+2];
			int len;
			findsubtour(_nb_targets+2, y_sol, &len, tour);
			if (len < _nb_targets+2) {
				GRBLinExpr expr = 0;
				for (int i = 0; i < len; i++) {
					for (int j = i + 1; j < len; j++) {
						expr += _var_y[tour[i]][tour[j]] + _var_y[tour[j]][tour[i]];
					}
				}
				addLazy(expr <= len - 1);
				_CB_nb_subtour_cuts++;
			}
			else {
				// find a feasible sequence. Add Benders optimality cuts by solving shortest path problem
				for (int i = 0; i < _nb_targets + 1; i++) {
					_fseq.at(i) = tour[i];
				}
				// print_sequence(_fseq);
				if(_which_Bcut == 1){// Test the correctness of adding Benders cuts by solving the dual model			
					// print_ySol(y_sol);
					_formul_dual->set_objective(y_sol);
					_formul_dual->solve();
					GRBLinExpr expr = 0;
					_formul_dual->get_Benders_user_cut(expr, _var_y);
					addLazy(expr <= *_var_v);
					_CB_nb_Benders_cuts++;
				}

				if(_which_Bcut == 2){
					_SDS = new vector<vector<double>>(_nb_targets + 2);
					_partition->solve_shortestpath(_fseq, *_SDS);
					GRBLinExpr expr = 0;
					expr = generate_Benderscut_SP(&_fseq);
					addLazy(expr >= 0);
					_CB_nb_Benders_cuts++;
				}

				if(_which_Bcut == 3){
					_SDS = new vector<vector<double>>(_nb_targets + 2);
					_partition->solve_shortestpath(_fseq, *_SDS);
					GRBLinExpr expr = 0;
					expr = generate_StrongBenderscut(&_fseq, false);
					addLazy(expr >= 0);
					_CB_nb_Benders_cuts++;
				}

				if(_which_Bcut == 4){
					_SDS = new vector<vector<double>>(_nb_targets + 2);
					_partition->solve_shortestpath(_fseq, *_SDS);
					GRBLinExpr expr = 0;
					expr = generate_StrongBenderscut(&_fseq, true);
					addLazy(expr >= 0);
					_CB_nb_Benders_cuts++;
				}

			}

			for (int i = 0; i < _nb_targets+2; i++)
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
GRBLinExpr BendersCuts::generate_Benderscut_SP(vector<int> * fseq_) {
	int idx_circle, idxmat_1, idxmat_2, circle_from, circle_to;
	double coef, dist;
	double sd = (*_SDS)[_nb_targets + 1][0]; // sink node (depot)
	vector<tuple<int, int, double>>  CoefSet;
	double smallest_coef = INFINITY;
	//  node weight from 0 to all nodes in each circle
	for (int to = 2; to <= _nb_targets; to++) {
		coef = 0.0;
		idx_circle = fseq_->at(to);
		idxmat_1 = (idx_circle - 1) * 2 * _nb_dstzn + 1;
		for (int i = 0; i < _nb_dstzn; i++) {
			if((*_G)[0][idxmat_1 + i].first == true){
				dist = (*_G)[0][idxmat_1 + i].second;
				coef += max(0.0, (*_SDS)[to][i] - dist);
			//	coef = max(coef, max(0.0, (*SDS)[to][i] - dist));
				if (coef >= sd) {
					coef = sd;
					break;
				}
			}
		}
		if (coef < smallest_coef) {
			smallest_coef = coef;
		}
		CoefSet.push_back(make_tuple(0, idx_circle, coef));
	}

	// (2)  node weights from circle to circle
	// int circle_from, circle_to, from, to;
	for (int from = 1; from <= _nb_targets - 2; from++) {
		circle_from = fseq_->at(from);
		idxmat_1 = (circle_from - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
		for (int to = from + 2; to <= _nb_targets; to++) {
			circle_to = fseq_->at(to);
			idxmat_2 = (circle_to - 1) * 2 * _nb_dstzn + 1;
			coef = 0.0;
			for (int i = 0; i < _nb_dstzn; i++) {
				for (int j = 0; j < _nb_dstzn; j++) {
					if((*_G)[idxmat_1 + i][idxmat_2 + j].first == true){
						dist = (*_G)[idxmat_1 + i][idxmat_2 + j].second;
						coef += max(0.0, (*_SDS)[to][j] - (*_SDS)[from][i] - dist);
				//		coef = max(coef, max(0.0, (*SDS)[to][j] - (*SDS)[from][i] - dist));
					}
					
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
	for (int to = 1; to <= _nb_targets - 1; to++) {
		coef = 0.0;
		idx_circle = fseq_->at(to);
		idxmat_1 = (idx_circle - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
		for (int i = 0; i < _nb_dstzn; i++) {
			if((*_G)[idxmat_1 + i][0].first == true){
				dist = (*_G)[idxmat_1 + i][0].second;
				coef += max(0.0, (*_SDS)[_nb_targets + 1][0] - (*_SDS)[to][i] - dist);
			//	coef = max(coef, max(0.0, (*SDS)[_nb_targets + 1][0] - (*SDS)[to][i] - dist));
				if (coef >= sd) {
					coef = sd;
					break;
				}
			}			
		}
		if (coef < smallest_coef) {
			smallest_coef = coef;
		}
		if (coef > 0.000001)
			CoefSet.push_back(make_tuple(idx_circle, _nb_targets + 1, coef));
	}
	GRBLinExpr expr = 0;
	for (unsigned int i = 0; i < CoefSet.size(); i++) {
		expr += get<2>(CoefSet[i]) * _var_y[get<0>(CoefSet[i])][get<1>(CoefSet[i])];
	//	cout << get<2>(CoefSet[i]) << "*" << "y_" << get<0>(CoefSet[i]) << get<1>(CoefSet[i]) << " + ";
	}
//	cout << endl;
	expr = expr + *_var_v - sd;
	delete _SDS;
	return expr;
}

GRBLinExpr BendersCuts::generate_StrongBenderscut(vector<int> * fseq, bool supercut_on) {
	int i, j, idx_circle, idxmat_1, idxmat_2;
	double coef, dist;
	double sd = (*_SDS)[_nb_targets + 1][0]; // sink node (depot)
	vector<tuple<int, int, double>>  CoefSet;
	double smallest_coef = INFINITY;
	// (1) node weight from 0 to all nodes in each circle
	for (int to = 2; to <= _nb_targets; to++) {
		coef = 0.0;
		idx_circle = fseq->at(to);
		idxmat_1 = (idx_circle - 1) * 2 * _nb_dstzn + 1;
		for (i = 0; i < _nb_dstzn; i++) {
			if((*_G)[0][idxmat_1 + i].first == true){
				dist = (*_G)[0][idxmat_1 + i].second;
				coef = max(coef, max(0.0, (*_SDS)[to][i] - dist));
				if (coef >= sd) {
					coef = sd;
					break;
				}
			}
		}
		if (coef < smallest_coef) {
			smallest_coef = coef;
		}
		CoefSet.push_back(make_tuple(0, idx_circle, coef));
	}

	// (2)  node weights from circle to circle
	int circle_from, circle_to, from, to;
	for (from = 1; from <= _nb_targets - 2; from++) {
		circle_from = fseq->at(from);
		idxmat_1 = (circle_from - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
		for (to = from + 2; to <= _nb_targets; to++) {
			circle_to = fseq->at(to);
			idxmat_2 = (circle_to - 1) * 2 * _nb_dstzn + 1;
			coef = 0.0;
			for (i = 0; i < _nb_dstzn; i++) {
				for (j = 0; j < _nb_dstzn; j++) {
					if((*_G)[idxmat_1 + i][idxmat_2 + j].first == true){
						dist = (*_G)[idxmat_1 + i][idxmat_2 + j].second;
						coef = max(coef, max(0.0, (*_SDS)[to][j] - (*_SDS)[from][i] - dist));
					}					
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
	for (to = 1; to <= _nb_targets - 1; to++) {
		coef = 0.0;
		idx_circle = fseq->at(to);
		idxmat_1 = (idx_circle - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
		for (i = 0; i < _nb_dstzn; i++) {
			if((*_G)[idxmat_1 + i][_partition->_size_G -1].first == true){
				dist = (*_G)[idxmat_1 + i][_partition->_size_G -1].second;
				coef = max(coef, max(0.0, (*_SDS)[_nb_targets + 1][0] - (*_SDS)[to][i] - dist));
				if (coef >= sd) {
					coef = sd;
					break;
				}
			}
		}
		if (coef < smallest_coef) {
			smallest_coef = coef;
		}
		if (coef > 0.000001)
			CoefSet.push_back(make_tuple(idx_circle, _nb_targets + 1, coef));
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

	
	if (supercut_on) {
		int s, t;
		double gain;
		if(_idx_supercut <= _max_supercuts){
			for (unsigned int i = 0; i < CoefSet.size(); i++) {
				s = get<0>(CoefSet[i]);
				t = get<1>(CoefSet[i]);
				cout <<"====>> fix edge: " << s << " " << t << " = " << get<2>(CoefSet[i]) << endl;
				/* solve the super cut formulation */
				_formul_supercut->fix_edge(s, t);
				_formul_supercut->set_objective(sd, CoefSet);
				// if(i == 20)
				// 	_formul_supercutsupercut->_model->write("./ret/supercutformul.lp");
				double **sol2 = new double*[_formul_supercut->_size_var_x];
				for (int ii = 0; ii < _formul_supercut->_size_var_x; ii++) {
					sol2[ii] = new double[_formul_supercut->_size_var_x];
					// sol2[ii] = getSolution(_formul_supercut->_var_x[ii], _nb_targets+2);
				}
				sol2 = _formul_supercut->solve_tmp();
				gain = _formul_supercut->solve();

				_formul_supercut->free_edge(s, t);
				cout << "gain: " << gain << endl;
				/******************************TEST *********************************8*/
				cout << "===>> Sequence 1: ";
				print_sequence(fseq);
				double L1 = sd;
				cout << "====>> L1 = " << L1 << endl;
				
				double delta_1 = 0.0;
				for (unsigned int kk = 0; kk < CoefSet.size(); kk++){
					delta_1 += get<2>(CoefSet[kk]) * sol2[get<0>(CoefSet[kk])][get<1>(CoefSet[kk])];
				}
				cout << "====>> delta_1: " << delta_1 << endl;

				cout << "====>> sequence 2 matrix: " << '\n';
				for (int ki = 0; ki < _formul_supercut->_size_var_x; ki++) {
					for (int kj = 0; kj < _formul_supercut->_size_var_x; kj++) {
						if(abs(sol2[ki][kj]) < 1e-6){
								cout << 0 << "  ";
						}
						if(abs(sol2[ki][kj]-1) < 1e-6){
							cout << 1 << "  ";
						}
						// cout << sol[i][j] << "   ";	
					}
					cout << endl;
				}
				int len;
				int *tour2 = new int[_formul_supercut->_size_var_x];
				BendersCuts::findsubtour(_formul_supercut->_size_var_x, sol2, &len, tour2);
				vector<int> fseq22;
				for (int mi = 0; mi < _nb_targets + 1; mi++) {
					fseq22.push_back(tour2[mi]);
				}

				for (int mi = 0; mi < _nb_targets; mi++) {
					for(int jj= 0; jj < CoefSet.size(); jj++){
						if(get<0>(CoefSet[jj]) ==  fseq22[mi] && get<1>((CoefSet[jj])) == fseq22[mi+1]){
							cout << fseq22[mi] << " -> " << fseq22[mi+1] << " : " << get<2>(CoefSet[jj]) << endl;
						}
					}
					
				}
				cout << "==========================" << endl;
				for(int jj= 0; jj < CoefSet.size(); jj++){
					cout << get<0>(CoefSet[jj]) << " -> " << get<1>(CoefSet[jj]) << ": " << get<2>(CoefSet[jj]) << endl;
				}



				vector<vector<double>> * SDS2 = new vector<vector<double>>(_nb_targets + 2);
				_partition->solve_shortestpath(fseq22, (*SDS2));
				double L2 = (*SDS2)[_nb_targets+1][0];
				cout << "====>> Sequence 2: ";
				print_sequence(fseq22);
				// for (int ni = 0; ni < len; ni++) {
				// 	cout << tour2[ni] << "  ";
				// }
				// cout << '\n';
				cout << "====>> L2 = " << L2 << endl;

				cout << "L1 - L2 = " << L1 - L2 << endl;
				for (int ni = 0; ni < _formul_supercut->_size_var_x; ni++)
					delete[] sol2[ni];
				delete[] sol2;

				exit(0);
				// cout << i << ": " << s << "-->" << t << "   alpha= " << sd <<  "\t old coef:" << get<2>(CoefSet[i]) << " reduction: " << gain << endl;
				if (gain > 0) { 
					gain = max(0.0, get<2>(CoefSet[i]) - gain);
					CoefSet.at(i) = make_tuple(s, t, gain);
				}
			}
			_idx_supercut++;
		}
		
	}
	// cout << endl;
	
	GRBLinExpr expr = 0;
	for (unsigned int i = 0; i < CoefSet.size(); i++) {
		expr += get<2>(CoefSet[i]) * _var_y[get<0>(CoefSet[i])][get<1>(CoefSet[i])];
	//	cout << get<2>(CoefSet[i]) << "*" << "y_" << get<0>(CoefSet[i]) << get<1>(CoefSet[i]) << " + ";
	}
//	cout << endl;
	expr = expr + *_var_v - sd;
	delete _SDS;
	return expr;
}


double BendersCuts::improve_coef(int s, int t, double val_alpha_arrival, vector<tuple<int, int, double>> & CoefSet) {
	_formul_supercut->fix_edge(s, t);
	_formul_supercut->set_objective(val_alpha_arrival, CoefSet);
	_formul_supercut->_model->write("./ret/supercutformul.lp");

	double gain = _formul_supercut->solve();
	_formul_supercut->free_edge(s, t);
	return gain;
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
	for (int i = 0; i < _nb_targets+2; i++) {
		for (int j = 0; j < _nb_targets+2; j++) {
			cout << y_sol[i][j] << "   ";
		}
		cout << endl;
	}

}
