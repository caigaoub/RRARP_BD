#include "BendersCuts.h"
#include "SubtourCuts.h"
#include <tuple>
#include <assert.h>
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
	this->_var_y = y_;
	this->_var_v = v_;
	this->_partition = partition_;
	this->_nb_targets = _partition->_dataset->_nb_targets;
	this->_nb_dstzn = _partition->_nb_dstzn;
	this->_G = &(_partition->_G);
	this->_formul_dual = dual_;
	this->_formul_supercut = supercut_formul_;
	this->_which_Bcut = which_cut;

	for (int i = 0; i <= _nb_targets + 1; i++)
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
				for (int i = 0; i <= _nb_targets + 1; i++) {
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
					vector<vector<double>>  SDS(_nb_targets + 2);
					_partition->solve_shortestpath(_fseq, SDS);
					GRBLinExpr expr = 0;
					expr = generate_Benderscut_SP(&_fseq, SDS);
					addLazy(expr >= 0);
					_CB_nb_Benders_cuts++;
				}

				if(_which_Bcut == 3){
					vector<vector<double>>  SDS(_nb_targets + 2);
					_partition->solve_shortestpath(_fseq, SDS);
					GRBLinExpr expr = 0;
					expr = generate_StrongCut(&_fseq, SDS); 
					addLazy(expr >= 0);
					_CB_nb_Benders_cuts++;
				}

				if(_which_Bcut == 4){
					vector<vector<double>>  SDS(_nb_targets + 2);
					_partition->solve_shortestpath(_fseq, SDS);
					GRBLinExpr expr = 0;
					expr = generate_SuperCut(&_fseq, SDS);
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
GRBLinExpr BendersCuts::generate_Benderscut_SP(vector<int> * fseq_, vector<vector<double>> & SDS) {
	int idx_circle, idxmat_1, idxmat_2, circle_from, circle_to;
	double coef, dist;
	double sd = SDS[_nb_targets + 1][0]; // sink node (depot)
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
				coef += max(0.0, SDS[to][i] - dist);
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
						coef += max(0.0, SDS[to][j] - SDS[from][i] - dist);
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
				coef += max(0.0, SDS[_nb_targets + 1][0] - SDS[to][i] - dist);
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
		// cout << get<2>(CoefSet[i]) << "*" << "y_" << get<0>(CoefSet[i]) << get<1>(CoefSet[i]) << " + ";
	}
	// cout << endl;
	expr = expr + *_var_v - sd;
	return expr;
}

GRBLinExpr BendersCuts::generate_StrongCut(vector<int> * fseq, vector<vector<double>> & SDS) {
	int i, j, idx_circle, idxmat_1, idxmat_2;
	double coef, dist;
	double sd = SDS[_nb_targets + 1][0]; // sink node (depot)

	vector<pair<pair<int,int>, double>> Coefs;
	double smallest_coef = INFINITY;
	// (1) node weight from 0 to all nodes in each circle
	for (int to = 2; to <= _nb_targets; to++) {
		coef = 0.0;
		idx_circle = fseq->at(to);
		idxmat_1 = (idx_circle - 1) * 2 * _nb_dstzn + 1;
		for (i = 0; i < _nb_dstzn; i++) {
			if((*_G)[0][idxmat_1 + i].first == true){
				dist = (*_G)[0][idxmat_1 + i].second;
				coef = max(coef, max(0.0, SDS[to][i] - dist));
				if (coef >= sd) {
					coef = sd;
					break;
				}
			}
		}
		if (coef < smallest_coef) {
			smallest_coef = coef;
		}

		Coefs.push_back(make_pair(make_pair(0, idx_circle), coef));
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
						coef = max(coef, max(0.0, SDS[to][j] - SDS[from][i] - dist));
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
				Coefs.push_back(make_pair(make_pair(circle_from, circle_to), coef));
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
				coef = max(coef, max(0.0, SDS[_nb_targets + 1][0] - SDS[to][i] - dist));
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
			Coefs.push_back(make_pair(make_pair(idx_circle, _nb_targets + 1), coef));
	}

	// generate the cut
	double delta = sd - smallest_coef;
	if (smallest_coef < sd * 0.5) {
		for (unsigned int i=0; i< Coefs.size(); i++) {
			if ( Coefs[i].second >= delta) {
				Coefs[i].second = delta;
			}
		}
	}
	else {
		for (unsigned int i=0; i< Coefs.size(); i++) {
			Coefs[i].second = sd * 0.5;
		}
	}

	GRBLinExpr expr2 = 0;
	// cout << sd;
	for (unsigned int i = 0; i < Coefs.size(); i++) {
		expr2 += Coefs[i].second * _var_y[Coefs[i].first.first][Coefs[i].first.second];
		// cout << " - " << Coefs[i].second << "*" << "x_" << Coefs[i].first.first <<','<< Coefs[i].first.second;
	}
	// cout << endl;
	expr2 = expr2 + (*_var_v) - sd;
	return expr2;
}




GRBLinExpr BendersCuts::generate_SuperCut(vector<int> * fseq, vector<vector<double>> & SDS) {
	int i, j, idx_circle, idxmat_1, idxmat_2;
	double coef, dist;
	double sd = SDS[_nb_targets + 1][0]; // sink node (depot)

	vector<pair<pair<int,int>, double>> Coefs;
	double smallest_coef = INFINITY;
	// (1) node weight from 0 to all nodes in each circle
	for (int to = 2; to <= _nb_targets; to++) {
		coef = 0.0;
		idx_circle = fseq->at(to);
		idxmat_1 = (idx_circle - 1) * 2 * _nb_dstzn + 1;
		for (i = 0; i < _nb_dstzn; i++) {
			if((*_G)[0][idxmat_1 + i].first == true){
				dist = (*_G)[0][idxmat_1 + i].second;
				coef = max(coef, max(0.0, SDS[to][i] - dist));
				if (coef >= sd) {
					coef = sd;
					break;
				}
			}
		}
		if (coef < smallest_coef) {
			smallest_coef = coef;
		}

		Coefs.push_back(make_pair(make_pair(0, idx_circle), coef));
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
						coef = max(coef, max(0.0, SDS[to][j] - SDS[from][i] - dist));
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
				Coefs.push_back(make_pair(make_pair(circle_from, circle_to), coef));
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
				coef = max(coef, max(0.0, SDS[_nb_targets + 1][0] - SDS[to][i] - dist));
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
			Coefs.push_back(make_pair(make_pair(idx_circle, _nb_targets + 1), coef));
	}

	// generate the cut
	double delta = sd - smallest_coef;
	if (smallest_coef < sd * 0.5) {
		for (unsigned int i=0; i< Coefs.size(); i++) {
			if ( Coefs[i].second >= delta) {
				Coefs[i].second = delta;
			}
		}
	}
	else {
		for (unsigned int i=0; i< Coefs.size(); i++) {
			Coefs[i].second = sd * 0.5;
		}
	}

	// cout << "-> " << supercut_on << endl;
	// if(false){
		// cout << "doing nothing " << endl;
	// }

	if(Coefs.size() != 0){
		sort(Coefs.begin(), Coefs.end(), [](const pair<pair<int,int>,double> & a, const pair<pair<int,int>,double> & b) -> bool{return a.second> b.second;});
	}
		
	int s, t;
	double gain = 0;	
	// cout << "coef size " << Coefs.size() << endl;
	// for(unsigned int i = 0; i < Coefs.size(); i++){
	// 	cout << Coefs[i].second << ": " << Coefs[i].first.first << ", " << Coefs[i].first.second << " **** ";
	// }

	if(_idx_supercut <= _max_supercuts){
		for (unsigned int i = 0; i < min((double)Coefs.size()/3.0,5.0); i++) {
			// cout << "cut: " << sd;
			// for (unsigned int j = 0; j < Coefs.size(); j++) {
			// 		cout << " " << Coefs[j].second << "*X_" << Coefs[j].first.first << "_" << Coefs[j].first.second;		
			// }
			// cout << endl;
			s = Coefs[i].first.first;
			t = Coefs[i].first.second;		
			// _formul_supercut->_model->reset(0);
//			cout << "size : " << Coefs.size() << endl;

			gain = _formul_supercut->get_gain(sd, &Coefs, s, t);	
			if(gain < 0){
				Coefs[i].second = max(0.0,  Coefs[i].second + gain);
			}else{
				// cout << "reduce coef " << i << "-th" << endl;
				break;
			}	
		}
		_idx_supercut++;
		// cout << " done with super cut " << _idx_supercut - 1 << endl;			
	}
	
	GRBLinExpr expr2 = 0;
	// cout << sd;
	for (unsigned int i = 0; i < Coefs.size(); i++) {
		expr2 += Coefs[i].second * _var_y[Coefs[i].first.first][Coefs[i].first.second];
		// cout << " - " << Coefs[i].second << "*" << "x_" << Coefs[i].first.first <<','<< Coefs[i].first.second;
	}
	// cout << endl;
	expr2 = expr2 + (*_var_v) - sd;
	return expr2;
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


BendersCuts::~BendersCuts(){
	;
}
