
#include "DualFormulation.h"
#define INF numeric_limits<double>::infinity()


void DualFormulation::create_variables(GRBModel * model_dual, PartitionScheme * partition_) {

	this->_model_dual = model_dual;
	this->_G = &(partition_->_G);
	this->_nb_dstzn = partition_->_nb_dstzn;
	this->_nb_targets = partition_->_dataset->_nb_targets;

	/* add dual variables [_alpha] associated with the flow balance constraints */
	_size_alpha = 2 * _nb_targets * _nb_dstzn + 2;
	_alpha = new GRBVar[_size_alpha];
	for (int i = 0; i < _size_alpha; i++) {
		_alpha[i] = _model_dual->addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "alpha_" + itos(i));
	}
	_model_dual->update();

	_beta = new GRBVar*[_size_alpha];
	for (int i = 0; i < _size_alpha; i++) {
		_beta[i] = new GRBVar[_size_alpha];
	}
	int idxmat_1, idxmat_2;
	for (int s = 1; s <= _nb_targets; s++) {
		idxmat_1 = (s - 1) * 2 * _nb_dstzn + 1;
		for (int i = 0; i < _nb_dstzn; i++) {
			if((*_G)[0][idxmat_1+i].first == true){
				_beta[0][idxmat_1 + i] = _model_dual->addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "beta_" + itos(0) + "," + itos(idxmat_1 + i));
			}
		}
		_model_dual->update();
	}
	for (int s = 1; s <= _nb_targets; s++) {
		idxmat_1 = (s - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
		for (int t = 1; t <= _nb_targets; t++) {
			idxmat_2 = (t - 1) * 2 * _nb_dstzn + 1;
			if (s != t) {
				for (int i = 0; i < _nb_dstzn; i++) {
					for (int j = 0; j < _nb_dstzn; j++) {
						if((*_G)[idxmat_1 + i][idxmat_2 + j].first == true){
							_beta[idxmat_1 + i][idxmat_2 + j] = _model_dual->addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "beta_" + itos(idxmat_1 + i) + "," + itos(idxmat_2 + j));
						}
					}
					_model_dual->update();
				}
			}
		}
	}
	_model_dual->update();
	for (int t = 1; t <= _nb_targets; t++) {
		idxmat_1 = (t - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
		for (int i = 0; i < _nb_dstzn; i++) {
			if((*_G)[idxmat_1 + i][_size_alpha - 1].first == true){
				_beta[idxmat_1 + i][_size_alpha - 1] = _model_dual->addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "beta_" + itos(idxmat_1 + i) + "," + itos(_size_alpha - 1));
			}
		}
		_model_dual->update();
	}

	for (int s = 1; s <= _nb_targets; s++) {
		idxmat_1 = (s - 1) * 2 * _nb_dstzn + 1;
		idxmat_2 = (s - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
		for (int i = 0; i < _nb_dstzn; i++) {
			for (int j = 0; j < _nb_dstzn; j++) {
				if ((*_G)[idxmat_1 + i][idxmat_2 + j].first == true){
					_beta[idxmat_1 + i][idxmat_2 + j] = _model_dual->addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "beta_" + itos(idxmat_1 + i) + "," + itos(idxmat_2 + j));
				}
			}
			_model_dual->update();
		}
	}

	// _size_constrs = 2 * _nb_dstzn + (_nb_targets - 1) * _nb_dstzn * _nb_dstzn + _nb_targets * (_nb_dstzn * _nb_dstzn - _nb_dstzn);
}

void DualFormulation::set_objective(double** val_y) {
	double small_value = 0.00001;
	GRBLinExpr expr_obj = 0;
	// // i) dual variable alpha_u''
	expr_obj += 1.0 * _alpha[_size_alpha - 1];
	int idxmat_1, idxmat_2;
	// ii) E_1
	for (int s = 1; s <= _nb_targets; s++) {
		idxmat_1 = (s - 1) * 2 * _nb_dstzn + 1;
		if (val_y[0][s] > small_value) {
			for (int i = 0; i < _nb_dstzn; i++) {
				if((*_G)[0][idxmat_1 + i].first == true){
					expr_obj -= val_y[0][s] * _beta[0][idxmat_1 + i];
				}
			}
		}
	}
	// iii) e in E_3
	for (int s = 1; s <= _nb_targets; s++) {
		idxmat_1 = (s - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
		for (int t = 1; t <= _nb_targets; t++) {
			idxmat_2 = (t - 1) * 2 * _nb_dstzn + 1;
			if (s != t && val_y[s][t] > small_value) {
				for (int i = 0; i < _nb_dstzn; i++) {
					for (int j = 0; j < _nb_dstzn; j++) {
						if((*_G)[idxmat_1 + i][idxmat_2 + j].first == true){
							expr_obj -= val_y[s][t] * _beta[idxmat_1 + i][idxmat_2 + j];
						}
					}
				}
			}
		}
	}

	// iv) e in E_4
	for (int t = 1; t <= _nb_targets; t++) {
		idxmat_1 = (t - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
		if (val_y[t][_nb_targets + 1] > small_value) {
			for (int i = 0; i < _nb_dstzn; i++) {
				if((*_G)[idxmat_1 + i][_size_alpha - 1].first == true){
					expr_obj -= val_y[t][_nb_targets + 1] * _beta[idxmat_1 + i][_size_alpha - 1];
				}
			}
		}
	}

	_model_dual->setObjective(expr_obj, GRB_MAXIMIZE);
	_model_dual->update();
}

void DualFormulation::remove_constraints() {
	for (int i = 0; i < _size_constrs; i++) {
		_model_dual->remove(_model_dual->getConstrs()[i]);
	}
	_size_constrs = 0;
}

void DualFormulation::set_constraints( ) {
	GRBLinExpr expr_constr = 0.0;
	int idxmat_1, idxmat_2;
	for (int s = 1; s <= _nb_targets; s++) {
		idxmat_1 = (s - 1) * 2 * _nb_dstzn + 1;
		for (int i = 0; i < _nb_dstzn; i++) {
			if((*_G)[0][idxmat_1 + i].first == true){
				expr_constr = _alpha[idxmat_1 + i] - _alpha[0] - _beta[0][idxmat_1 + i] - (*_G)[0][idxmat_1 + i].second;
				_model_dual->addConstr(expr_constr <= 0, "C_" + itos(0) + "," + itos(idxmat_1 + i));
				_size_constrs++;
			}
		}
		_model_dual->update();
	}

	for (int s = 1; s <= _nb_targets; s++) {
		idxmat_1 = (s - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
		for (int t = 1; t <= _nb_targets; t++) {
			idxmat_2 = (t - 1) * 2 * _nb_dstzn + 1;
			if (s != t) {
				for (int i = 0; i < _nb_dstzn; i++) {
					for (int j = 0; j < _nb_dstzn; j++) {
						if((*_G)[idxmat_1 + i][idxmat_2 + j].first == true){
							expr_constr = _alpha[idxmat_2 + j] - _alpha[idxmat_1 + i] - _beta[idxmat_1 + i][idxmat_2 + j] - (*_G)[idxmat_1 + i][idxmat_2 + j].second;
							_model_dual->addConstr(expr_constr <= 0, "C_" + itos(idxmat_1 + i) + "," + itos(idxmat_2 + j));
							_size_constrs++;
						}	
					}
				}
				_model_dual->update();
			}
		}
	}
	idxmat_2 = _size_alpha - 1;
	for (int t = 1; t <= _nb_targets; t++) {
		idxmat_1 = (t - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
		for (int i = 0; i < _nb_dstzn; i++) {
			if((*_G)[idxmat_1 + i][idxmat_2].first == true){
				expr_constr = _alpha[idxmat_2] - _alpha[idxmat_1 + i] - _beta[idxmat_1 + i][idxmat_2] - (*_G)[idxmat_1 + i][idxmat_2].second;
				_model_dual->addConstr(expr_constr <= 0, "C_" + itos(idxmat_1 + i) + "," + itos(idxmat_2));
				_size_constrs++;
			}
		}	
		_model_dual->update();
	}

	for (int s = 1; s <= _nb_targets; s++) {
		idxmat_1 = (s - 1) * 2 * _nb_dstzn + 1;
		idxmat_2 = (s - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
		for (int i = 0; i < _nb_dstzn; i++) {
			for (int j = 0; j < _nb_dstzn; j++) {
				if((*_G)[idxmat_1 + i][idxmat_2 + j].first == true){
					expr_constr = _alpha[idxmat_2 + j] - _alpha[idxmat_1 + i] - (*_G)[idxmat_1 + i][idxmat_2 + j].second;
					_model_dual->addConstr(expr_constr <= 0, "C_" + itos(idxmat_1 + i) + "," + itos(idxmat_2 + j));
					_size_constrs++;
				}
			}
		}
	}
	_model_dual->addConstr(_alpha[0] == 0.0);
	_model_dual->update();
}

double DualFormulation::solve() {
	try {
		_model_dual->optimize();
		if (_model_dual->get(GRB_IntAttr_Status) == GRB_OPTIMAL){
			_status = 0;
			double obj_val = _model_dual->get(GRB_DoubleAttr_ObjVal);
			return	obj_val;
		}
	}
	catch (GRBException e) {
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during optimization" << endl;
	}
	return INF;
}

void DualFormulation::get_Benders_user_cut(GRBLinExpr & expr, GRBVar** var_y) {

	int idxmat_1, idxmat_2;
	double coef = 0.0;
	expr += _alpha[_size_alpha - 1].get(GRB_DoubleAttr_X);;
	for (int s = 1; s <= _nb_targets; s++) {
		idxmat_1 = (s - 1) * 2 * _nb_dstzn + 1;
		coef = 0.0;
		for (int j = 0; j < _nb_dstzn; j++) {
			if((*_G)[0][idxmat_1 + j].first == true)
				coef += _beta[0][idxmat_1 + j].get(GRB_DoubleAttr_X);
		}
		expr -= coef * var_y[0][s];
	}
	for (int s = 1; s <= _nb_targets; s++) {
		idxmat_1 = (s - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
		for (int t = 1; t <= _nb_targets; t++) {
			if (s != t) {
				idxmat_2 = (t - 1) * 2 * _nb_dstzn + 1;
				coef = 0.0;
				for (int i = 0; i < _nb_dstzn; i++) {
					for (int j = 0; j < _nb_dstzn; j++) {
						if((*_G)[idxmat_1 + i][idxmat_2 + j].first == true)
							coef += _beta[idxmat_1 + i][idxmat_2 + j].get(GRB_DoubleAttr_X);
					}
				}
				expr -= coef * var_y[s][t];
			}
		}
	}
	for (int t = 1; t <= _nb_targets; t++) {
		idxmat_2 = (t - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
		coef = 0.0;
		for (int j = 0; j < _nb_dstzn; j++) {
			if((*_G)[idxmat_2 + j][_size_alpha - 1].first == true)	
				coef += _beta[idxmat_2 + j][_size_alpha - 1].get(GRB_DoubleAttr_X);
		}
		expr -= coef * var_y[t][_nb_targets + 1];
	}
}

DualFormulation::~DualFormulation() {
	delete _alpha;
	for (int i = 0; i < _size_alpha; i++) {
		delete[] _beta[i];
	}
	delete _beta;
}
