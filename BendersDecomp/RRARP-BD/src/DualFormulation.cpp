
// #include "DualFormulation.h"
// #define INF numeric_limits<double>::infinity()


// DualFormulation::build_formul(GRBModel * model_dual, PartitionScheme * partition_) {

// 	this->_model_dual = model_dual;
// 	this->_G = partition_->_G;
// 	this->_nb_dstzn = partition_->_nb_dstzn;
// 	this->_nb_targets = partition_->_dataset->_nb_targets;

// 	/* add dual variables [alpha] associated with the flow balance constraints */
// 	_size_alpha = 2 * _nb_targets * _nb_dstzn + 2;
// 	_alpha = new GRBVar[_size_alpha];
// 	for (int i = 0; i < _size_alpha; i++) {
// 		_alpha[i] = _model_dual->addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "alpha_" + itos(i));
// 	}
// 	_model_dual->update();

// 	_beta = new GRBVar*[_size_alpha];
// 	for (int i = 0; i < _size_alpha; i++) {
// 		_beta[i] = new GRBVar[_size_alpha];
// 	}
// 	int idxmat_1, idxmat_2;
// 	for (int s = 1; s <= _nb_targets; s++) {
// 		idxmat_1 = (s - 1) * 2 * _nb_dstzn + 1;
// 		for (int i = 0; i < _nb_dstzn; i++) {
// 			_beta[0][idxmat_1 + i] = _model_dual->addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "beta_" + itos(0) + "_" + itos(idxmat_1 + i));
// 		}
// 		_model_dual->update();
// 	}
// 	for (int s = 1; s <= _nb_targets; s++) {
// 		idxmat_1 = (s - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
// 		for (int t = 1; t <= _nb_targets; t++) {
// 			idxmat_2 = (t - 1) * 2 * _nb_dstzn + 1;
// 			if (s != t) {
// 				for (int i = 0; i < _nb_dstzn; i++) {
// 					for (int j = 0; j < _nb_dstzn; j++) {
// 						_beta[idxmat_1 + i][idxmat_2 + j] = _model_dual->addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "beta_" + itos(idxmat_1 + i) + "_" + itos(idxmat_2 + j));
// 					}
// 					_model_dual->update();
// 				}
// 			}
// 		}
// 	}
// 	_model_dual->update();
// 	for (int t = 1; t <= _nb_targets; t++) {
// 		idxmat_1 = (t - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
// 		for (int i = 0; i < _nb_dstzn; i++) {
// 			_beta[idxmat_1 + i][_size_alpha - 1] = _model_dual->addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "beta_" + itos(idxmat_1 + i) + "_" + itos(_size_alpha - 1));
// 		}
// 		_model_dual->update();
// 	}
// 	for (int s = 1; s <= _nb_targets; s++) {
// 		idxmat_1 = (s - 1) * 2 * _nb_dstzn + 1;
// 		idxmat_2 = (s - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
// 		for (int i = 0; i < _nb_dstzn; i++) {
// 			for (int j = 0; j < _nb_dstzn; j++) {
// 				if (_G[idxmat_1 + i][idxmat_2 + j] < INF)
// 					_beta[idxmat_1 + i][idxmat_2 + j] = _model_dual->addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "beta_" + itos(idxmat_1 + i) + "_" + itos(idxmat_2 + j));
// 			}
// 			_model_dual->update();
// 		}
// 	}

// 	_size_constrs = 2 * _nb_dstzn + (_nb_targets - 1) * _nb_dstzn * _nb_dstzn + _nb_targets * (_nb_dstzn * _nb_dstzn - _nb_dstzn);
// }

// void DualFormulation::set_objective(double** y) {
// 	double small_value = 0.00001;
// 	GRBLinExpr expr_obj = 0;

// 	// i) dual variable alpha_u''
// 	expr_obj += 1.0 * alpha[_size_alpha - 1];
// 	int idxmat_1, idxmat_2;
// 	// ii) E_1
// 	for (int s = 1; s <= _nb_targets; s++) {
// 		idxmat_1 = (s - 1) * 2 * _nb_dstzn + 1;
// 		if (y[0][s] > small_value) {
// 			for (int i = 0; i < _nb_dstzn; i++) {
// 				expr_obj -= y[0][s] * beta[0][idxmat_1 + i];
// 			}
// 		}
// 	}
// 	// iii) e in E_3
// 	for (int s = 1; s <= _nb_targets; s++) {
// 		idxmat_1 = (s - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
// 		for (int t = 1; t <= _nb_targets; t++) {
// 			idxmat_2 = (t - 1) * 2 * _nb_dstzn + 1;
// 			if (s != t && y[s][t] > small_value) {
// 				for (int i = 0; i < _nb_dstzn; i++) {
// 					for (int j = 0; j < _nb_dstzn; j++) {
// 						expr_obj -= y[s][t] * beta[idxmat_1 + i][idxmat_2 + j];
// 					}
// 				}
// 			}
// 		}
// 	}
// 	// iv) e in E_4
// 	for (int t = 1; t <= _nb_targets; t++) {
// 		idxmat_1 = (t - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
// 		if (y[t][_nb_targets + 1] > small_value) {
// 			for (int i = 0; i < _nb_dstzn; i++) {
// 				expr_obj -= y[t][_nb_targets + 1] * beta[idxmat_1 + i][_size_alpha - 1];
// 			}
// 		}
// 	}
// 	_model_dual->setObjective(expr_obj, GRB_MAXIMIZE);
// 	_model_dual->update();
// }

// void DualFormulation::remove_constraints() {
// 	for (int i = 0; i < _size_constrs; i++) {
// 		_model_dual->remove(_model_dual->getConstrs()[i]);
// 	}
// }

// void DualFormulation::set_constraints( ) {
// 	GRBLinExpr expr_constr = 0.0;
// 	int idxmat_1, idxmat_2;
// 	for (int s = 1; s <= _nb_targets; s++) {
// 		idxmat_1 = (s - 1) * 2 * _nb_dstzn + 1;
// 		for (int i = 0; i < _nb_dstzn; i++) {
// 			expr_constr = alpha[idxmat_1 + i] - alpha[0] - beta[0][idxmat_1 + i] - G[0][idxmat_1 + i];
// 			_model_dual->addConstr(expr_constr <= 0, "C_" + itos(0) + "_" + itos(idxmat_1 + i));
// 		}
// 		_model_dual->update();
// 	}
// 	for (int s = 1; s <= _nb_targets; s++) {
// 		idxmat_1 = (s - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
// 		for (int t = 1; t <= _nb_targets; t++) {
// 			idxmat_2 = (t - 1) * 2 * _nb_dstzn + 1;
// 			if (s != t) {
// 				for (int i = 0; i < _nb_dstzn; i++) {
// 					for (int j = 0; j < _nb_dstzn; j++) {
// 						expr_constr = alpha[idxmat_2 + j] - alpha[idxmat_1 + i] - beta[idxmat_1 + i][idxmat_2 + j] - G[idxmat_1 + i][idxmat_2 + j];
// 						_model_dual->addConstr(expr_constr <= 0, "C_" + itos(idxmat_1 + i) + "_" + itos(idxmat_2 + j));
// 					}
// 				}
// 				_model_dual->update();
// 			}
// 		}
// 	}
// 	idxmat_2 = _size_alpha - 1;
// 	for (int t = 1; t <= _nb_targets; t++) {
// 		idxmat_1 = (t - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
// 		for (int i = 0; i < _nb_dstzn; i++) {
// 			expr_constr = alpha[idxmat_2] - alpha[idxmat_1 + i] - beta[idxmat_1 + i][idxmat_2] - G[idxmat_1 + i][idxmat_2];
// 			_model_dual->addConstr(expr_constr <= 0, "C_" + itos(idxmat_1 + i) + "_" + itos(idxmat_2));
// 		}
// 		_model_dual->update();
// 	}

// 	for (int s = 1; s <= _nb_targets; s++) {
// 		idxmat_1 = (s - 1) * 2 * _nb_dstzn + 1;
// 		idxmat_2 = (s - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
// 		for (int i = 0; i < _nb_dstzn; i++) {
// 			for (int j = 0; j < _nb_dstzn; j++) {
// 				expr_constr = alpha[idxmat_2 + j] - alpha[idxmat_1 + i] - G[idxmat_1 + i][idxmat_2 + j];
// 				_model_dual->addConstr(expr_constr <= 0, "C_" + itos(idxmat_1 + i) + "_" + itos(idxmat_2 + j));
// 			}
// 		}
// 	}
// 	_model_dual->addConstr(alpha[0] == 0.0);
// 	_model_dual->update();
// }

// double DualFormulation::solve() {
// 	try {
// //		model_dual->write("dual.lp");
// 		_model_dual->optimize();
// 		if (_model_dual->get(GRB_IntAttr_Status) == GRB_OPTIMAL){
// 			_status_dual = 0;
// 			double obj_val = _model_dual->get(GRB_DoubleAttr_ObjVal);
// 			return	obj_val;
// 		}
// 	}
// 	catch (GRBException e) {
// 		cout << "Error number: " << e.getErrorCode() << endl;
// 		cout << e.getMessage() << endl;
// 	}
// 	catch (...) {
// 		cout << "Error during optimization" << endl;
// 	}
// 	return INF;
// }

// void DualFormulation::get_Benders_user_cut(GRBLinExpr & expr, GRBVar** y) {

// 	int idxmat_1, idxmat_2;
// 	double coef;
// 	double beta_endDepot = alpha[_size_alpha - 1].get(GRB_DoubleAttr_X);
// 	expr += beta_endDepot;
// 	for (int s = 1; s <= _nb_targets; s++) {
// 		idxmat_1 = (s - 1) * 2 * _nb_dstzn + 1;
// 		coef = 0.0;
// 		for (int j = 0; j < _nb_dstzn; j++) {
// 			coef += beta[0][idxmat_1 + j].get(GRB_DoubleAttr_X);
// 		}
// 		expr -= coef * y[0][s];
// 	}
// 	for (int s = 1; s <= _nb_targets; s++) {
// 		idxmat_1 = (s - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
// 		for (int t = 1; t <= _nb_targets; t++) {
// 			if (s != t) {
// 				idxmat_2 = (t - 1) * 2 * _nb_dstzn + 1;
// 				coef = 0.0;
// 				for (int i = 0; i < _nb_dstzn; i++) {
// 					for (int j = 0; j < _nb_dstzn; j++) {
// 						coef += beta[idxmat_1 + i][idxmat_2 + j].get(GRB_DoubleAttr_X);
// 					}
// 				}
// 				expr -= coef * y[s][t];
// 			}
// 		}
// 	}
// 	for (int t = 1; t <= _nb_targets; t++) {
// 		idxmat_2 = (t - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
// 		coef = 0.0;
// 		for (int j = 0; j < _nb_dstzn; j++) {
// 			coef += beta[idxmat_2 + j][_size_alpha - 1].get(GRB_DoubleAttr_X);
// 		}
// 		expr -= coef * y[t][_nb_targets + 1];
// 	}
// }

// DualFormulation::~DualFormulation() {
// 	delete alpha;
// 	for (int i = 0; i < _size_alpha; i++) {
// 		delete[] beta[i];
// 	}
// 	delete beta;
// }
