
#include "DualFormulation.h"
#define INF numeric_limits<double>::infinity()


DualFormulation::DualFormulation(GRBModel * m_dual, PartitionScheme * ps, int k) {

	this->model_dual = m_dual;
	this->G = ps->get_G();
	this->num_dstzn = k;
	this->num_targets = ps->get_num_targets();

	/* add dual variables [alpha] associated with the flow balance constraints */
	size_alpha = 2 * num_targets * num_dstzn + 2;
	alpha = new GRBVar[size_alpha];
	for (int i = 0; i < size_alpha; i++) {
		alpha[i] = model_dual->addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "alpha_" + itos(i));
	}
	model_dual->update();

	beta = new GRBVar*[size_alpha];
	for (int i = 0; i < size_alpha; i++) {
		beta[i] = new GRBVar[size_alpha];
	}
	int idxmat_1, idxmat_2;
	for (int s = 1; s <= num_targets; s++) {
		idxmat_1 = (s - 1) * 2 * num_dstzn + 1;
		for (int i = 0; i < num_dstzn; i++) {
			beta[0][idxmat_1 + i] = model_dual->addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "beta_" + itos(0) + "_" + itos(idxmat_1 + i));
		}
		model_dual->update();
	}
	for (int s = 1; s <= num_targets; s++) {
		idxmat_1 = (s - 1) * 2 * num_dstzn + num_dstzn + 1;
		for (int t = 1; t <= num_targets; t++) {
			idxmat_2 = (t - 1) * 2 * num_dstzn + 1;
			if (s != t) {
				for (int i = 0; i < num_dstzn; i++) {
					for (int j = 0; j < num_dstzn; j++) {
						beta[idxmat_1 + i][idxmat_2 + j] = model_dual->addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "beta_" + itos(idxmat_1 + i) + "_" + itos(idxmat_2 + j));
					}
					model_dual->update();
				}
			}
		}
	}
	model_dual->update();
	for (int t = 1; t <= num_targets; t++) {
		idxmat_1 = (t - 1) * 2 * num_dstzn + num_dstzn + 1;
		for (int i = 0; i < num_dstzn; i++) {
			beta[idxmat_1 + i][size_alpha - 1] = model_dual->addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "beta_" + itos(idxmat_1 + i) + "_" + itos(size_alpha - 1));
		}
		model_dual->update();
	}
	for (int s = 1; s <= num_targets; s++) {
		idxmat_1 = (s - 1) * 2 * num_dstzn + 1;
		idxmat_2 = (s - 1) * 2 * num_dstzn + num_dstzn + 1;
		for (int i = 0; i < num_dstzn; i++) {
			for (int j = 0; j < num_dstzn; j++) {
				if (G[idxmat_1 + i][idxmat_2 + j] < INF)
					beta[idxmat_1 + i][idxmat_2 + j] = model_dual->addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "beta_" + itos(idxmat_1 + i) + "_" + itos(idxmat_2 + j));
			}
			model_dual->update();
		}
	}

	size_constrs = 2 * num_dstzn + (num_targets - 1) * num_dstzn * num_dstzn + num_targets * (num_dstzn * num_dstzn - num_dstzn);
}

void DualFormulation::set_objective(double** y) {
	double small_value = 0.00001;
	GRBLinExpr expr_obj = 0;

	// i) dual variable alpha_u''
	expr_obj += 1.0 * alpha[size_alpha - 1];
	int idxmat_1, idxmat_2;
	// ii) E_1
	for (int s = 1; s <= num_targets; s++) {
		idxmat_1 = (s - 1) * 2 * num_dstzn + 1;
		if (y[0][s] > small_value) {
			for (int i = 0; i < num_dstzn; i++) {
				expr_obj -= y[0][s] * beta[0][idxmat_1 + i];
			}
		}
	}
	// iii) e in E_3
	for (int s = 1; s <= num_targets; s++) {
		idxmat_1 = (s - 1) * 2 * num_dstzn + num_dstzn + 1;
		for (int t = 1; t <= num_targets; t++) {
			idxmat_2 = (t - 1) * 2 * num_dstzn + 1;
			if (s != t && y[s][t] > small_value) {
				for (int i = 0; i < num_dstzn; i++) {
					for (int j = 0; j < num_dstzn; j++) {
						expr_obj -= y[s][t] * beta[idxmat_1 + i][idxmat_2 + j];
					}
				}
			}
		}
	}
	// iv) e in E_4
	for (int t = 1; t <= num_targets; t++) {
		idxmat_1 = (t - 1) * 2 * num_dstzn + num_dstzn + 1;
		if (y[t][num_targets + 1] > small_value) {
			for (int i = 0; i < num_dstzn; i++) {
				expr_obj -= y[t][num_targets + 1] * beta[idxmat_1 + i][size_alpha - 1];
			}
		}
	}
	model_dual->setObjective(expr_obj, GRB_MAXIMIZE);
	model_dual->update();
}

void DualFormulation::remove_constraints() {
	for (int i = 0; i < size_constrs; i++) {
		model_dual->remove(model_dual->getConstrs()[i]);
	}
}

void DualFormulation::set_constraints( ) {
	GRBLinExpr expr_constr = 0.0;
	int idxmat_1, idxmat_2;
	for (int s = 1; s <= num_targets; s++) {
		idxmat_1 = (s - 1) * 2 * num_dstzn + 1;
		for (int i = 0; i < num_dstzn; i++) {
			expr_constr = alpha[idxmat_1 + i] - alpha[0] - beta[0][idxmat_1 + i] - G[0][idxmat_1 + i];
			model_dual->addConstr(expr_constr <= 0, "C_" + itos(0) + "_" + itos(idxmat_1 + i));
		}
		model_dual->update();
	}
	for (int s = 1; s <= num_targets; s++) {
		idxmat_1 = (s - 1) * 2 * num_dstzn + num_dstzn + 1;
		for (int t = 1; t <= num_targets; t++) {
			idxmat_2 = (t - 1) * 2 * num_dstzn + 1;
			if (s != t) {
				for (int i = 0; i < num_dstzn; i++) {
					for (int j = 0; j < num_dstzn; j++) {
						expr_constr = alpha[idxmat_2 + j] - alpha[idxmat_1 + i] - beta[idxmat_1 + i][idxmat_2 + j] - G[idxmat_1 + i][idxmat_2 + j];
						model_dual->addConstr(expr_constr <= 0, "C_" + itos(idxmat_1 + i) + "_" + itos(idxmat_2 + j));
					}
				}
				model_dual->update();
			}
		}
	}
	idxmat_2 = size_alpha - 1;
	for (int t = 1; t <= num_targets; t++) {
		idxmat_1 = (t - 1) * 2 * num_dstzn + num_dstzn + 1;
		for (int i = 0; i < num_dstzn; i++) {
			expr_constr = alpha[idxmat_2] - alpha[idxmat_1 + i] - beta[idxmat_1 + i][idxmat_2] - G[idxmat_1 + i][idxmat_2];
			model_dual->addConstr(expr_constr <= 0, "C_" + itos(idxmat_1 + i) + "_" + itos(idxmat_2));
		}
		model_dual->update();
	}

	for (int s = 1; s <= num_targets; s++) {
		idxmat_1 = (s - 1) * 2 * num_dstzn + 1;
		idxmat_2 = (s - 1) * 2 * num_dstzn + num_dstzn + 1;
		for (int i = 0; i < num_dstzn; i++) {
			for (int j = 0; j < num_dstzn; j++) {
				expr_constr = alpha[idxmat_2 + j] - alpha[idxmat_1 + i] - G[idxmat_1 + i][idxmat_2 + j];
				model_dual->addConstr(expr_constr <= 0, "C_" + itos(idxmat_1 + i) + "_" + itos(idxmat_2 + j));
			}
		}
	}
	model_dual->addConstr(alpha[0] == 0.0);
	model_dual->update();
}

double DualFormulation::solve() {
	try {
//		model_dual->write("dual.lp");
		model_dual->optimize();
		if (model_dual->get(GRB_IntAttr_Status) == GRB_OPTIMAL){
			status_dual = 0;
			double obj_val = model_dual->get(GRB_DoubleAttr_ObjVal);
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

void DualFormulation::get_Benders_user_cut(GRBLinExpr & expr, GRBVar** y) {

	int idxmat_1, idxmat_2;
	double coef;
	double beta_endDepot = alpha[size_alpha - 1].get(GRB_DoubleAttr_X);
	expr += beta_endDepot;
	for (int s = 1; s <= num_targets; s++) {
		idxmat_1 = (s - 1) * 2 * num_dstzn + 1;
		coef = 0.0;
		for (int j = 0; j < num_dstzn; j++) {
			coef += beta[0][idxmat_1 + j].get(GRB_DoubleAttr_X);
		}
		expr -= coef * y[0][s];
	}
	for (int s = 1; s <= num_targets; s++) {
		idxmat_1 = (s - 1) * 2 * num_dstzn + num_dstzn + 1;
		for (int t = 1; t <= num_targets; t++) {
			if (s != t) {
				idxmat_2 = (t - 1) * 2 * num_dstzn + 1;
				coef = 0.0;
				for (int i = 0; i < num_dstzn; i++) {
					for (int j = 0; j < num_dstzn; j++) {
						coef += beta[idxmat_1 + i][idxmat_2 + j].get(GRB_DoubleAttr_X);
					}
				}
				expr -= coef * y[s][t];
			}
		}
	}
	for (int t = 1; t <= num_targets; t++) {
		idxmat_2 = (t - 1) * 2 * num_dstzn + num_dstzn + 1;
		coef = 0.0;
		for (int j = 0; j < num_dstzn; j++) {
			coef += beta[idxmat_2 + j][size_alpha - 1].get(GRB_DoubleAttr_X);
		}
		expr -= coef * y[t][num_targets + 1];
	}
}

DualFormulation::~DualFormulation() {
	delete alpha;
	for (int i = 0; i < size_alpha; i++) {
		delete[] beta[i];
	}
	delete beta;
}
