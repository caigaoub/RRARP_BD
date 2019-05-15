
#include "DualFormulation.h"
#define INF numeric_limits<double>::infinity()


DualFormulation::DualFormulation(GRBModel * m_dual, PartitionScheme * ps, int k) {

	this->model_dual = m_dual;
	this->G = ps->get_G();
	this->num_dstzn = k;
	this->num_targets = ps->get_num_targets();


	size_beta = 2 * num_targets * num_dstzn + 2; // number of nodes in  graph for Shortest-Path problem
//	size_alpha = num_targets * (num_targets - 1) * k * k + num_targets * k * k + 2 * k*num_targets;

	// add dual variables [beta] associated with the flow balance constraints
	beta = new GRBVar[size_beta];
//	beta = model_dual->addVars(size_beta);
//	model_dual->update();
	for (int i = 0; i < size_beta; i++) {
		beta[i] = model_dual->addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "beta_" + itos(i));
	}
	model_dual->update();

	alpha = new GRBVar*[size_beta];
	for (int i = 0; i < size_beta; i++) {
		alpha[i] = new GRBVar[size_beta];
	}

	int idxmat_1, idxmat_2;
	for (int s = 1; s <= num_targets; s++) {
		idxmat_1 = (s - 1) * 2 * num_dstzn + 1;
		for (int i = 0; i < num_dstzn; i++) {
			alpha[0][idxmat_1 + i] = model_dual->addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "alpha_" + itos(0) + "_" + itos(idxmat_1 + i));
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
						alpha[idxmat_1 + i][idxmat_2 + j] = model_dual->addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "alpha_" + itos(idxmat_1 + i) + "_" + itos(idxmat_2 + j));
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
			alpha[idxmat_1 + i][size_beta - 1] = model_dual->addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "alpha_" + itos(idxmat_1 + i) + "_" + itos(size_beta - 1));
		}
		model_dual->update();
	}


	for (int s = 1; s <= num_targets; s++) {
		idxmat_1 = (s - 1) * 2 * num_dstzn + 1;
		idxmat_2 = (s - 1) * 2 * num_dstzn + num_dstzn + 1;
		for (int i = 0; i < num_dstzn; i++) {
			for (int j = 0; j < num_dstzn; j++) {
				if (G[idxmat_1 + i][idxmat_2 + j] < INF)
					alpha[idxmat_1 + i][idxmat_2 + j] = model_dual->addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "alpha_" + itos(idxmat_1 + i) + "_" + itos(idxmat_2 + j));
			}
			model_dual->update();
		}
	}

}

void DualFormulation::set_objective(double** y) {
	double small_value = 0.000000001;
	GRBLinExpr expr_obj = 0;

	// add the last dual variable beta_p0^\prime
	expr_obj += 1.0 * beta[size_beta - 1];

	int idxmat_1, idxmat_2;
	for (int s = 1; s <= num_targets; s++) {
		idxmat_1 = (s - 1) * 2 * num_dstzn + 1;
		if (y[0][s] > small_value) {
			for (int i = 0; i < num_dstzn; i++) {
				expr_obj -= y[0][s] * alpha[0][idxmat_1 + i];
			}
		}
	}

	for (int s = 1; s <= num_targets; s++) {
		idxmat_1 = (s - 1) * 2 * num_dstzn + num_dstzn + 1;
		for (int t = 1; t <= num_targets; t++) {
			idxmat_2 = (t - 1) * 2 * num_dstzn + 1;
			if (s != t && y[s][t] > small_value) {
				for (int i = 0; i < num_dstzn; i++) {
					for (int j = 0; j < num_dstzn; j++) {
						expr_obj -= y[s][t] * alpha[idxmat_1 + i][idxmat_2 + j];
					}
				}
			}
		}
	}

	for (int t = 1; t <= num_targets; t++) {
		idxmat_1 = (t - 1) * 2 * num_dstzn + num_dstzn + 1;
		if (y[t][num_targets + 1] > small_value) {
			for (int i = 0; i < num_dstzn; i++) {
				expr_obj -= y[t][num_targets + 1] * alpha[idxmat_1 + i][size_beta - 1];
			}
		}
	}

	for (int s = 1; s <= num_targets; s++) {
		idxmat_1 = (s - 1) * 2 * num_dstzn + 1;
		idxmat_2 = idxmat_1 + num_dstzn;
		for (int i = 0; i < num_dstzn; i++) {
			for (int j = 0; j < num_dstzn; j++) {
				if(G[idxmat_1 + i][idxmat_2 + j] < INF)
					expr_obj -= 1.0 * alpha[idxmat_1 + i][idxmat_2 + j];
			}
		}
	}

	model_dual->setObjective(expr_obj, GRB_MAXIMIZE);
//	model_dual->write("dual.lp");
	model_dual->update();

}
void DualFormulation::remove_constraints() {
	int size_constrs = 2 * num_dstzn + (num_targets - 1) * num_dstzn * num_dstzn + num_targets *  (num_dstzn * num_dstzn - num_dstzn);
	for (int i = 0; i < size_constrs; i++) {
		model_dual->remove(model_dual->getConstrs()[i]);
	}

}

void DualFormulation::set_constraints( ) {

//	double small_value = 0.0000001;

	GRBLinExpr expr_constr = 0.0;
	/*--------------------------------------------------------------------------*/
	int idxmat_1, idxmat_2;
	for (int s = 1; s <= num_targets; s++) {
		idxmat_1 = (s - 1) * 2 * num_dstzn + 1;
		for (int i = 0; i < num_dstzn; i++) {
			expr_constr = beta[idxmat_1 + i] - beta[0] - alpha[0][idxmat_1 + i] - G[0][idxmat_1 + i];
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
						expr_constr = beta[idxmat_2 + j] - beta[idxmat_1 + i] - alpha[idxmat_1 + i][idxmat_2 + j] - G[idxmat_1 + i][idxmat_2 + j];
						model_dual->addConstr(expr_constr <= 0, "C_" + itos(idxmat_1 + i) + "_" + itos(idxmat_2 + j));
					}
				}
				model_dual->update();
			}

		}
	}

	idxmat_2 = size_beta - 1;
	for (int t = 1; t <= num_targets; t++) {
		idxmat_1 = (t - 1) * 2 * num_dstzn + num_dstzn + 1;
		for (int i = 0; i < num_dstzn; i++) {
			expr_constr = beta[idxmat_2] - beta[idxmat_1 + i] - alpha[idxmat_1 + i][idxmat_2] - G[idxmat_1 + i][0];
			model_dual->addConstr(expr_constr <= 0, "C_" + itos(idxmat_1 + i) + "_" + itos(idxmat_2));
		}
		model_dual->update();
	}

	/*--------------------------------------------------------------------------*/

	for (int s = 1; s <= num_targets; s++) {
		idxmat_1 = (s - 1) * 2 * num_dstzn + 1;
		idxmat_2 = (s - 1) * 2 * num_dstzn + num_dstzn + 1;
		for (int i = 0; i < num_dstzn; i++) {
			for (int j = 0; j < num_dstzn; j++) {
				if (G[idxmat_1 + i][idxmat_2 + j] < INF) {
					expr_constr = beta[idxmat_2 + j] - beta[idxmat_1 + i] - alpha[idxmat_1 + i][idxmat_2 + j] - G[idxmat_1 + i][idxmat_2 + j];
					model_dual->addConstr(expr_constr <= 0, "C_" + itos(idxmat_1 + i) + "_" + itos(idxmat_2 + j));
				}
			}
			model_dual->update();
		}

	}
	model_dual->update();
	model_dual->addConstr(beta[0] == 0.0);
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

	double beta_sink = beta[size_beta - 1].get(GRB_DoubleAttr_X);
	expr += beta_sink;

	int idxmat_1, idxmat_2;
	double coef;
	/*--------------------------------------------------------------------------*/
	for (int s = 1; s <= num_targets; s++) {
		idxmat_1 = (s - 1) * 2 * num_dstzn + 1;
		coef = 0.0;
		for (int j = 0; j < num_dstzn; j++) {
			coef += alpha[0][idxmat_1 + j].get(GRB_DoubleAttr_X);
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
						coef += alpha[idxmat_1 + i][idxmat_2 + j].get(GRB_DoubleAttr_X);
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
			coef += alpha[idxmat_2 + j][size_beta - 1].get(GRB_DoubleAttr_X);
		}
		expr -= coef * y[t][num_targets + 1];
	}
	/*--------------------------------------------------------------------------*/
	coef = 0.0;
	for (int s = 1; s <= num_targets; s++) {
		idxmat_1 = (s - 1) * 2 * num_dstzn + 1;
		idxmat_2 = (s - 1) * 2 * num_dstzn + num_dstzn + 1;
		for (int i = 0; i < num_dstzn; i++) {
			for (int j = 0; j < num_dstzn; j++) {
				if (G[idxmat_1 + i][idxmat_2 + j] < INF)
					coef +=  alpha[idxmat_1 + i][idxmat_2 + j].get(GRB_DoubleAttr_X);
			}
		}
	}
	expr -= coef;

}

DualFormulation::~DualFormulation() {

	delete beta;
	for (int i = 0; i < size_beta; i++) {
		delete[] alpha[i];
	}
	delete alpha;

}
