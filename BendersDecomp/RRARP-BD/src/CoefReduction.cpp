#include "CoefReduction.h"
#include "SubtourCuts.h"

CoefReduction::CoefReduction(GRBModel * tsp2, int N) {
	this->N = N;
	this->model = tsp2;
	model->getEnv().set(GRB_IntParam_LazyConstraints, 1);
	w = vector<vector<GRBVar>>(N);
	for (int i = 0; i < N; i++) {
		w[i].resize(N);
	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			w[i][j] = tsp2->addVar(0.0, 1.0, 0.0, GRB_BINARY, "w_" + itos(i) + "_" + itos(j));
		}
	}
}

void CoefReduction::set_objective(double alpha_endDepot, vector<tuple<int, int, double>> & CoefSet) {
	GRBLinExpr expr = 0;
	expr += alpha_endDepot;
	for (unsigned int i = 0; i < CoefSet.size(); i++) {
		expr -= get<2>(CoefSet[i]) * w[get<0>(CoefSet[i])][get<1>(CoefSet[i])];
//		cout << get<0>(CoefSet[i]) << " " << get<1>(CoefSet[i]) << endl;
	}
	model->setObjective(expr, GRB_MAXIMIZE);
}


void CoefReduction::set_constraints() {
	GRBLinExpr expr1, expr2;
	for (int i = 0; i < N; i++) {
		expr1 = 0;
		expr2 = 0;
		for (int j = 0; j < N; j++) {
			expr1 += w[i][j];
			expr2 += w[j][i];
		}
		model->addConstr(expr1 == 1, "Deg1_Row" + itos(i));
		model->addConstr(expr2 == 1, "Deg1_Col" + itos(i));
		model->update();
	}

	w[N - 1][0].set(GRB_DoubleAttr_LB, 1);
	model->update();
	for (int i = 0; i < N; i++) {
		w[i][i].set(GRB_DoubleAttr_UB, 0);
	}
	model->update();
}


double CoefReduction::solve() {
	try {
		SubtourCuts * cb = new SubtourCuts(w, N);
		model->setCallback(cb);

		model->optimize();
//		model->write("tsp3.lp");
		if (model->get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
			double obj_val = model->get(GRB_DoubleAttr_ObjVal);
			delete cb;
			return obj_val;
		}else{
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
	return -1;
}

void CoefReduction::fix_edge(int i, int j) {
	w[i][j].set(GRB_DoubleAttr_LB, 1);
	model->update();
}

void CoefReduction::free_edge(int i, int j) {
	w[i][j].set(GRB_DoubleAttr_LB, 0);
	model->update();
}
