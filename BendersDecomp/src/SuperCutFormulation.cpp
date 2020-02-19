#include "SuperCutFormulation.h"
#include "SubtourCuts.h"

void SuperCutFormulation::create_variables(GRBModel * model_supercut, int N) {
	this->_size_var_x = N;
	this->_model = model_supercut;
	
	this->_model->getEnv().set(GRB_IntParam_LazyConstraints, 1);
	_var_x = new GRBVar*[_size_var_x];
	for (int i = 0; i < _size_var_x; i++) {
		_var_x[i] = new GRBVar[_size_var_x] ;
	}
	for (int i = 0; i < _size_var_x; i++) {
		for (int j = 0; j < _size_var_x; j++) {
			_var_x[i][j] = _model->addVar(0.0, 1.0, 0.0, GRB_BINARY, "x_" + itos(i) + "," + itos(j));
		}
	}

}

void SuperCutFormulation::set_objective(double alpha_arrival, vector<tuple<int, int, double>> & CoefSet) {
	GRBLinExpr expr = 0;
	expr += alpha_arrival;
	for (unsigned int i = 0; i < CoefSet.size(); i++) {
		expr -= get<2>(CoefSet[i]) * _var_x[get<0>(CoefSet[i])][get<1>(CoefSet[i])];
//		cout << get<0>(CoefSet[i]) << " " << get<1>(CoefSet[i]) << endl;
	}
	_model->setObjective(expr, GRB_MAXIMIZE);
}


void SuperCutFormulation::set_constraints() {
	GRBLinExpr expr1, expr2;
	for (int i = 0; i < _size_var_x; i++) {
		expr1 = 0;
		expr2 = 0;
		for (int j = 0; j < _size_var_x; j++) {
			expr1 += _var_x[i][j];
			expr2 += _var_x[j][i];
		}
		_model->addConstr(expr1 == 1, "Deg1_Row" + itos(i));
		_model->addConstr(expr2 == 1, "Deg1_Col" + itos(i));
		_model->update();
	}

	_var_x[_size_var_x - 1][0].set(GRB_DoubleAttr_LB, 1);
	_model->update();
	for (int i = 0; i < _size_var_x; i++) {
		_var_x[i][i].set(GRB_DoubleAttr_UB, 0);
	}
	_model->update();
}


double SuperCutFormulation::solve() {
	try {
		// _model->write("./ret/supercutformul.lp");
		SubtourCuts * cb = new SubtourCuts(_var_x, _size_var_x);
		_model->setCallback(cb);
		_model->optimize();
		if (_model->get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
			double obj_val = _model->get(GRB_DoubleAttr_ObjVal);
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

void SuperCutFormulation::fix_edge(int i, int j) {
	_var_x[i][j].set(GRB_DoubleAttr_LB, 1);
	_model->update();
}

void SuperCutFormulation::free_edge(int i, int j) {
	_var_x[i][j].set(GRB_DoubleAttr_LB, 0);
	_model->update();
}
