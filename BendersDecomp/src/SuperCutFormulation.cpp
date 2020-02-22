#include "SuperCutFormulation.h"
#include "SubtourCuts.h"

void SuperCutFormulation::create_variables(GRBModel * model_supercut, int N) {
	this->_size_var_x = N;
	this->_model = model_supercut;
	
	this->_model->getEnv().set(GRB_IntParam_LazyConstraints, 1);


	_var_z = new GRBVar();
	*_var_z = _model->addVar(-1000.0, 1000.0, 0.0, GRB_CONTINUOUS, "z");

	_var_a = new GRBVar();
	*_var_a = _model->addVar(0.0, 1.0, 0.0, GRB_BINARY, "a");

	_var_x = new GRBVar*[_size_var_x];
	for (int i = 0; i < _size_var_x; i++) {
		_var_x[i] = new GRBVar[_size_var_x] ;
	}
	for (int i = 0; i < _size_var_x; i++) {
		for (int j = 0; j < _size_var_x; j++) {
			_var_x[i][j] = _model->addVar(0.0, 1.0, 0.0, GRB_BINARY, "x_" + itos(i) + "," + itos(j));
		}
	}

	GRBLinExpr expr = 0;
	expr += 1.0 * (*_var_z);
	_model->setObjective(expr, GRB_MAXIMIZE);
	_model->update();

}

void SuperCutFormulation::set_objective(double alpha_arrival, vector<pair<pair<int, int>, double>> & Coefs) {
	GRBLinExpr expr = 0;
	expr += alpha_arrival;
	for (unsigned int i =0; i < Coefs.size(); i++) {
		expr -= Coefs[i].second * _var_x[Coefs[i].first.first][Coefs[i].first.second];
	}
	 _model->setObjective(expr, GRB_MAXIMIZE);
}


void SuperCutFormulation::set_objective() {
	GRBLinExpr expr = 0;
	expr += 1.0 * (*_var_z);
	 _model->setObjective(expr, GRB_MAXIMIZE);
	 _model->update();
}

void SuperCutFormulation::update_coefs(double alpha_arrival, vector<pair<pair<int, int>, double>> & Coefs) {
	GRBLinExpr expr = 0;
	_var_a->set(GRB_DoubleAttr_LB, alpha_arrival);
	_var_a->set(GRB_DoubleAttr_UB, alpha_arrival+0.000001);

	for (int i = 0; i < _size_var_x; i++) {
		for (int j = 0; j < _size_var_x; j++) {
			_model->chgCoeff(_model->getConstrByName("CGain"), _var_x[i][j], 0.0);
		}
	}	
	_model->update();
	for (unsigned int i =0; i < Coefs.size(); i++) {
		_model->chgCoeff(_model->getConstrByName("CGain"), _var_x[Coefs[i].first.first][Coefs[i].first.second], -1.0*Coefs[i].second);
	}
	 _model->update();
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

	expr1 = 0;
	expr1 += 1.0 * (*_var_a);
	for (int i = 0; i < _size_var_x; i++) {
		for (int j = 0; j < _size_var_x; j++) {
			expr1 -= 1.0 * _var_x[i][j];
		}
	}
	_model->addConstr(expr1 == 1.0 * (*_var_z), "CGain" );
	_model->update();
}


double SuperCutFormulation::solve() {
	try {
		SubtourCuts * cb = new SubtourCuts(_var_x, _size_var_x);
		_model->setCallback(cb);
		_model->optimize();
		// _model->write("./ret/supercutformul_" + to_string(++_idx)+".lp");
		cout << "outer flag ***"<<endl;

		// exit(0);
		if (_model->get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
			cout << "inner flag ***"<<endl;
			double obj_val = _model->get(GRB_DoubleAttr_ObjVal);
			cout << obj_val << endl;
			// cout << "inner flag 2 ---" << endl;	
			// double **sol2 = new double*[_size_var_x];
			// for (int ii = 0; ii < _size_var_x; ii++) {
			// 	sol2[ii] = new double[_size_var_x];
			// 	cout << _size_var_x << endl;
			// 	// sol2[ii] = getSolution(_var_x[ii], _size_var_x);
			// }
			// for (int ki = 0; ki < _size_var_x; ki++) {
			// 	for (int kj = 0; kj < _size_var_x; kj++) {
			// 		cout << _var_x[ki][kj].get(GRB_DoubleAttr_X) << "  ";
			// 	}
			// 	cout << endl;
			// }
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

SuperCutFormulation::~SuperCutFormulation(){
	for (int i = 0; i < _size_var_x; i++)
		delete[] _var_x[i];
	delete[] _var_x;
}
