#include "SubtourCuts.h"

SubtourCuts::SubtourCuts(GRBVar*** xVars, int size_xVars) {
	this->_var_x = xVars;
	this->_size_var_x = size_xVars;
}

void SubtourCuts::callback() {
	try {
			// cout << "where = " << where <<  endl;
		if (where == GRB_CB_MIPSOL) {
			double **x_sol = new double*[_size_var_x];
			for (int i = 0; i < _size_var_x; i++) {
				x_sol[i] = new double[_size_var_x];			
				x_sol[i] = getSolution((*_var_x)[i],_size_var_x);
			}
			int *tour = new int[_size_var_x];
			int len;
			BendersCuts::findsubtour(_size_var_x, x_sol, &len, tour);
			if (len < _size_var_x) {
				// Add subtour elimination constraint
				GRBLinExpr expr = 0;
				for (int i = 0; i < len; i++) {
					for (int j = i + 1; j < len; j++) {
						expr += (*_var_x)[tour[i]][tour[j]] + (*_var_x)[tour[j]][tour[i]];					
					}
				}
				addLazy(expr <= len - 1);									
			}
			for (int i = 0; i < _size_var_x; i++)
				delete[] x_sol[i];
			delete[] x_sol;
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
