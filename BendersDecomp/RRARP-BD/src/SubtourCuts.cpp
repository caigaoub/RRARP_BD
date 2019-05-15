#include "SubtourCuts.h"

SubtourCuts::SubtourCuts(vector<vector<GRBVar>> & wVar, int N) {
	this->w = wVar;
	this->N = N;

}

void SubtourCuts::callback() {
	try {
		if (where == GRB_CB_MIPSOL) {
			double **w_sol = new double*[N];
			for (int i = 0; i < N; i++) {
				w_sol[i] = new double[N];
				for (int j = 0; j < N; j++) {				
					w_sol[i][j] = getSolution(w[i][j]);
				}				
			}
			/*
			for (int i = 0; i < N; i++) {
				for (int j = 0; j < N; j++) {
					cout << w_sol[i][j] << "  ";
				}
				cout << endl;
			}
			*/
			int *tour = new int[N];
			int len;
			BendersCuts::findsubtour(N, w_sol, &len, tour);
			if (len < N) {
				// Add subtour elimination constraint
				GRBLinExpr expr = 0;
				for (int i = 0; i < len; i++) {
					for (int j = i + 1; j < len; j++) {
						expr += w[tour[i]][tour[j]] + w[tour[j]][tour[i]];					
					}
				}
				addLazy(expr <= len - 1);									
			}

			for (int i = 0; i < N; i++)
				delete[] w_sol[i];
			delete[] w_sol;
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
