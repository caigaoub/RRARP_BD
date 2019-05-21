#include <stdlib.h>
#include "gurobi_c++.h"
#include "DataHandler.h"
#include "myNameClass.h"
#include "PartitionScheme.h"
#include "STEFormulation.h"
// #include "CoefReduction.h"
#include "BendersCuts.h"
#include <ctime>
#include <chrono>

using namespace std;
void Fischetti_method(int , STEFormulation & );
int main(int argc, const char* argv[]) {

	argc = argc; // just for avoid warning: unused argc
	const char* filename = argv[1];
    const int num_dstzn = atoi(argv[2]);
	try {
//		int num_dstzn = 4;
//		const char* filename = "RRARP_n_7_E_3.txt";
		auto start = chrono::system_clock::now();
		DataHandler instance(filename);
		PartitionScheme ps(num_dstzn, instance);

		GRBEnv * evn_dual = new GRBEnv();
		GRBModel model_dual = GRBModel(*evn_dual);
		model_dual.getEnv().set(GRB_IntParam_OutputFlag, 0);
		DualFormulation DualForm(&model_dual, &ps, num_dstzn);
		DualForm.set_constraints();

		GRBEnv * evn_MP = new GRBEnv();
		GRBModel model_MP = GRBModel(*evn_MP);
		model_MP.getEnv().set(GRB_IntParam_OutputFlag, 0);
		STEFormulation STEForm(&model_MP, &ps, &DualForm);

		int algorithm = 1;
		if (algorithm == 1) {
	     pair<double, double> ret	=	STEForm.solve_IP_TSP();
	//	 STEForm.printSol(&model_MP);
	//   STEForm.print_num_Benders_cuts();
	//	 STEForm.print_num_subtour_cuts();
			auto end = chrono::system_clock::now();
			chrono::duration<double> elapsed_seconds = end-start;
			time_t end_time = chrono::system_clock::to_time_t(end);
		  fstream fs;
	    fs.open("./ret/table2.dat", fstream::app | fstream::out);
	    fs << elapsed_seconds.count() << '\t' << model_MP.get(GRB_DoubleAttr_MIPGap) << '\t' \
			   << STEForm.get_num_subtour_cuts() << '\t' << STEForm.get_num_Benders_cuts() << '\n';
	    fs.close();

		}
		if (algorithm == 2) {
			Fischetti_method(ps.get_num_targets() + 2, STEForm);
			auto end = chrono::system_clock::now();
			chrono::duration<double> elapsed_seconds = end-start;
			time_t end_time = std::chrono::system_clock::to_time_t(end);

			fstream fs;
			fs.open("./ret/table2.dat", fstream::app | fstream::out);
			fs << elapsed_seconds.count() << '\t' << model_MP.get(GRB_DoubleAttr_MIPGap) << '\t' \
			   << STEForm.get_num_subtour_cuts() << '\t' << STEForm.get_num_Benders_cuts() << '\n';
	    fs.close();
		}
	}
	catch (const GRBException& ex) {
		cout << "Error number: " << ex.getErrorCode() << endl;
		cout << ex.getMessage() << endl;
	}
	catch (...) {
		cerr << "Error" << endl;
	}
//	system("pause");
	return 0;
}


void Fischetti_method(int N, STEFormulation & stef) {

	stef.set_model_LP(); // solve LP-TSP model

	// create variables
	double** y_tilde = new double*[N]; // stablizer point
	double** y_star = new double*[N]; // optimal solution of current TSP_LP
	double** y_bar = new double*[N];
	for (int i = 0; i < N; i++) {
		y_tilde[i] = new double[N];
		y_star[i] = new double[N];
		y_bar[i] = new double[N];
	}

	//initialize y tilde to represent the feasible sequence (0 -> 1 -> 2 ->...->N-1)
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			y_tilde[i][j] = 0;
		}
	}
	for (int i = 0; i < N - 1; i++) {
		y_tilde[i][i + 1] = 1;
	}
	y_tilde[N - 1][0] = 1;


//	double val_MP,
//	double objVal, v_Val;
	stef.solve_LP_TSP();
	stef.get_optimal_sol(y_star);

	/*
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			cout << y_star[i][j] << '\t';
		}
		cout << '\n';
	}
	cout << "--------------------------------------------" << endl;

	*/
	double alpha = 0.2;
	double UB = INFINITY;
	double LB = -INFINITY;

	double dual_obj;
	double old_LB = -1;
	int iteration = 0;
	bool is_disconnected = true;


	while (is_disconnected || UB - LB > 0.0000001) {
		cout << "UB: " << UB << "   " << "LB: " << LB << '\n';

		is_disconnected = stef.add_SECs(y_star); //if one subtour elimi constrait is added, flag is
		if (!is_disconnected) {
			for (int i = 0; i < N; i++) { // update the stablizer
				for (int j = 0; j < N; j++) {
					y_tilde[i][j] = 0.5 * (y_star[i][j] + y_tilde[i][j]);
				}
			}
			for (int i = 0; i < N; i++) { // new solution for the dual model
				for (int j = 0; j < N; j++) {
					y_bar[i][j] = alpha * y_star[i][j] + (1 - alpha) * y_tilde[i][j];
				}
			}
			dual_obj = stef.add_USER_cuts(y_bar);
		}

		pair<double, double> result_tsp =stef.solve_LP_TSP();
		stef.get_optimal_sol(y_star);

		if (!is_disconnected) {
			LB = result_tsp.first;
			if (old_LB == LB) {
				iteration++;
				if (iteration >= 5) {
					alpha = 1;
				}
			}
			old_LB = LB;
			UB = min(UB, result_tsp.first - result_tsp.second + dual_obj);

		}
	}
	cout << "UB: " << UB << "   " << "LB: " << LB << '\n';

	cout << "Gap is closed!" << '\n';


	stef.set_model_MIP(); // solve IP-TSP
	stef.solve_IP_TSP(); // add Benders cuts

}
