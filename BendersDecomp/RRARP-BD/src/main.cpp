#include <stdlib.h>
#include "gurobi_c++.h"
#include "DataHandler.h"
#include "PartitionScheme.h"
#include "STEFormulation.h"
// #include "CoefReduction.h"
#include "BendersCuts.h"
#include <ctime>
#include <chrono>

using namespace std;
// void Fischetti_method(int , STEFormulation & );
int main(int argc, const char* argv[]) {

	argc = argc; // just for avoid warning: unused argc
	const int nb_dstzn = atoi(argv[1]);
	string filename = argv[2];    
	try {
		auto start = chrono::system_clock::now();
		DataHandler dataset_;
		dataset_.parse(filename);
		// dataset_.print();
		PartitionScheme network_;
		network_.build(dataset_, nb_dstzn);
		
		vector<int> fseq_;
		fseq_.push_back(0);
		fseq_.push_back(1);
		fseq_.push_back(2);
		fseq_.push_back(3);
		fseq_.push_back(4);
		fseq_.push_back(5);
		fseq_.push_back(6);
		vector<vector<double>> * SDS = new vector<vector<double>>(dataset_._nb_targets + 2);
		network_.solve_shortestpath_v2(fseq_, *SDS);


		// GRBEnv * evn_dual_ = new GRBEnv();
		// GRBModel model_dual_ = GRBModel(*evn_dual_);
		// model_dual_.getEnv().set(GRB_IntParam_OutputFlag, 0);
	
		// DualFormulation dualform_(&model_dual, &ps, num_dstzn);
		// dualform_.set_constraints();

		GRBEnv * evn_MP_ = new GRBEnv();
		GRBModel model_MP_ = GRBModel(*evn_MP_);
		// model_MP_.getEnv().set(GRB_DoubleParam_TimeLimit, 7200);
//		model_MP_.getEnv().set(GRB_IntParam_OutputFlag, 0);
		STEFormulation formul_master;
		formul_master.build_formul(&model_MP_, &network_);

		int algorithm = 0;
		if (algorithm == 1) {
			auto end = chrono::system_clock::now();
			formul_master.solve_IP_TSP();
			formul_master.printSol(&model_MP_);
			chrono::duration<double> elapsed_seconds = end-start;
		  	// fstream fs;
	    // 	fs.open("./ret/table2.dat", fstream::app | fstream::out);
	    // 	fs << "Time(secs):" << elapsed_seconds.count() << '\t' << model_MP_.get(GRB_DoubleAttr_MIPGap) << '\t' << formul_master._total_nb_subtour_cuts << '\t' << formul_master._total_nb_Benders_cuts << '\n';
	    // 	fs.close();
		}
		if (algorithm == 2) {
			// Fischetti_method(ps.get_num_targets() + 2, STEForm);
			// auto end = chrono::system_clock::now();
			// chrono::duration<double> elapsed_seconds = end-start;
			// fstream fs;
			// fs.open("./ret/table2.dat", fstream::app | fstream::out);
	    	// fs << elapsed_seconds.count() << '\t' << model_MP_.get(GRB_DoubleAttr_MIPGap) << '\t' << formul_master._total_nb_subtour_cuts << '\t' << formul_master._total_nb_Benders_cuts << '\n';
	      //  		fs.close();
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


// void Fischetti_method(int N, STEFormulation & stef) {

// 	stef.set_model_LP(); // solve LP-TSP model

// 	// create variables
// 	double** y_tilde = new double*[N]; // stablizer point
// 	double** y_star = new double*[N]; // optimal solution of current TSP_LP
// 	double** y_bar = new double*[N];
// 	for (int i = 0; i < N; i++) {
// 		y_tilde[i] = new double[N];
// 		y_star[i] = new double[N];
// 		y_bar[i] = new double[N];
// 	}

// 	//initialize y tilde to represent the feasible sequence (0 -> 1 -> 2 ->...->N-1)
// 	for (int i = 0; i < N; i++) {
// 		for (int j = 0; j < N; j++) {
// 			y_tilde[i][j] = 0;
// 		}
// 	}
// 	for (int i = 0; i < N - 1; i++) {
// 		y_tilde[i][i + 1] = 1;
// 	}
// 	y_tilde[N - 1][0] = 1;


// //	double val_MP,
// //	double objVal, v_Val;
// 	stef.solve_LP_TSP();
// 	stef.get_optimal_sol(y_star);

// 	/*
// 	for (int i = 0; i < N; i++) {
// 		for (int j = 0; j < N; j++) {
// 			cout << y_star[i][j] << '\t';
// 		}
// 		cout << '\n';
// 	}
// 	cout << "--------------------------------------------" << endl;

// 	*/
// 	double alpha = 0.2;
// 	double UB = INFINITY;
// 	double LB = -INFINITY;

// 	double dual_obj;
// 	double old_LB = -1;
// 	int iteration = 0;
// 	bool is_disconnected = true;


// 	while (is_disconnected || UB - LB > 0.0000001) {
// 		cout << "UB: " << UB << "   " << "LB: " << LB << '\n';

// 		is_disconnected = stef.add_SECs(y_star); //if one subtour elimi constrait is added, flag is
// 		if (!is_disconnected) {
// 			for (int i = 0; i < N; i++) { // update the stablizer
// 				for (int j = 0; j < N; j++) {
// 					y_tilde[i][j] = 0.5 * (y_star[i][j] + y_tilde[i][j]);
// 				}
// 			}
// 			for (int i = 0; i < N; i++) { // new solution for the dual model
// 				for (int j = 0; j < N; j++) {
// 					y_bar[i][j] = alpha * y_star[i][j] + (1 - alpha) * y_tilde[i][j];
// 				}
// 			}
// 			dual_obj = stef.add_USER_cuts(y_bar);
// 		}

// 		pair<double, double> result_tsp =stef.solve_LP_TSP();
// 		stef.get_optimal_sol(y_star);

// 		if (!is_disconnected) {
// 			LB = result_tsp.first;
// 			if (old_LB == LB) {
// 				iteration++;
// 				if (iteration >= 5) {
// 					alpha = 1;
// 				}
// 			}
// 			old_LB = LB;
// 			UB = min(UB, result_tsp.first - result_tsp.second + dual_obj);

// 		}
// 	}
// 	cout << "UB: " << UB << "   " << "LB: " << LB << '\n';

// 	cout << "Gap is closed!" << '\n';


// 	stef.set_model_MIP(); // solve IP-TSP
// 	stef.solve_IP_TSP(); // add Benders cuts

// }
