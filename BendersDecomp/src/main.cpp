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
void Fischetti_method(int , STEFormulation & );



int main(int argc, const char* argv[]) {

	argc = argc; // just for avoid warning: unused argc
	const int nb_dstzn = atoi(argv[1]);
	string filename = argv[2];   
	// int Fischetti_on = atoi(argv[3]); 
	int which_cut = atoi(argv[3]);
	try {
		DataHandler dataset_;
		dataset_.parse(filename);
		// dataset_.print();
		PartitionScheme network_;
		network_.build(dataset_, nb_dstzn);
		
		/* test */
		// vector<int> fseq_;
		// fseq_.push_back(0);
		// fseq_.push_back(3);
		// fseq_.push_back(4);
		// fseq_.push_back(5);
		// fseq_.push_back(1);
		// fseq_.push_back(2);
		// fseq_.push_back(6);
		// vector<vector<double>> * SDS = new vector<vector<double>>(dataset_._nb_targets + 2);
		// network_.solve_shortestpath_v2(fseq_, *SDS);

		GRBEnv * evn_MP_ = new GRBEnv();
		GRBModel model_MP_ = GRBModel(*evn_MP_);
		STEFormulation formul_master;
		formul_master.build_formul(&model_MP_, &network_);


		/* Create dual formulation*/ 
		GRBEnv * evn_dual_ = new GRBEnv();
		GRBModel model_dual_ = GRBModel(*evn_dual_);
		model_dual_.getEnv().set(GRB_IntParam_OutputFlag, 0);
		DualFormulation formul_dual_;
		formul_dual_.create_variables(&model_dual_, &network_);
		formul_dual_.set_constraints();

		formul_master.add_dualformul(&formul_dual_);
	
		auto start = chrono::system_clock::now();
		formul_master.solve_formul_wCB(which_cut);
		auto end = chrono::system_clock::now();
		formul_master.print_solution(&model_MP_);
		chrono::duration<double> elapsed_seconds = end-start;
		cout << "====>>> Algorithm: " << which_cut << " time: " << std::chrono::duration<double>(elapsed_seconds).count()  << endl;		

		
		// if(Fischetti_on){
		// 	auto start_fischetti= chrono::system_clock::now();
		// 	Fischetti_method(dataset_._nb_targets + 2, formul_master);
		// 	auto endt_fischetti = chrono::system_clock::now();
		// 	chrono::duration<double> elapsed_seconds_fischetti = endt_fischetti-start_fischetti;
		// 	cout << "====>>> Fischetti_method: " << std::chrono::duration<double>(elapsed_seconds_fischetti).count()  << endl;		
		// }else{
		// 	auto start_no_fischetti= chrono::system_clock::now();
		// 	formul_master.set_vars_integer();
		// 	formul_master.solve_formul_woCB();
		// 	auto endt_no_fischetti = chrono::system_clock::now();
		// 	chrono::duration<double> elapsed_seconds_no_fischetti = endt_no_fischetti-start_no_fischetti;
		// 	cout << "====>>> NO Fischetti_method: " << std::chrono::duration<double>(elapsed_seconds_no_fischetti).count()  << endl;		
		// }
		
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


void Fischetti_method(int N, STEFormulation & formul_master) {

	/* create variables */
	double** y_tilde = new double*[N]; // stablizer point
	double** y_star = new double*[N]; // optimal solution of current TSP_LP
	double** y_bar = new double*[N];
	for (int i = 0; i < N; i++) {
		y_tilde[i] = new double[N];
		y_star[i] = new double[N];
		y_bar[i] = new double[N];
	}

	/*initialize y_tilde to represent the feasible sequence (0 -> 1 -> 2 ->...->N-1) */
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			y_tilde[i][j] = 0;
		}
	}
	for (int i = 0; i < N - 1; i++) {
		y_tilde[i][i + 1] = 1;
	}
	y_tilde[N - 1][0] = 1;


	formul_master.set_vars_continuous(); // solve LP-TSP model
	formul_master.solve_formul_woCB();
	formul_master.get_optimal_sol(y_star);

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

		is_disconnected = formul_master.add_SECs(y_star); //if one subtour elimi constrait is added, flag is
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
			dual_obj = formul_master.add_USER_cuts(y_bar);
		}

		pair<double, double> result_tsp =formul_master.solve_formul_woCB();
		formul_master.get_optimal_sol(y_star);

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


	formul_master.set_vars_integer(); // solve IP-TSP
	formul_master.solve_formul_woCB(); // add Benders cuts

}
