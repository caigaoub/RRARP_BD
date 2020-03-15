#include <stdlib.h>
#include "gurobi_c++.h"
#include "DataHandler.h"
#include "PartitionScheme.h"
#include "STEFormulation.h"
#include "TSPModels.h"
#include "SuperCutFormulation.h"
#include "BendersCuts.h"
#include <ctime>
#include <chrono>
#include <boost/lexical_cast.hpp>

using namespace std;
pair<int,int> Fischetti_method(int , STEFormulation & );
pair<int,int> improve_root(int N, STEFormulation & formul_master);


int main(int argc, const char* argv[]) {
	argc = argc; // just for avoid warning: unused argc
	int which_BDCut = atoi(argv[1]);
	// bool fischetti_on = boost::lexical_cast<bool>(argv[2]);
	const int fischetti_on = atoi(argv[2]);
	const int nb_dstzn = atoi(argv[3]);
	const int type_trajc = atoi(argv[4]);
	string configfile = argv[5];

	fstream file(configfile);
	if (!file) {
		cerr << "ERROR: could not open config '" << configfile << "' for reading'" << endl;
		throw(-1);
	}
	string instance_name_only;
	file >> instance_name_only;
	file.close();
	try {
		//	string cur_dir  = boost::filesystem::current_path().string();
		//	auto pos = cur_dir.find_last_of("/");
    	//      cur_dir = cur_dir.substr(0, pos);
		// string cur_dir = "/projects/academic/josewalt/caigao/RRARP_BD/BendersDecomp/dat/";
		string cur_dir = "/home/caigao/Dropbox/Box_Research/Github/RRARP_BD/BendersDecomp/dat/";
		// string cur_dir = "/home/cai/Dropbox/Box_Research/Github/RRARP_BD/BendersDecomp/dat/";
		struct stat buffer;
	  	if(stat (cur_dir.c_str(), &buffer) != 0){
	  		cerr << " Path of instances does not exist!! (in main.cpp:line 42) " << endl;
	  	}
		string instance_wPath = cur_dir + instance_name_only;
		DataHandler dataset_;
		dataset_.parse(instance_wPath);
		// dataset_.print();
		PartitionScheme network_;
		network_.build(dataset_, nb_dstzn, type_trajc);
		if(false){ 
			/* solve the TSP sequence */
			network_.calc_risk_C2C();
			TSPModel_STE tspsol;
			tspsol.init_edge_weights(network_._dataset->_nb_targets+2, network_._risk_C2C);
			tspsol.create_formula();
			tspsol.solve();
		}
		// exit(0);
		/*Gurobi model for master problem */
		GRBEnv * evn_MP_ = new GRBEnv();
		GRBModel model_MP_ = GRBModel(*evn_MP_);
		STEFormulation formul_master;
		// model_MP_.getEnv().set(GRB_IntParam_OutputFlag, 0);

		formul_master.build_formul(&model_MP_, &network_);

		/*Gurobi model for dual formulation */ 
		GRBEnv * evn_dual_ = new GRBEnv();
		GRBModel model_dual_ = GRBModel(*evn_dual_);
		model_dual_.getEnv().set(GRB_IntParam_OutputFlag, 0);
		DualFormulation formul_dual_;
		formul_dual_.create_variables(&model_dual_, &network_);
		formul_dual_.set_constraints();

		/*Gurobi model for SuperCut formulation */
		GRBEnv * evn_supercut_ = new GRBEnv();
		GRBModel model_supercut_ = GRBModel(*evn_supercut_);
		model_supercut_.getEnv().set(GRB_IntParam_OutputFlag, 0);
		SuperCutFormulation formul_supercut_;
		formul_supercut_.add_model(&model_supercut_, dataset_._nb_targets+2);

		formul_master.add_dualformul(&formul_dual_);
		formul_master.add_SuperCutformul(&formul_supercut_);

		if(fischetti_on == 1){
			auto start_fischetti= chrono::system_clock::now();
			Fischetti_method(dataset_._nb_targets + 2, formul_master);
			auto endt_fischetti = chrono::system_clock::now();
			chrono::duration<double> elapsed_seconds_fischetti = endt_fischetti-start_fischetti;
			cout << "====>>> Total time of Fischetti_method: " << std::chrono::duration<double>(elapsed_seconds_fischetti).count()  << endl;		
		}
		if(fischetti_on == 2){
			auto start_no_fischetti= chrono::system_clock::now();
			improve_root(dataset_._nb_targets + 2, formul_master);
			auto endt_no_fischetti = chrono::system_clock::now();
			chrono::duration<double> elapsed_seconds_no_fischetti = endt_no_fischetti-start_no_fischetti;
			cout << "====>>> Total time of NO Fischetti_method: " << std::chrono::duration<double>(elapsed_seconds_no_fischetti).count()  << endl;		
		}

		if(true){
			formul_master.solve_formul_wCB(which_BDCut);
			// formul_master.print_solution();
			formul_master.write_solution(dataset_._name, which_BDCut);
		}
		
		delete evn_MP_;
		delete evn_dual_;
		delete evn_supercut_;
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


pair<int,int> Fischetti_method(int N, STEFormulation & formul_master) {

	int nb_Subtour_Cuts = 0;
	int nb_USER_Cuts = 0;
	auto s_= chrono::system_clock::now();
	// auto e_= chrono::system_clock::now();
	chrono::duration<double> total_time_SP = s_ - s_;
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


	formul_master.set_vars_continuous(); // solve linear relaxation 
	formul_master.solve_formul_woCB();
	formul_master.get_optimal_sol(y_star);

	double alpha = 0.2;
	// double UB = INFINITY;
	double LB = -INFINITY;

	double dual_obj = 0;
	double old_LB = -1;
	int iteration = 0;
	pair<bool,int> ret_SEC;
	pair<double, double> result_tsp = make_pair(INFINITY, INFINITY);

	while (true ) {//is_disconnected.first || UB - LB > 0.0000001
		// cout << "*UB: " << UB << "   " << "LB: " << LB << '\n';

		ret_SEC = formul_master.add_SECs(y_star); // add subtour elimi constraint
		nb_Subtour_Cuts += ret_SEC.second;
		if (ret_SEC.first) {
			for (int i = 0; i < N; i++) { // update the stablizer
				for (int j = 0; j < N; j++) {
					y_tilde[i][j] = 0.5 *(y_star[i][j] + y_tilde[i][j]);
					// y_tilde[i][j] = y_star[i][j];
				
				}
			}
			for (int i = 0; i < N; i++) { // new solution for the dual model
				for (int j = 0; j < N; j++) {
					y_bar[i][j] = alpha * y_star[i][j] + (1 - alpha) * y_tilde[i][j];
				}
			}
			auto start_= chrono::system_clock::now();
			dual_obj = formul_master.add_USER_cuts(y_bar);
			auto end_= chrono::system_clock::now();
			total_time_SP += end_ - start_;			
			nb_USER_Cuts++;
		}

		auto start_= chrono::system_clock::now();
		result_tsp =formul_master.solve_formul_woCB();
		auto end_= chrono::system_clock::now();
		total_time_SP += end_ - start_;

		formul_master.get_optimal_sol(y_star);

		if (ret_SEC.first) {
			LB = result_tsp.first;
			if (old_LB == LB) {
				iteration++;
				if (iteration >= 5) {
					alpha = 1;
				}
			}
			old_LB = LB;
			if(abs(result_tsp.second - dual_obj) < 0.001){
				break;
			}
		}
	}
	// cout << "UB: " << UB << "   " << "LB: " << LB << '\n';
	cout << "====>>> FISCHETTI GIVES: " << endl;
	cout << "====>>> total time of solving dual problem: " << chrono::duration<double>(total_time_SP).count()  << endl;		
	cout << "====>>> nb of subtour cuts added: " << nb_Subtour_Cuts  << endl;		
	cout << "====>>> nb of user cuts added: " << nb_USER_Cuts << endl;
	cout << "====>>> optimal: " << LB << endl;
	// cout << "====>>> Time: " << formul_master._model->get(GRB_DoubleAttr_Runtime) << endl;
	formul_master.set_vars_integer(); // solve IP-TSP
	// formul_master.solve_formul_woCB(); // add Benders cuts
	return make_pair(nb_Subtour_Cuts, nb_USER_Cuts);
}


pair<int,int> improve_root(int N, STEFormulation & formul_master) {

	int nb_Subtour_Cuts = 0;
	int nb_USER_Cuts = 0;
	auto s_= chrono::system_clock::now();
	// auto e_= chrono::system_clock::now();
	chrono::duration<double> total_time_SP = s_ - s_;
	
	double** y_star = new double*[N]; 
	for (int i = 0; i < N; i++) {
		y_star[i] = new double[N];
	}

	formul_master.set_vars_continuous(); // solve linear relaxation 
	formul_master.solve_formul_woCB();
	formul_master.get_optimal_sol(y_star);

	// double LB = -INFINITY;
	double dual_obj = INFINITY;
	pair<bool,int> ret_SEC;
	pair<double, double> result_tsp;
	while (true) { 
		// cout << "LB: " << LB << '\n';
		// for (int i = 0; i < N; i++) {
		// 	for(int j = 0; j < N; j++){
		// 		cout << y_star[i][j] << " ";
		// 	}
		// 	cout << endl;
		// }

		ret_SEC = formul_master.add_SECs(y_star); // add subtour elimi constraint
		// cout << "connected: " <<  ret_SEC.first << endl;
		nb_Subtour_Cuts += ret_SEC.second;
		if (ret_SEC.first) {
			auto start_= chrono::system_clock::now();
			dual_obj = formul_master.add_USER_cuts(y_star);
			// exit(0);
			auto end_= chrono::system_clock::now();
			total_time_SP += end_ - start_;
			nb_USER_Cuts++;
			// cout << "adding user cuts !!!!" << endl;
		}

		auto start_= chrono::system_clock::now();
		result_tsp =formul_master.solve_formul_woCB();
		auto end_= chrono::system_clock::now();
		total_time_SP += end_ - start_;

		formul_master.get_optimal_sol(y_star);

		// LB = result_tsp.first;
			// old_LB = LB;
		// cout << "LB: " << result_tsp.first << " ";
		// cout << result_tsp.second << " " << dual_obj << endl;
		if(abs(result_tsp.second - dual_obj) < 0.001){
			break;
		}
	}
	// cout << "UB: " << UB << "   " << "LB: " << LB << '\n';
	cout << "====>>> NO FISCHETTI GIVES: " << endl;
	cout << "====>>> total time of solving dual problem: " << chrono::duration<double>(total_time_SP).count()  << endl;		
	cout << "====>>> nb of subtour cuts added: " << nb_Subtour_Cuts  << endl;		
	cout << "====>>> nb of user cuts added: " << nb_USER_Cuts << endl;
	cout << "====>>> optimal: " << result_tsp.first << endl;
	// cout << "====>>> Time: " << formul_master._model->get(GRB_DoubleAttr_Runtime) << endl;

	formul_master.set_vars_integer(); // solve IP-TSP
	// formul_master.solve_formul_woCB(); // add Benders cuts
	return make_pair(nb_Subtour_Cuts, nb_USER_Cuts);
}
