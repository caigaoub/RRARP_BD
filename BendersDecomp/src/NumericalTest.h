#include "TSPModels.h"


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


void compare_tspSol_vs_optSol(PartitionScheme & network_, int fischetti_on){
   /* find TSP solution*/
	network_.calc_risk_C2C();
	TSPModel_STE tsp_formul;
	tsp_formul.init_edge_weights(network_._dataset->_nb_targets+2, network_._risk_C2C);
	tsp_formul.create_formula();
	tsp_formul.solve();
    vector<int>& tspSeq = tsp_formul._opt_seq;
    vector<vector<double>>  SDS(network_._dataset->_nb_targets + 2);
    network_.solve_shortestpath_v2(tspSeq, SDS);
    double risk_tspSol = network_.calc_withdrawal_risk(tspSeq) + SDS[network_._dataset->_nb_targets+1][0]; 

	/**find R2ARP solution */
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
	formul_supercut_.add_model(&model_supercut_, network_._dataset->_nb_targets+2);

	formul_master.add_dualformul(&formul_dual_);
	formul_master.add_SuperCutformul(&formul_supercut_);


    int which_BDCut = 3; // using strong cut
	if(fischetti_on == 1){
		Fischetti_method(network_._dataset->_nb_targets+2, formul_master);
	}
	auto ret = formul_master.solve_formul_wCB(which_BDCut);
    vector<int>& optSeq = formul_master._opt_seq;


    vector<vector<double>>  SDS2(network_._dataset->_nb_targets + 2);
    network_.solve_shortestpath_v3(optSeq, SDS2);

    double risk_optSol = ret.first;


	delete evn_MP_;
	delete evn_dual_;
	delete evn_supercut_;
    cout << "Instance name: " << network_._dataset->_name << endl;
    cout << "TSP: " << risk_tspSol << endl;
    cout << "Seq: " << endl;
    for(unsigned int i = 0; i < tspSeq.size(); i++){
        cout << tspSeq.at(i) << ' ';
    }
    cout << '\n';
    cout << "OPT: " << risk_optSol << endl;
    cout << "Seq: " << endl;
    for(unsigned int i = 0; i < tspSeq.size(); i++){
        cout << optSeq.at(i) << ' ';
    }
    cout << '\n';
    cout << "Gap: " << (risk_tspSol - risk_optSol)/risk_optSol << '\n';
}

void test_K(PartitionScheme & network_){
	for(int i = 4; i < 16; i++){
		;
	}
}

void compare_dualcut_vs_shortestpathcut(PartitionScheme & network_){

	// vector<int> size_testInst = {6, 8, 10}


}