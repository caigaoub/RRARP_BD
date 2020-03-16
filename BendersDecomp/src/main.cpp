#include <stdlib.h>
#include "gurobi_c++.h"
#include "DataHandler.h"
#include "PartitionScheme.h"
#include "STEFormulation.h"
#include "TSPModels.h"
#include "SuperCutFormulation.h"
#include "BendersCuts.h"
#include "NumericalTest.h"
#include <ctime>
#include <chrono>
#include <boost/lexical_cast.hpp>

using namespace std;

int main(int argc, const char* argv[]) {
	argc = argc; // get rid of warning: unused argc
	int which_BDCut = atoi(argv[1]);
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
		compare_tspSol_vs_optSol(network_);
		exit(0);

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
