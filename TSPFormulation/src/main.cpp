#include <stdlib.h>
#include "gurobi_c++.h"
#include "CostInfo.h"
#include "STEFormulation.h"
#include "DataHandler.h"
#include "myNameClass.h"
#include "PartiScheme.h"

using namespace std;
int main(int argc, char * argv[]){
	argc = argc;
	const char *  filename =  argv[1];
	int k = stoi(argv[2]);
//	const char* filename = "RRARP_instance_n_7_E_1.txt";
	DataHandler instance(filename);
  CostInfo cost(&instance);

	/* build the env&model in gurobi */
	GRBEnv* env = new GRBEnv();
	GRBModel model = GRBModel(*env);
	model.getEnv().set(GRB_DoubleParam_TimeLimit,1000);
	model.getEnv().set(GRB_IntParam_OutputFlag, 0);
   /*run the instance */
	STEFormulation STEForm(&cost, &model);
	vector<int> fseq(instance.get_num_targets() + 2);
	STEForm.solve(&model, fseq);

	PartiScheme ps(k, instance);
	vector<vector<double>> SDS(instance.get_num_targets() + 2);
	double val = ps.solve_shortestpath(SDS, fseq);
	cout << "RRARP tour by TSP: " << val << endl;
	fstream fs;
	fs.open("./out/tsp_objval.dat", fstream::app);
	fs << val << '\n';
	fs.close();
//	system("pause");
	return 0;
}
