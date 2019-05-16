#include <stdlib.h>
#include "gurobi_c++.h"
#include "CostInfo.h"
#include "STEFormulation.h"
#include "DataHandler.h"
#include "myNameClass.h"
#include "PartiScheme.h"

using namespace std;
int main( int argc, char * argv[]){
	argc = argc;
	const char *  filename =  argv[1];
	int k = stoi(argv[2]);
	DataHandler instance(filename);
  CostInfo cost(&instance);
	/* build the env&model in gurobi */
	GRBEnv* env = new GRBEnv();
	GRBModel model = GRBModel(*env);
	model.getEnv().set(GRB_DoubleParam_TimeLimit,1000);
  /*run the instance */
	STEFormulation STEForm(&cost, &model);
	vector<int> fseq(instance.get_num_targets() + 2);
	STEForm.solve(&model, fseq);

	PartiScheme ps(k, instance);
	vector<vector<double>> SDS;
	double val = ps.solve_shortestpath(SDS, fseq);
	cout << "RRARP tour by TSP: " << val << endl;
	return 0;
}
