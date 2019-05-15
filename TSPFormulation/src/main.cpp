#include "gurobi_c++.h"
#include "costInfo.h"
#include "STEFormulation.h"

int main( int argc, char * argv[]){
	argc = argc;
	const char *  filename =  argv[1];
	costInfo cost(filename);
	/* build the env&model in gurobi */
	GRBEnv* env = new GRBEnv();
	GRBModel model = GRBModel(*env);
	model.getEnv().set(GRB_DoubleParam_TimeLimit,1000);
  /*run the instance */
	STEFormulation STEForm(&cost, &model);
	STEForm.solve(&model);
	STEForm.writeSol(&model);
//	STEForm.printSol(&model);
	return 0;
}
