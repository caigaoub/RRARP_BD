#include <iostream>
#include "gurobi_c++.h"
#include "InstanceGenerator.h"

using namespace std;
int main(int argc, char* argv[]) {
  cout << "Instance Generator: " << endl;
  int num_targets = atoi(argv[1]);
  int difflevel = atoi(argv[2]);

  GRBEnv* evn = new GRBEnv();
  GRBModel model = GRBModel(*evn);
  InstanceGenerator sample(num_targets, &model);
  sample.produce(difflevel);

  cout << "This is the end. Thank you!" << endl;
  return 0;

}
