#include <iostream>
#include "InstanceGenerator.h"

using namespace std;
int main() {

  cout << "Instance Generator: " << endl;
  InstanceGenerator sample(10);
  sample.produce(2);

  cout << "This is the end. Thank you!" << endl;
  return 0;

}
