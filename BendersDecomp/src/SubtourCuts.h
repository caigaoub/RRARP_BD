#ifndef _SUBTOURCUTS_H_
#define _SUBTOURCUTS_H_
#include <cmath>
#include <stdlib.h>
#include <queue>
#include "gurobi_c++.h"
#include "PartitionScheme.h"
#include "BendersCuts.h"
#include <sstream>
using namespace std;

class SubtourCuts : public GRBCallback {
private:
	int 						_size_var_x;
	GRBVar*** 					_var_x;
public:
	SubtourCuts(GRBVar***, int);
	string itos(int i) { stringstream s; s << i; return s.str(); }
protected:
	void callback();
};


#endif // !_SUBTOURCUTS_H_
