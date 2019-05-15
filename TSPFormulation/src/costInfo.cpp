#include <iostream>
#include <math.h>
#include "costInfo.h"
#include <algorithm>

using namespace std;

costInfo::costInfo(const char* filename) {
	this->filename = filename;
	fstream file(filename);
	if (!file) {
		cerr << "ERROR: could not open file '" << filename
			<< "' for reading" << endl;
		throw (-1);
	}
	file >> n;
	x = new double[n];
	y = new double[n];
	for (int i = 0; i < n; i++) {
		file >> x[i] >> y[i];
	}
	/* Cost Matrix */
	cost = new double*[n];
	for (int i = 0; i < n; i++){
		cost[i] = new double[n];
	}

	for (int i = 0; i < n; i++){
		for (int j = i + 1; j < n; j++) {
			cost[i][j] = sqrt(pow(x[i] - x[j], 2) + pow(y[i] - y[j], 2));
			cost[j][i] = cost[i][j];
		}
	}
}

double costInfo::getCost(int i, int j){
    if (i<j)
        return cost[i][j];
    else if(j<i)
        return cost[j][i];
    else
        return 0;
}

void costInfo::print(){
    for (int i=0; i<n; i++)
        for (int j=i+1; j<n; j++){
            cout<<i<<" "<<j<<" "<<cost[i][j]<<endl;
        }
}

costInfo::~costInfo(){

    for (int i=0; i<n; i++){
        delete(cost[i]);
    }
    delete x;
    delete y;
}
