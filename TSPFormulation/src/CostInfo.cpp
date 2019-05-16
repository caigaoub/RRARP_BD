#include <iostream>
#include <math.h>
#include "CostInfo.h"
#include <algorithm>

using namespace std;

CostInfo::CostInfo(DataHandler* dh_var) {
	this->DH = dh_var;
	this->n = DH->get_num_targets() + 2;
	x = new double[n];
	y = new double[n];
	for (int i = 0; i < n; i++) {
		if(i == 0){
			x[i] = DH->get_depot1_loc().x;
			y[i] = DH->get_depot1_loc().y;
		}
		if (i == n-1){
			x[i] = DH->get_depot2_loc().x;
			y[i] = DH->get_depot2_loc().y;
		}
		if(i != 0 && i != n-1){
			x[i] = DH->get_target_locs()[i].x;
			y[i] = DH->get_target_locs()[i].y;
		}
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

double CostInfo::getCost(int i, int j){
    if (i<j)
        return cost[i][j];
    else if(j<i)
        return cost[j][i];
    else
        return 0;
}

void CostInfo::print(){
    for (int i=0; i<n; i++)
        for (int j=i+1; j<n; j++){
            cout<<i<<" "<<j<<" "<<cost[i][j]<<endl;
        }
}

CostInfo::~CostInfo(){
    for (int i=0; i<n; i++){
        delete(cost[i]);
    }
    delete x;
    delete y;
}
