
#ifndef _COSTINFO_H_
#define _COSTINFO_H_

#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include "DataHandler.h"

using namespace std;

class CostInfo{

private:
    DataHandler * DH;
    int n;// number of nodes
    double* x;
    double* y;
    double** cost;
public:
    ~CostInfo();
    CostInfo(DataHandler*);
	  inline double** getCost() { return cost;};
    inline int getNumNodes(){return n;};
    double getCost (int i, int j);
    void print();
    inline static string itos(int i) {stringstream s; s << i; return s.str();};
    inline static string dtos(double i) {stringstream s; s << i; return s.str();};
};

#endif
