#ifndef _INSTANCEGENERATOR_H_
#define _INSTANCEGENERATOR_H_

#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
using namespace std;

class costInfo{

private:
    const char* filename;
    int n;// number of nodes
    double* x;
    double* y;
    double** cost;
public:
    ~costInfo();
    costInfo(const char* filename);
	  inline double** getCost() { return cost;};
    inline int getNumNodes(){return n;};
    inline const char* getFileName(){return filename;};
    double getCost (int i, int j);
    void print();
    inline static string itos(int i) {stringstream s; s << i; return s.str();};
    inline static string dtos(double i) {stringstream s; s << i; return s.str();};
};

#endif
