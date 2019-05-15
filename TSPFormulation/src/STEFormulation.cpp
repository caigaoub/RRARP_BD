#include "gurobi_c++.h"
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <fstream>
#include "STEFormulation.h"
#include "costInfo.h"
using namespace std;
extern void findsubtour(int  n, double** sol, int*  tourlenP, int*  tour);
STEFormulation::STEFormulation(costInfo* c, GRBModel * model)
{
	this->cost = c;
	this->filename = c->getFileName();
	N = cost->getNumNodes();
	// set up the model
	model->set(GRB_IntAttr_ModelSense, 1);
	model->getEnv().set(GRB_IntParam_LazyConstraints, 1);
	model->getEnv().set(GRB_IntParam_PreCrush, 1);
	model->getEnv().set(GRB_IntParam_MIPFocus, 1); // focus on feasibility
  int i, j;
	// add variables
	y = new GRBVar*[N];
	for (i = 0; i < N; i++)
		y[i] = new GRBVar[N];

	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			y[i][j] = model->addVar(0.0, 1.0, cost->getCost(i,j), GRB_BINARY, "y_" + itos(i) + "_" + itos(j));
		}
	}
	model->update();

	// Degree-1 constraints
	for (i = 0; i < N; i++)	{
		GRBLinExpr expr = 0;
		for (j = 0; j < N; j++)	{
			expr += 1.0* y[i][j];
		}
		model->addConstr(expr == 1, "deg1_row" + itos(i));
		GRBLinExpr expr2 = 0;
		for (j = 0; j < N; j++)	{
			expr2 += 1.0* y[j][i];
		}
		model->addConstr(expr2 == 1, "deg1_col" + itos(i));
	}

	// Forbid edge from node back to itself
	for (i = 0; i < N; i++)	{
		y[i][i].set(GRB_DoubleAttr_UB, 0);
	}
	model->update();
	y[N - 1][0].set(GRB_DoubleAttr_LB, 1);
	model->update();
}

void STEFormulation::solve(GRBModel *model)
{
	try
	{
		SubtourElimCuts cb = SubtourElimCuts(y, N);
		model->setCallback(&cb);
		// Optimize model
		model->optimize();
		numSubtourConst = cb.get_numSubtourCuts();
		cout << numSubtourConst << endl;
		if (model->get(GRB_IntAttr_Status) == GRB_OPTIMAL)
			status = 0;
	}
	catch (GRBException e) {
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during optimization" << endl;
	}
}

void STEFormulation::writeSol(GRBModel * model){
	if (model->get(GRB_IntAttr_SolCount) > 0)
	{
		fstream fs;
		fs.open("./out/TSPSols.dat", fstream::app);
	//	ofstream TSPSols("TSPSols.dat");
		fs << filename << ":" << '\t';
		double **sol = new double*[N];
		int i;
		for (i = 0; i < N; i++)
			sol[i] = model->get(GRB_DoubleAttr_X, y[i], N);

		int* tour = new int[N];
		int len;
		findsubtour(N, sol, &len, tour);
		for (i = 0; i < len; i++)
				fs << tour[i]<< '\t';
		fs << '\n';

		for (i = 0; i < N; i++)
			delete[] sol[i];
		delete[] sol;
		delete[] tour;
	}
}



void STEFormulation::printSol(GRBModel *model) {
	if (model->get(GRB_IntAttr_SolCount) > 0)
	{
		double **sol = new double*[N];
		int i;
		for (i = 0; i < N; i++)
			sol[i] = model->get(GRB_DoubleAttr_X, y[i], N);

		int* tour = new int[N];
		int len;

		findsubtour(N, sol, &len, tour);

		cout << "Tour: ";
		for (i = 0; i < len; i++)
			cout << tour[i] << " ";
		cout << endl;

		for (i = 0; i < N; i++)
			delete[] sol[i];
		delete[] sol;
		delete[] tour;
	}
}
STEFormulation::~STEFormulation(){
	for (int i = 0; i < N; i++)
		delete[] y[i];
	delete[] y;
}
