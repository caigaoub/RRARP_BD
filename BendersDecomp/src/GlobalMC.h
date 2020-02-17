#pragma once
#include <vector>
#include <iostream>
#include <set>
#include <map>
#include <unordered_set>
#include "Edge.h"
#include <tuple>
#include <boost/heap/fibonacci_heap.hpp>

class GlobalMC
{
public:
	unsigned int mc_n;
	unsigned int mc_m;
	std::set<int> mcVers;
	std::map<std::pair<int, int>, Edge*> Edges;
	std::map<std::pair<int, int>, Edge*> mcEdges;
	std::map<std::pair<int, int>, std::vector<std::pair<int, int>>> mc_mrg_edges;
	std::map<int, std::set<int>> mcAdjVers;
	std::vector<Edge> oriEdges;
	std::vector<Edge> copyEdges;
	std::vector<std::vector<int>> nodePartition;

	std::tuple<double, std::unordered_set<Edge*>, std::vector<int>> GlobalMinCut();
	void MergeVers(int v_mrg, int u_del);
	std::pair<int, int> getKey(int v1, int v2);

	GlobalMC(int n, int m, double** x);
	~GlobalMC();
};
