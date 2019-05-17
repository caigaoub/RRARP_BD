#include "GlobalMC.h"
#include <float.h>
#include<bits/stdc++.h>

using namespace std;

GlobalMC::GlobalMC(int n, int m, double** x)
{
	mc_n = n;
	mc_m = m;
	nodePartition.resize(n);
	for (unsigned int i = 0; i < mc_n; i++) {
		mcVers.insert(i);
		nodePartition[i].push_back(i);
	}
	oriEdges.resize(m);
	copyEdges.resize(m);
	int k;
	k = 0;
	for (unsigned int i = 0; i < mc_n; i++) {
		set<int> tmpNeib;
		mcAdjVers.insert(pair<int, set<int>>(i, tmpNeib));
		for (unsigned int j = 0; j < mc_n; j++) {
			if (x[i][j] >= 0) {
				mcAdjVers[i].insert(j);
				if (i < j) {
					Edge e_0(i, j, x[i][j]);
					Edge e_1(i, j, x[i][j]);
					e_0.eid = k;
					e_1.eid = k;
					oriEdges.at(k) = e_0;
					copyEdges.at(k) = e_1;
					Edges.insert(pair<pair<int, int>, Edge*>(pair<int, int>(e_0.a, e_0.b), &oriEdges.at(k)));
					mcEdges.insert(pair<pair<int, int>, Edge*>(pair<int,int>(e_0.a,e_0.b), &copyEdges.at(k)));
					k++;
				}
			}
		}
	}

}

pair<int, int> GlobalMC::getKey(int v1, int v2) {
	pair<int, int> eKey;
	if (v1 < v2) {
		eKey.first = v1;
		eKey.second = v2;
	}
	else {
		eKey.first = v2;
		eKey.second = v1;
	}
	return eKey;
}

void GlobalMC::MergeVers(int v_mrg, int u_del) {
	nodePartition[v_mrg].push_back(u_del);
//	int t1, t2;//end point
	bool ifFnd;
	pair<int, int> eKeyDel, eKeyMrg;
	for (int t1 : mcAdjVers.at(u_del)) {
		eKeyDel = getKey(u_del, t1);
		//Case 1: t1==v: delete edge, and adjver for t1
		if (t1 == v_mrg) {
			mcEdges.erase(eKeyDel);
			mcAdjVers.at(t1).erase(u_del);
			mc_m--;
			continue;
		}
		ifFnd = false;
		//case 2: t1 in v.adj: add weight to (v,t1), delete adjver for t1, delete edge
		for (int t2 : mcAdjVers.at(v_mrg)) {
			if (t2 == t1) {
				eKeyMrg = getKey(t2, v_mrg);
				//mc_mrg_edges.at(eKeyMrg).push_back(eKeyDel);
				mc_mrg_edges.at(eKeyMrg).insert(mc_mrg_edges.at(eKeyMrg).begin(), mc_mrg_edges.at(eKeyDel).begin(), mc_mrg_edges.at(eKeyDel).end());
				//mcEdges.at(eKeyMrg)->c += mcEdges.at(eKeyDel)->c;
				mcEdges.at(eKeyMrg)->w += mcEdges.at(eKeyDel)->w;
				mcEdges.erase(eKeyDel);
				mcAdjVers.at(t1).erase(u_del);
				mc_mrg_edges.erase(eKeyDel);
				mc_m--;
				ifFnd = true;
				break;
			}
		}
		//case 3: t1 not in v.adj: add new edge, add adjver for v, delete old edge, delete adjver for u
		if (ifFnd == false) {
			Edge* eTmp;
			eTmp = mcEdges.at(eKeyDel);
			pair<int, int> eKey;

			vector<pair<int, int>> mrg_E = mc_mrg_edges.at(eKeyDel);

			mcEdges.erase(eKeyDel);
			mcAdjVers.at(t1).erase(u_del);
			mc_mrg_edges.erase(eKeyDel);

			int v1 = t1 < v_mrg ? t1 : v_mrg;
			int v2 = t1 > v_mrg ? t1 : v_mrg;

			eTmp->a = v1;
			eTmp->b = v2;

			eKey.first = v1;
			eKey.second = v2;
			//update edges and adjlist

			mcEdges.insert(pair<pair<int, int>, Edge*>(eKey, eTmp));
			mcAdjVers.at(v_mrg).insert(t1);
			mcAdjVers.at(t1).insert(v_mrg);
			mc_mrg_edges.insert(pair<pair<int, int>, vector<pair<int, int>>>(eKey, mrg_E));
		}
	}
	mcAdjVers.erase(u_del);
	mcVers.erase(u_del);
	mc_n--;

}

std::tuple<double, std::unordered_set<Edge*>, std::vector<int>> GlobalMC::GlobalMinCut() {
	unordered_set<Edge*> min_cut;

	//initialize mc_mrg_edges;
	for (pair<pair<int, int>, Edge*> my_e : mcEdges) {
		mc_mrg_edges.insert(pair<pair<int, int>, vector<pair<int, int>>>(my_e.first, vector<pair<int, int>>()));
		mc_mrg_edges.at(my_e.first).push_back(my_e.first);
	}

	set<int> A;
	double gMinCut = DBL_MAX;
	double tmpMinCut = -1;
	vector<int> partition;

	//set<int> inA;
	int u = -1; //neighbor
//	int u1;
	int vm = -1;
//	int k;

	//vector<int> verID;
	vector<vector<pair<int, int>>> mySol;
	vector<pair<int, int>> tmp_cut, tmp_mc;
	map<int, int> revVerID;
	//minCut update phase

	//int cnt = 0;
	//int err_num = 3;
	while (mc_n > 1)
	{
		//cout << ++cnt << endl;
		//if (cnt == err_num){
		//	int tt00 = 0;
		//}
		A.clear(); //clear set A


				   //ini the vector of min-cut
		mySol.clear();
		mySol.assign(mc_n, vector<pair<int, int>>());


		//(v, cost)
		typedef pair<const int, double>* prim_node;

		struct mycomp {
			bool operator() (const prim_node v1, const prim_node v2) const
			{
				return v1->second < v2->second;
			};
		};
		typedef boost::heap::fibonacci_heap<prim_node, boost::heap::compare<mycomp>> my_heap;

		my_heap use_heap;
		map<int, my_heap::handle_type> keys;
		map<int, double> nodes_set;

		//initial the minimal element

		int ini = 0;
		for (set<int>::iterator it = mcVers.begin(); it != mcVers.end(); ++it) {
			if (ini == 0) {
				nodes_set.insert(pair<int, double>(*it, INT_MAX));
				ini = 1;
			}
			else {
				nodes_set.insert(pair<int, double>(*it, 0));
			}
		}

		//initialize the heap
		for (map<int, double>::iterator it = nodes_set.begin(); it != nodes_set.end(); ++it) {
			keys.insert(pair<int, my_heap::handle_type>(it->first, use_heap.push(&(*it))));
		}


		prim_node v_m = NULL;
		prim_node u_m = NULL;
		//minCut calc phase
		while (A.size() <= (mc_n - 1)) {

			//get mincut
			if (A.size() == (mc_n - 1)) {

				u_m = use_heap.top();
				use_heap.pop();
				keys.erase(u_m->first);

				//u_m = PQ->extractMin();
				tmpMinCut = u_m->second;

				u = u_m->first;
				tmp_cut.clear();
				for (int v : mcAdjVers.at(u)) {
					int v1 = v < u ? v : u;
					int v2 = v > u ? v : u;
					pair<int, int> my_p(v1, v2);
					tmp_cut.insert(tmp_cut.end(), mc_mrg_edges.at(my_p).begin(), mc_mrg_edges.at(my_p).end());
				}

				vm = v_m->first;
				break;
			}

			v_m = use_heap.top();
			use_heap.pop();
			keys.erase(v_m->first);

			//A.push_back(v_m->first);
			A.insert(v_m->first);
			//inA.at(v_m->first) = true;
			pair<int, int> eKey;
			//int cnt1 = 0;
			//int err_num1 = 17;
			for (int u_0 : mcAdjVers.at(v_m->first)) {
				//get u
				//cout << "\t"<<++cnt1 << endl;
				//if (cnt == err_num && cnt1 == err_num1){
				//	int tt1 = 0;
				//}
				//if not in set A
				if (A.find(u_0) == A.end()) {
					if (keys.find(u_0) != keys.end()) {
						eKey = getKey(v_m->first, u_0);
						nodes_set.at(u_0) += mcEdges.at(eKey)->w;
						use_heap.increase(keys.at(u_0));
					}

				}

			}

		}
		//delete PQ;
		//merge v_m->v and u
		MergeVers(vm, u);

		//update new global min-cut
		if (tmpMinCut < gMinCut) {
			gMinCut = tmpMinCut;
			tmp_mc = tmp_cut;
			partition = nodePartition[u];
		}
	}
	double check_mc = 0;
	for (unsigned int i = 0; i < tmp_mc.size(); i++) {
		min_cut.insert(Edges.at(tmp_mc.at(i)));
		check_mc += Edges.at(tmp_mc.at(i))->w;
	}

	return tuple<double, unordered_set<Edge*>,vector<int>>(gMinCut, min_cut, partition);
}


GlobalMC::~GlobalMC(){

}
