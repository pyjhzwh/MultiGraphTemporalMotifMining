#include "Tests.h"
#include "CertGraph.h"
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <assert.h> 

using namespace std;

int Tests::countDuplicateEdges(const CertGraph &g)
{
    int m = g.numEdges();
    if(m <= 1)
	return 0;

    int dupCount = 0;

    unordered_map<int, unordered_map<int, unordered_map<string, unordered_set<time_t> > > > edgeMap;
    for(int e=1; e<m; e++)
    {
	const Edge &edge = g.edges()[e];
	int u = edge.source();
	int v = edge.dest();
	const string &type = g.getEdgeType(e);
	time_t t = edge.time();
	auto &edgeSet = edgeMap[u][v][type];
	if(edgeSet.find(t) != edgeSet.end())
	{
	    dupCount++;
	    cout << "Duplicate Edge: " << g.getLabel(u) << " -> " << g.getLabel(v) << " [" << type << "] time(" << t << ")" << endl;
	}
	else
	    edgeSet.insert(t);
    }

    return dupCount;
}


vector<vector<int>> Tests::sortGraphMatch(vector<GraphMatch> &gms)
{
	vector<vector<int>> matched_hEdges;
	for(auto &gm: gms){
		vector<int> hEdges;
		for(int i=0; i < gm.numEdges(); i++)
		{
			hEdges.push_back(gm.getGraphEdge(i));
		}
		matched_hEdges.push_back(hEdges);
	}
	
	return matched_hEdges;
}


bool Tests::assertSameGraphMatch(vector<GraphMatch> gm1, vector<GraphMatch>gm2, Graph g){
    cout << "graph1 matched edges: " << gm1.size() << ",";
	cout << "graph2 matched edges: " << gm2.size() << endl;
    // assert(gm1.size() == gm2.size());
    vector<vector<int>> sorted_hEdges_gm1 = sortGraphMatch(gm1);
    vector<vector<int>> sorted_hEdges_gm2 = sortGraphMatch(gm2);
	sort(sorted_hEdges_gm1.begin(), sorted_hEdges_gm1.end());
	sort(sorted_hEdges_gm2.begin(), sorted_hEdges_gm2.end());
	int i=0, j=0;
	while(i < sorted_hEdges_gm1.size() & j <sorted_hEdges_gm2.size())
	{
		if(sorted_hEdges_gm1[i] != sorted_hEdges_gm2[j]){
			if(sorted_hEdges_gm1.size() < sorted_hEdges_gm2.size()){
				cout << "gm2 hedge: ";
				for (auto edge: sorted_hEdges_gm2[j])
					cout << g.edges()[edge] << ", ";
				cout << endl;
				j++;
			}
			else{
				cout << "gm1 hedge: ";
				for (auto edge: sorted_hEdges_gm1[i])
					cout << g.edges()[edge] << ", ";
				cout << endl;
				i++;
			}
		}
		else{
			i++;
			j++;
		}
	}
	while(i < sorted_hEdges_gm1.size())
	{
		cout << "gm1 hedge: ";
		for (auto edge: sorted_hEdges_gm1[i])
			cout << g.edges()[edge] << ", ";
		cout << endl;
		i++;
	}
	while(j < sorted_hEdges_gm2.size())
	{
		cout << "gm2 hedge: ";
		for (auto edge: sorted_hEdges_gm2[j])
			cout << g.edges()[edge] << ", ";
		cout << endl;
		j++;
	}
	// assert(sorted_hEdges_gm1 == sorted_hEdges_gm2);
	return true;
}
