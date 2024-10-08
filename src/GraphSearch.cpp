#include <limits.h>
#include <unordered_set>
#include <unordered_map>
#include <stack>
#include <iostream>
#include "GraphSearch.h"
#include "Graph.h"
#include "Edge.h"
#include <limits.h>
#include "timer.h"
#include <vector>
#include <algorithm>
#include <omp.h>
#include "DataGraph.h"
#include <random>
#include "timer.h"
#include <iostream>
#include <fstream>
#include <queue>

// #define SAVE_SEARCH_TREE
#define MAX_ST_CNT 1000

using namespace std;


// https://stackoverflow.com/questions/9345087/choose-m-elements-randomly-from-a-vector-containing-n-elements
template<class BidiIter >
BidiIter GraphSearch::random_unique(BidiIter begin, BidiIter end, size_t num_random, unsigned int seed) {
    size_t left = std::distance(begin, end);
    while (num_random--) {
        BidiIter r = begin;
        std::advance(r, rand_r(&seed)%left);
        std::swap(*begin, *r);
        ++begin;
        --left;
    }
    return begin;
}

vector<int> GraphSearch::random_select_n(vector<int>& list, size_t num_random, mt19937& gen) {
    vector<int> result;

    uniform_int_distribution<> dis(0, list.size()-1);
    while(num_random--) {
        // https://stackoverflow.com/questions/2129705/why-is-rand-anything-always-0-in-c
        // int rand_num = rand_r(&seed);
        int r = dis(gen);
        result.push_back(list[r]);
        list.erase(list.begin() + r);
    }

    return result;
}

void GraphSearch::DFSUtil(const Graph &h,
    int node_id, // current node id to be visited
    vector<bool>& visited, vector<int>& currentTree,
    vector<vector<int>>& allSpanningTrees)
{
    visited[node_id] = true;
    for(int e_id: h.nodes()[node_id].edges()) // edges that incident to node_id
    {
        const Edge& e = h.edges()[e_id];
        int next_node_id = e.source() == node_id ? e.dest() : e.source();
        if(!visited[next_node_id])
        {
            currentTree.push_back(e_id);
            // must contain the edge 0
            if (currentTree.size() == h.numNodes() - 1 && find(currentTree.begin(), currentTree.end(), 0) != currentTree.end())
            {
                allSpanningTrees.push_back(currentTree);
            }
            else
            {
                DFSUtil(h, next_node_id, visited, currentTree, allSpanningTrees);
            }
            currentTree.pop_back();
        }
    }
    visited[node_id] = false;
}

void GraphSearch::getStatsSPTree(
    const vector<vector<int>>& allSpanningTrees, int numEdges,
    vector<int>& imbalances, vector<int>& increRegions, vector<int>& absDiffSums)
{
    imbalances.resize(allSpanningTrees.size());
    increRegions.resize(allSpanningTrees.size());
    absDiffSums.resize(allSpanningTrees.size());
    int best_idx = -1;
    for(int j = 0; j < allSpanningTrees.size(); j++)
    {
        int increRegion = 0;
        int absDiffSum = 0;
        for(int i = 1 ; i < allSpanningTrees[j].size(); i++) // we prefer the tree with less increase region
        {
            if (allSpanningTrees[j][i] - allSpanningTrees[j][i-1] < 0)
            {
                increRegion++;
            }
            absDiffSum += abs(allSpanningTrees[j][i] - allSpanningTrees[j][i-1] - 1);
        }
        int minDiff = numEdges;
        int maxDiff = 0;
        vector<int> sorted_tree = allSpanningTrees[j];
        std::sort(sorted_tree.begin(), sorted_tree.end());
        sorted_tree.push_back(numEdges - 1);
        for(int i = 1; i < sorted_tree.size(); i++)
        {
            int diff = sorted_tree[i] - sorted_tree[i-1] - 1;
            if (diff < 0)
                continue;
            minDiff = std::min(minDiff, diff);
            maxDiff = std::max(maxDiff, diff);
        }
        int imba = maxDiff - minDiff;
        imbalances[j] = imba;
        increRegions[j] = increRegion;
        absDiffSums[j] = absDiffSum;
    }
    return;
}


int GraphSearch::getBestSPTree(const vector<int>& imbalances, const vector<int>& increRegions, const vector<int>& absDiffSums)
{
    // int minImbalance = INT_MAX;
    // int minIncreRegion = INT_MAX;
    // int minabsDiffSum = INT_MAX;
    // int best_idx = -1;
    // for(int j = 0; j < imbalances.size(); j++)
    // {
    //     if (imbalances[j] < minImbalance
    //         || (imbalances[j] == minImbalance && increRegions[j] < minIncreRegion)
    //         || (imbalances[j] == minImbalance && increRegions[j] == minIncreRegion && absDiffSums[j] < minabsDiffSum))
    //     {
    //         minImbalance = imbalances[j];
    //         minIncreRegion = increRegions[j];
    //         minabsDiffSum = absDiffSums[j];
    //         best_idx = j;
    //     }
    // }
    // return best_idx;
    vector<float> scores(imbalances.size(), 0);
    for(int i = 0; i < imbalances.size(); i++)
    {
        
        scores[i] = imbalances[i] + increRegions[i] * 1.5 + absDiffSums[i] * 0.2;
    }
    return std::min_element(scores.begin(), scores.end()) - scores.begin();
}

vector<int> GraphSearch::analyzeSPTreeBT(const Graph& h)
{
    vector<vector<int>> allSpanningTrees;
    vector<bool> visited(h.numNodes(), false);
    vector<int> currentTree;
    
    for(int i = 0; i < h.numNodes(); i++)
    {
        DFSUtil(h, i, visited, currentTree, allSpanningTrees);
    }

    // get best spannin tree so the edges span the edge list almost evenly
    vector<int> imbalances;
    vector<int> increRegions;
    vector<int> absDiffSums;
    getStatsSPTree(allSpanningTrees, h.numEdges(), imbalances, increRegions, absDiffSums);
    int best_sptree_idx = getBestSPTree(imbalances, increRegions, absDiffSums);

    // print the best spanning tree
    cout << "Choose spanning tree: ";
    for(auto e_id: allSpanningTrees[best_sptree_idx])
    {
        cout << e_id << ", ";
    }
    cout << endl;
    // return the best spanning tree
    return allSpanningTrees[best_sptree_idx];

}

vector<vector<int>> GraphSearch::analyzeSPTreeSample(const Graph& h)
{
    vector<vector<int>> allSpanningTrees;
    vector<bool> visited(h.numNodes(), false);
    vector<int> currentTree;
    
    for(int i = 0; i < h.numNodes(); i++)
    {
        DFSUtil(h, i, visited, currentTree, allSpanningTrees);
    }

    // get best spannin tree so the edges span the edge list almost evenly
    vector<int> imbalances;
    vector<int> increRegions;
    vector<int> absDiffSums;
    getStatsSPTree(allSpanningTrees, h.numEdges(), imbalances, increRegions, absDiffSums);
    int best_sptree_idx = getBestSPTree(imbalances, increRegions, absDiffSums);

    vector<int>& best_sptree = allSpanningTrees[best_sptree_idx];
    // get the ordering of preprocessing
    // topological sort , where we treat edge as node, and node as edge
    vector<int> degree(h.numNodes(), 0);
    for(int e_id: best_sptree)
    {
        const Edge& e = h.edges()[e_id];
        degree[e.source()]++;
        degree[e.dest()]++;
    }
    vector<int> remaining_edges = best_sptree;
    // find zero degree nodes and find the corresponding edges, then push to a queue
    queue<int> q;
    for(int e_id: best_sptree)
    {
        const Edge& e = h.edges()[e_id];
        if (degree[e.source()] == 1 || degree[e.dest()] == 1)
        {
            q.push(e.index());
            remaining_edges.erase(std::remove(remaining_edges.begin(), remaining_edges.end(), e.index()), remaining_edges.end());
        }
    }
    vector<vector<int>> topo_order;
    while(!q.empty()) //TODO: only one center edge...
    {
        int size = q.size();
        vector<int> cur_level;
        while(size--!=0)
        {
            int e_id = q.front();
            q.pop();
            cur_level.push_back(e_id);
            const Edge& e = h.edges()[e_id];
            degree[e.source()]--;
            degree[e.dest()]--;
            // find the edge in spanning tree whose degree becomes 1 and push to queue
            for(int e_id: remaining_edges)
            {
                const Edge& e = h.edges()[e_id];
                if (degree[e.source()] == 1 || degree[e.dest()] == 1)
                {
                    q.push(e.index());
                    remaining_edges.erase(std::remove(remaining_edges.begin(), remaining_edges.end(), e_id), remaining_edges.end());
                }
            }
        }
        topo_order.push_back(cur_level);
    }
    // make sure that the  last level has only one center edge;
    // if not, we need to move one center edge to the new level
    if (topo_order.back().size() > 1)
    {
        int last_level_size = topo_order.back().size();
        int center_edge = topo_order.back()[last_level_size - 1];
        topo_order.back().pop_back();
        topo_order.push_back(vector<int>{center_edge});
    }

    // print the best spanning tree
    cout << "Choose spanning tree: ";
    for(auto e_id: allSpanningTrees[best_sptree_idx])
    {
        cout << e_id << ", ";
    }
    cout << endl;
    // print the topo order
    cout << "Choose spanning tree: ";
    for(auto level: topo_order)
    {
        cout << "{";
        for(auto e_id: level)
        {
            cout << e_id << ", ";
        }
        cout << "}, ";
    }
    cout << endl;
    // return the best spanning tree
    return topo_order;

}

vector<GraphMatch> GraphSearch::findAllSubgraphs(const Graph &g, const Graph &h, long long int limit)
{
    // If no criteria specified, just use the "dummy" criteria, that accepts everything.
    MatchCriteria criteria;
    return this->findAllSubgraphs(g, h, criteria, limit);
}

vector<GraphMatch> GraphSearch::findAllSubgraphs(const Graph &g, const Graph &h, const MatchCriteria &criteria, long long int limit)
{
    // Store class data structures
    _g = &g;
    _h = &h;
    _criteria = &criteria;

    bool debugOutput = false;

    // Stores the matching subgraphs as list of edge indices
    vector<GraphMatch> results;

    int n = _g->numNodes();
    int m = _g->numEdges();

    // Create lists of nodes that could be mapped to the given nodes
    vector<unordered_set<int>> h2gPossible = this->mapPossibleNodes();

    // Tables for mapping nodes and edges between the two graphs
    // -1 means no match has been assigned yet
    _h2gNodes.clear();
    _h2gNodes.resize(h.numNodes(), -1);
    _g2hNodes.clear();
    _g2hNodes.resize(n, -1);

    // Perform subgraph search, storing results along the way
    int numAssigned = 0;
    this->search(numAssigned, h2gPossible, results);

    return results;
}

vector<unordered_set<int>> GraphSearch::mapPossibleNodes()
{
    vector<unordered_set<int>> possible(_h->numNodes());

    int n = _g->numNodes();
    // Look at each vertex in H
    for (int h_v = 0; h_v < _h->numNodes(); h_v++)
    {
        // cout << "Testing possible nodes for " << h_v << endl;
        //  Check each vertex in G for a possible match
        for (int g_v = 0; g_v < n; g_v++)
        {
            // cout << "Trying " << g_v << endl;
            //  If it passes the criteria, than add it to the list of possible
            if (_criteria->isNodeMatch(*_g, g_v, *_h, h_v))
            {
                possible[h_v].insert(g_v);
            }
        }
    }
    return possible;
}

bool GraphSearch::search(int &numAssigned, vector<unordered_set<int>> &h2gPossible)
{
    // Test if nodes/edges match so far
    if (!matchesSoFar(numAssigned))
        return false;

    // Stop if we've reach the end
    if (numAssigned == _h->numNodes())
        return true;

    // Performs recursive DFS for matches
    int h_v = numAssigned;
    const auto &possible = h2gPossible[h_v];
    for (int g_v : possible)
    {
        if (_g2hNodes[g_v] < 0)
        {
            _h2gNodes[h_v] = g_v;
            _g2hNodes[g_v] = h_v;
            numAssigned++;
            if (search(numAssigned, h2gPossible))
                return true;
            _g2hNodes[g_v] = -1;
            _h2gNodes[h_v] = -1;
            numAssigned--;
        }
    }

    return false;
}

bool GraphSearch::search(int &numAssigned, vector<unordered_set<int>> &h2gPossible, vector<GraphMatch> &results)
{
    // Test if nodes/edges match so far
    if (!matchesSoFar(numAssigned))
        return false;

    // We've found a match if we've reached the end
    if (numAssigned == _h->numNodes())
    {
        // cout << "Found match!" << endl;
        //  Find the matching edges for the assignment
        // vector<int> matchingEdgeIndexes;
        GraphMatch matchingEdges;
        for (const Edge &hEdge : _h->edges())
        {
            int h_u = hEdge.source();
            int h_v = hEdge.dest();
            int g_u = _h2gNodes[h_u];
            int g_v = _h2gNodes[h_v];
            const vector<int> &gEdges = _g->getEdgeIndexes(g_u, g_v);
            for (int g_e : gEdges)
            {
                if (_criteria->isEdgeMatch(*_g, g_e, *_h, hEdge.index()))
                {
                    // matchingEdgeIndexes.push_back(g_e);
                    matchingEdges.addEdge(_g->edges()[g_e], hEdge);
                }
            }
        }
        results.push_back(matchingEdges);
        return true;
    }

    // Performs recursive DFS for matches
    int h_v = numAssigned;
    const auto &possible = h2gPossible[h_v];
    for (int g_v : possible)
    {
        if (_g2hNodes[g_v] < 0)
        {
            _h2gNodes[h_v] = g_v;
            _g2hNodes[g_v] = h_v;
            numAssigned++;
            search(numAssigned, h2gPossible, results);
            _g2hNodes[g_v] = -1;
            _h2gNodes[h_v] = -1;
            numAssigned--;
        }
    }

    return false;
}

bool GraphSearch::matchesSoFar(int numAssigned)
{
    // Check to see if every edge between the currently assigned vertices from
    // the subgraph is a matching edge in our graph
    for (const Edge &edge : _h->edges())
    {
        int h_u = edge.source();
        int h_v = edge.dest();
        if (h_u < numAssigned && h_v < numAssigned)
        {
            int g_u = _h2gNodes[h_u];
            int g_v = _h2gNodes[h_v];

            bool hasEdge = _g->hasEdge(g_u, g_v);
            if (hasEdge)
            {
                // Check to make sure they match the criteria
                const vector<int> edges = _g->getEdgeIndexes(g_u, g_v);
                for (int e : edges)
                {
                    if (_criteria->isEdgeMatch(*_g, e, *_h, edge.index()))
                    {
                        hasEdge = true;
                        break;
                    }
                }
            }
            if (!hasEdge)
            {
                return false;
                break;
            }
        }
    }
    return true;
}

long long int GraphSearch::findOrderedSubgraphs(const Graph &g, const Graph &h, long long int limit, int delta)
{
    // If no criteria specified, just use the "dummy" criteria, that accepts everything.
    MatchCriteria criteria;
    return findOrderedSubgraphs(g, h, criteria, limit, delta);
}

long long int GraphSearch::findOrderedSubgraphs(const Graph &g, const Graph &h, const MatchCriteria &criteria,
                                                long long int limit, int delta)
{
    // Store class data structures
    _g = &g;
    _h = &h;
    _criteria = &criteria;
    _delta = delta;

    bool debugOutput = false;

    // Stores the matching subgraphs as list of edge indices
    long long int results = 0;

    int n = _g->numNodes();
    int m = _g->numEdges();

    // List of all edge indexes
    _allEdges.resize(m);
    for (int i = 0; i < m; i++)
        _allEdges[i] = i;

    // Tables for mapping nodes and edges between the two graphs
    // -1 means no match has been assigned yet
    _h2gNodes.clear();
    _h2gNodes.resize(h.numNodes(), -1);
    _g2hNodes.clear();
    _g2hNodes.resize(n, -1);

    // Keeps track of the number of search edges mapped to a particular
    // node, so we can know if we need to reset its mapping when removing
    // edges from a search trail.
    _numSearchEdgesForNode.clear();
    _numSearchEdgesForNode.resize(n, 0);

    Timer t_findnext;

    // Stores all edges found that match our query.
    // Stack used to backtrack when a particular search ends up a dead-end.
    while (_sg_edgeStack.empty() == false) // Make sure it's empty to start
    {
        _sg_edgeStack.pop();
        //_h_edgeStack.pop();
    }

    // The edge from H we are trying to match in G
    int h_i = 0;
    // The current edge from G we are testing out (yes, should start at -1)
    int g_i = 0;

    // Loop until we can account for all subgraphs matching our edges
    while (true)
    {
        // Time of the last matched edge, for delta comparison
        // time_t lastEdgeTime = 0;
        // if(_sg_edgeStack.empty()==false)
        // lastEdgeTime = g.edges()[_sg_edgeStack.top()].time();

        // Time of current edge, for delta comparison
        time_t curEdgeTime = 0;
        if (g_i < m)
            curEdgeTime = g.edges()[g_i].time();

        // If we've run out of edges, we need to pop the last edge used,
        // and start the search back at the edge after that one
        while (g_i >= m || (_sg_edgeStack.empty() == false && g.edges()[g_i].time() - _firstEdgeTime > delta))
        {
            // If the edge stack is empty, then we have no options left
            // and need to give up.
            if (_sg_edgeStack.empty())
            {
                return results;
            }

            // Pop the stack
            int last_g_i = _sg_edgeStack.top();
            _sg_edgeStack.pop();

            // Get edge object
            const Edge &g_edge = _g->edges()[last_g_i];
            if (debugOutput)
                cout << "Giving up on edge " << last_g_i << ": " << g_edge.source() << ", " << g_edge.dest() << endl;

            // Decrement the number of edges for the nodes in g_i
            _numSearchEdgesForNode[g_edge.source()]--;
            _numSearchEdgesForNode[g_edge.dest()]--;

            // If any of them reach zero, then we need to remove the
            // node mapping for that node, since none of our edges are
            // currently using it (making it free to be re-assigned).
            if (_numSearchEdgesForNode[g_edge.source()] == 0)
            {
                int old_h_u = _g2hNodes[g_edge.source()];
                _h2gNodes[old_h_u] = -1;
                _g2hNodes[g_edge.source()] = -1;
            }
            if (_numSearchEdgesForNode[g_edge.dest()] == 0)
            {
                int old_h_v = _g2hNodes[g_edge.dest()];
                _h2gNodes[old_h_v] = -1;
                _g2hNodes[g_edge.dest()] = -1;
            }

            // Decrement h_i, so that we can find a new one
            h_i--;
            // Make sure we start the search immediately after the failed edge
            g_i = last_g_i + 1;

        }

        // Get query edge
        const Edge &h_edge = _h->edges()[h_i];
        int h_u = h_edge.source();
        int h_v = h_edge.dest();

        // Find matching edge, if possible
        // t_findnext.Start();
        g_i = this->findNextMatch(h_i, g_i);
        // t_findnext.Stop();

        if (g_i < m)
        {
            // Test to see if whole graph is found
            if (h_i + 1 == _h->numEdges())
            {
                // Create GraphMatch object
                // GraphMatch match = convert(_sg_edgeStack, g_i);
                // for(int i=0; i < match.edges().size(); i++)
                //     cout << _g->edges()[match.edges()[i]] << ",";
                // cout << endl;
                // Add new subgraph to the results
                // results.push_back(match);
                results++;
                // Don't increment h_i (or perform mappings), because we want
                // to find if there are other alternative subgraphs for that edge.
                // cout << "Found subgraph #" << results.size() << endl;
                // Test if we've reached our limit, and stop if we have.
                if (results >= limit)
                    return results;
            }
            // Otherwise, add the edge and mappings to the subgraph search
            // and continue on to find next edges.
            else
            {
                // Get matched edge
                const Edge &g_edge = _g->edges()[g_i];
                int g_u = g_edge.source();
                int g_v = g_edge.dest();

                // Set the first edge time, if needed
                if (_sg_edgeStack.empty())
                {
                    _firstEdgeTime = g_edge.time();
                    _firstEdgeid = g_edge.index();
                    _upper_time = _firstEdgeTime + _delta;
                }

                // Map the nodes from each graph
                _h2gNodes[h_u] = g_u;
                _h2gNodes[h_v] = g_v;
                _g2hNodes[g_u] = h_u;
                _g2hNodes[g_v] = h_v;

                // Increment number of search edges for each node in our G edge
                _numSearchEdgesForNode[g_u]++;
                _numSearchEdgesForNode[g_v]++;

                // Add it to the stack
                _sg_edgeStack.push(g_i);
                //_h_edgeStack.push(h_i);

                // Increment to next edge to find
                h_i++;
            }
        }
        // Increment the edge to test
        g_i++;
    }
    return results;
}

long long int GraphSearch::findOrderedSubgraphsMultiThread(const Graph &g, const Graph &h, const MatchCriteria &criteria,
                                                long long int limit, int delta, int num_of_threads, int partition_per_thread)
{
    long long int sum_results = 0;
    // prepare omp
    vector<GraphSearch *> searches;
    vector<long long int *> results;
    for (int i = 0; i < num_of_threads * partition_per_thread; i++)
    {
        searches.push_back(new GraphSearch());
        results.push_back(new long long int());
    }
    vector<int> start_index_vec;
    vector<int> end_index_vec;
    for (int i = 0; i < num_of_threads * partition_per_thread; i++)
    {
        int start_index = int(g.numEdges() * 1.0 * i / (num_of_threads * partition_per_thread));
        int end_index = int(g.numEdges() * 1.0 * (i + 1) / (num_of_threads * partition_per_thread));
        while (start_index > 0 && start_index < g.numEdges() && g.edges()[start_index - 1].source() == g.edges()[start_index - 1].dest())
        {
            start_index++;
        }
        while (end_index > 0 && end_index < g.numEdges() && g.edges()[end_index - 1].source() == g.edges()[end_index - 1].dest())
        {
            end_index++;
        }
        start_index_vec.push_back(start_index);
        end_index_vec.push_back(end_index);
    }
    double avg_time;
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < num_of_threads * partition_per_thread; i++)
    {
        searches[i]->findOrderedSubgraphs(results[i], start_index_vec[i], end_index_vec[i], g, h, criteria, limit, delta);
    }
    for (int i = 0; i < num_of_threads * partition_per_thread; i++)
    {
        sum_results += *(results[i]);
    }
    return sum_results;
}

void GraphSearch::findOrderedSubgraphs(
    long long int *results, int start_edge_idx, int end_edge_idx,
    const Graph &g, const Graph &h, const MatchCriteria &criteria,
    long long int limit, int delta)
{
    // Store class data structures
    _g = &g;
    _h = &h;
    _criteria = &criteria;
    _delta = delta;

    bool debugOutput = false;

    // Stores the matching subgraphs as list of edge indices

    int n = _g->numNodes();
    int m = _g->numEdges();

    // List of all edge indexes
    _allEdges.resize(m);
    for (int i = 0; i < m; i++)
        _allEdges[i] = i;

    // Tables for mapping nodes and edges between the two graphs
    // -1 means no match has been assigned yet
    _h2gNodes.clear();
    _h2gNodes.resize(h.numNodes(), -1);
    _g2hNodes.clear();
    _g2hNodes.resize(n, -1);

    // Keeps track of the number of search edges mapped to a particular
    // node, so we can know if we need to reset its mapping when removing
    // edges from a search trail.
    _numSearchEdgesForNode.clear();
    _numSearchEdgesForNode.resize(n, 0);


    Timer t_findnext;

    // Stores all edges found that match our query.
    // Stack used to backtrack when a particular search ends up a dead-end.
    while (_sg_edgeStack.empty() == false) // Make sure it's empty to start
    {
        _sg_edgeStack.pop();
        //_h_edgeStack.pop();
    }

    // The edge from H we are trying to match in G
    int h_i = 0;
    // The current edge from G we are testing out (yes, should start at -1)
    int g_i = start_edge_idx;
    // Loop until we can account for all subgraphs matching our edges
    while (true)
    {
        // Time of the last matched edge, for delta comparison
        // time_t lastEdgeTime = 0;
        // if(_sg_edgeStack.empty()==false)
        // lastEdgeTime = g.edges()[_sg_edgeStack.top()].time();

        // Time of current edge, for delta comparison
        time_t curEdgeTime = 0;
        if (g_i < m)
            curEdgeTime = g.edges()[g_i].time();

        // If we've run out of edges, we need to pop the last edge used,
        // and start the search back at the edge after that one
        while (g_i >= m || _sg_edgeStack.empty() || (_sg_edgeStack.empty() == false && g.edges()[g_i].time() - _firstEdgeTime > delta))
        {
            // If the edge stack is empty, then we have no options left
            // and need to give up.
            if (_sg_edgeStack.empty())
            {
                if (g_i >= end_edge_idx)
                    return;
                else
                    break;
            }

            // Pop the stack
            int last_g_i = _sg_edgeStack.top();
            _sg_edgeStack.pop();

            // Get edge object
            const Edge &g_edge = _g->edges()[last_g_i];
            if (debugOutput)
                cout << "Giving up on edge " << last_g_i << ": " << g_edge.source() << ", " << g_edge.dest() << endl;

            // Decrement the number of edges for the nodes in g_i
            _numSearchEdgesForNode[g_edge.source()]--;
            _numSearchEdgesForNode[g_edge.dest()]--;

            // If any of them reach zero, then we need to remove the
            // node mapping for that node, since none of our edges are
            // currently using it (making it free to be re-assigned).
            if (_numSearchEdgesForNode[g_edge.source()] == 0)
            {
                int old_h_u = _g2hNodes[g_edge.source()];
                _h2gNodes[old_h_u] = -1;
                _g2hNodes[g_edge.source()] = -1;
            }
            if (_numSearchEdgesForNode[g_edge.dest()] == 0)
            {
                int old_h_v = _g2hNodes[g_edge.dest()];
                _h2gNodes[old_h_v] = -1;
                _g2hNodes[g_edge.dest()] = -1;
            }

            // Decrement h_i, so that we can find a new one
            h_i--;
            // Make sure we start the search immediately after the failed edge
            g_i = last_g_i + 1;
        }

        // Get query edge
        const Edge &h_edge = _h->edges()[h_i];
        int h_u = h_edge.source();
        int h_v = h_edge.dest();

        // Find matching edge, if possible
        // t_findnext.Start();
        g_i = this->findNextMatch(h_i, g_i);
        // t_findnext.Stop();

        if (g_i < m)
        {
            // Test to see if whole graph is found
            if (h_i + 1 == _h->numEdges())
            {
                // Convert stack to vector
                // vector<int> edges = convert(_sg_edgeStack);
                // Add the last edge to the list
                // edges.push_back(g_i);
                // Add new subgraph to the results
                // results.push_back(edges);

                // Create GraphMatch object
                // GraphMatch match = convert(_sg_edgeStack, g_i);
                // for(int i=0; i < match.edges().size(); i++)
                //     cout << match.edges()[i] << ",";
                // cout << endl;
                // Add new subgraph to the results
                // results.push_back(match);
                (*results)++;
                // Don't increment h_i (or perform mappings), because we want
                // to find if there are other alternative subgraphs for that edge.
                // cout << "Found subgraph #" << results.size() << endl;
                // Test if we've reached our limit, and stop if we have.
                if ((*results) >= limit)
                    return;
            }
            // Otherwise, add the edge and mappings to the subgraph search
            // and continue on to find next edges.
            else
            {
                // Get matched edge
                const Edge &g_edge = _g->edges()[g_i];
                int g_u = g_edge.source();
                int g_v = g_edge.dest();

                // Set the first edge time, if needed
                if (_sg_edgeStack.empty())
                {
                    _firstEdgeTime = g_edge.time();
                    _firstEdgeid = g_edge.index();
                    _upper_time = _firstEdgeTime + _delta;
                }

                // Map the nodes from each graph
                _h2gNodes[h_u] = g_u;
                _h2gNodes[h_v] = g_v;
                _g2hNodes[g_u] = h_u;
                _g2hNodes[g_v] = h_v;

                // Increment number of search edges for each node in our G edge
                _numSearchEdgesForNode[g_u]++;
                _numSearchEdgesForNode[g_v]++;

                // Add it to the stack
                _sg_edgeStack.push(g_i);
                //_h_edgeStack.push(h_i);

                // Increment to next edge to find
                h_i++;
            }
        }
        // Increment the edge to test
        g_i++;
    }
    return;
}

int GraphSearch::findNextMatch(int h_i, int g_i)
{
    bool debugOutput = false;

    // Get query edge
    const Edge &h_edge = _h->edges()[h_i];
    int h_u = h_edge.source();
    int h_v = h_edge.dest();

    // Default is to search over all edges starting at g_i
    const vector<int> *searchEdges = &_allEdges;

    // Look to see if nodes are already mapped, and just use those
    // node edges, if so. (Much faster!)
    if (_h2gNodes[h_u] >= 0 && _h2gNodes[h_v] >= 0)
    {
        const vector<int> &uEdges = _g->nodes()[_h2gNodes[h_u]].outEdges();
        const vector<int> &vEdges = _g->nodes()[_h2gNodes[h_v]].inEdges();
        if (uEdges.size() < vEdges.size())
            searchEdges = &uEdges;
        else
            searchEdges = &vEdges;
    }
    else if (_h2gNodes[h_u] >= 0)
    {
        searchEdges = &_g->nodes()[_h2gNodes[h_u]].outEdges();
    }
    else if (_h2gNodes[h_v] >= 0)
    {
        searchEdges = &_g->nodes()[_h2gNodes[h_v]].inEdges();
    }

    // Find starting place in the list
    int start = findStart(g_i, *searchEdges);

    // If no starting place can be found, just return that it's not possible
    if (start == searchEdges->size())
    {
        return _g->numEdges();
    }
    /*if(debugOutput)
    {
        cout << "Looking for an edge >= " << g_i << endl;
        cout << "Found this one: " <<
    }*/

    // Perform search
    int fnm = findNextMatch(h_i, *searchEdges, start);

    return fnm;
}

int GraphSearch::findStart(int g_i, const std::vector<int> &edgeIndexes)
{
    // If it's the original edges, just return g_i
    if (edgeIndexes.size() > g_i && edgeIndexes[g_i] == g_i)
        return g_i;

    // Test if any edge will work
    if (edgeIndexes.empty())
        return edgeIndexes.size();
    if (edgeIndexes.back() < g_i)
        return edgeIndexes.size();
    if (edgeIndexes.front() >= g_i)
        return 0;

    /*for(int i=0; i<edgeIndexes.size(); i++)
    {
        if(edgeIndexes[i] >= g_i)
            return i;
    }
    // If nothing was found, return the size of our list
    return edgeIndexes.size();*/

    // Otherwise, perform binary search
    int left = 0, right = edgeIndexes.size() - 1;
    while (true)
    {
        if (right <= left)
        {
            return left;
        }
        int i = (right + left) / 2;
        int ei = edgeIndexes[i];
        if (ei == g_i)
            return i;
        if (ei >= g_i && i == left)
            return i;
        if (ei < g_i)
            left = i + 1;
        else
        {
            if (edgeIndexes[i - 1] < g_i)
            {
                return i;
            }
            right = i - 1;
        }
    }
}

int GraphSearch::findNextMatch(int h_i, const std::vector<int> &edgesToSearch, int startIndex)
{
    bool debugOutput = false;

    // Get query edge
    const Edge &h_edge = _h->edges()[h_i];
    int h_u = h_edge.source();
    int h_v = h_edge.dest();

    // Check the time against the previous matched edge, if any exist
    bool checkTime = _sg_edgeStack.empty() == false;

    // Loop over all the edges to search
    for (int i = startIndex; i < edgesToSearch.size(); i++)
    {
        // Get the index of our edge in G
        int g_i = edgesToSearch[i];

        // Get original edge
        const Edge &g_edge = _g->edges()[g_i];
        int g_u = g_edge.source();
        int g_v = g_edge.dest();

        // If we've gone past our delta, stop the search
        // if(checkTime && g_edge.time() - _firstEdgeTime > _delta)
        if (checkTime && g_edge.time() > _upper_time)
        {
            return _g->numEdges();

            if (debugOutput)
            {
                cout << "Trying edge " << g_i << ": " << g_u << ", " << g_v << "    ";
                cout << "Need to match edge " << h_i << ": " << h_u << ", " << h_v << endl;
                cout << "   g[" << g_u << "]=" << _g2hNodes[g_u] << " g[" << g_v << "]=" << _g2hNodes[g_v] << endl;
                cout << "   h[" << h_u << "]=" << _h2gNodes[h_u] << " h[" << h_v << "]=" << _h2gNodes[h_v] << endl;
            }
        }
        // Make sure if the edge is a self-loop or not
        if ((h_u == h_v && g_u == g_v) || (h_u != h_v && g_u != g_v))
        {
            // Test if source nodes match, or both are unassigned
            if (_h2gNodes[h_u] == g_u || (_h2gNodes[h_u] < 0 && _g2hNodes[g_u] < 0))
            {
                // Test if destination nodes match, or both are unassigned
                if (_h2gNodes[h_v] == g_v || (_h2gNodes[h_v] < 0 && _g2hNodes[g_v] < 0))
                {
                    // Test if metadata criteria is a match
                    if (_criteria->isEdgeMatch(*_g, g_i, *_h, h_i))
                    {
                        if (debugOutput)
                            cout << "Edge " << g_i << ": " << g_u << ", " << g_v << " is a match" << endl;
                        return g_i;
                    }
                }
            }
        }
    }
    // If no match found, return the number of edges
    return _g->numEdges();
}

vector<int> GraphSearch::convert(stack<int> s)
{
    vector<int> v(s.size());
    for (int i = v.size() - 1; i >= 0; i--)
    {
        v[i] = s.top();
        s.pop();
    }
    return v;
}

GraphMatch GraphSearch::convert(const std::stack<int> &s, int g_lastEdge)
{
    GraphMatch gm;
    vector<int> gEdges = convert(s);
    gEdges.push_back(g_lastEdge);
    for (int h_i = 0; h_i < gEdges.size(); h_i++)
    {
        int g_i = gEdges[h_i];
        gm.addEdge(_g->edges()[g_i], _h->edges()[h_i]);
    }
    return gm;
}

long long int GraphSearch::findOrderedSubgraphsSpanningTreeMultiThread(
    const Graph &g, const Graph &h, const MatchCriteria &criteria,
    const vector<int> &spanning_tree, long long int limit, int delta,
    int num_of_threads, int partition_per_thread)
{
    // store the lower and upperbound respect to spanning tree edges of the extra edges that are not in the spanning tree
    // SPTreeRangeEdges: (SPtree edge_id 1, SPtree edge_id 2, vector<int> h_edge_ids of extra edges)
    // the edges in the vector<int> have the index range of (SPtree edge_id 1, SPtree edge_id 2)
    // edges with same range are counted by n-pointer technique (linear time complexity)
    vector<int> sorted_spanning_tree = spanning_tree;
    sort(sorted_spanning_tree.begin(), sorted_spanning_tree.end());
    map<pair<int, int>, vector<int>> sp_tree_range_edges = analyze_exatra_edges(h, sorted_spanning_tree);
    cout << "sp_tree_range_edges: " << endl;
    for(auto &[k, vs]: sp_tree_range_edges)
    {
        cout << "[ " << k.first << "," << k.second << " ] : ";
        for(auto v: vs)
        {
            cout << v << ", ";
        }
        cout << endl;
    }
    // split the spanning_tree into non-deceasing sublist
    // For example, given matching order as {1, 2, 3, 0}, return value is {{1, 2, 3}, {0}}.
    vector<vector<int>> matching_regions_order = getMatchingRegions(spanning_tree);
    cout << "matching_regions_order: " << endl;
    for(auto &v: matching_regions_order)
    {
        cout << "{";
        for(auto vv: v)
        {
            cout << vv << ", ";
        }
        cout << "},";
    }
    cout << endl;


    long long int sum_results = 0;
    // prepare omp
    vector<GraphSearch *> searches;
    vector<long long int > results;
    for (int i = 0; i < num_of_threads * partition_per_thread; i++)
    {
        searches.push_back(new GraphSearch());
        results.push_back(0);
    }
    vector<int> start_index_vec;
    vector<int> end_index_vec;
    for (int i = 0; i < num_of_threads * partition_per_thread; i++)
    {
        int start_index = int(g.numEdges() * 1.0 * i / (num_of_threads * partition_per_thread));
        int end_index = int(g.numEdges() * 1.0 * (i + 1) / (num_of_threads * partition_per_thread));
        while (start_index > 0 && start_index < g.numEdges() && g.edges()[start_index - 1].source() == g.edges()[start_index - 1].dest())
        {
            start_index++;
        }
        while (end_index > 0 && end_index < g.numEdges() && g.edges()[end_index - 1].source() == g.edges()[end_index - 1].dest())
        {
            end_index++;
        }
        start_index_vec.push_back(start_index);
        end_index_vec.push_back(end_index);
    }
    double avg_time;
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < num_of_threads * partition_per_thread; i++)
    {
        results[i] = searches[i]->findOrderedSubgraphsSpanningTreeInner(
            start_index_vec[i], end_index_vec[i], g, h, criteria, 
            sorted_spanning_tree, sp_tree_range_edges, matching_regions_order,
            limit, delta);
    }
    for (int i = 0; i < num_of_threads * partition_per_thread; i++)
    {
        // cout << i << " results[i]: " << results[i] << endl;
        sum_results += (results[i]);
    }
    for (int i = 0; i < num_of_threads * partition_per_thread; i++)
    {
        delete searches[i];
    }
    return sum_results;

}

long long int GraphSearch::findOrderedSubgraphsSpanningTree(
    const Graph &g, const Graph &h, const MatchCriteria &criteria,
    const vector<int> &spanning_tree, long long int limit, int delta)
{
    // store the lower and upperbound respect to spanning tree edges of the extra edges that are not in the spanning tree
    // SPTreeRangeEdges: (SPtree edge_id 1, SPtree edge_id 2, vector<int> h_edge_ids of extra edges)
    // the edges in the vector<int> have the index range of (SPtree edge_id 1, SPtree edge_id 2)
    // edges with same range are counted by n-pointer technique (linear time complexity)
    vector<int> sorted_spanning_tree = spanning_tree;
    sort(sorted_spanning_tree.begin(), sorted_spanning_tree.end());
    map<pair<int, int>, vector<int>> sp_tree_range_edges = analyze_exatra_edges(h, sorted_spanning_tree);
    cout << "sp_tree_range_edges: " << endl;
    for(auto &[k, vs]: sp_tree_range_edges)
    {
        cout << "[ " << k.first << "," << k.second << " ] : ";
        for(auto v: vs)
        {
            cout << v << ", ";
        }
        cout << endl;
    }
    // split the spanning_tree into non-deceasing sublist
    // For example, given matching order as {1, 2, 3, 0}, return value is {{1, 2, 3}, {0}}.
    vector<vector<int>> matching_regions_order = getMatchingRegions(spanning_tree);
    cout << "matching_regions_order: " << endl;
    for(auto &v: matching_regions_order)
    {
        cout << "{";
        for(auto vv: v)
        {
            cout << vv << ", ";
        }
        cout << "},";
    }
    cout << endl;
    return findOrderedSubgraphsSpanningTreeInner(
        0, g.numEdges(), g, h, criteria, sorted_spanning_tree, sp_tree_range_edges, matching_regions_order, limit, delta);
}

long long int GraphSearch::findOrderedSubgraphsSpanningTreeInner(
    int start_edge_idx, int end_edge_idx,
    const Graph &g, const Graph &h, const MatchCriteria &criteria,
    const vector<int> &spanning_tree, const map<pair<int, int>, vector<int>>& sp_tree_range_edges,
    const vector<vector<int>>& matching_regions_order,
    long long int limit, int delta)
{

    /*
    spanning_tree are the list of edges in motif h that selected to be the spanning tree
    */
    // Store class data structures
    _g = &g;
    _h = &h;
    _criteria = &criteria;
    _delta = delta;
    _matching_regions_order = matching_regions_order;

    bool debugOutput = false;

    // Stores the matching subgraphs as list of edge indices
    long long int results = 0;

    int n = _g->numNodes();
    int m = _g->numEdges();

    // List of all edge indexes
    _allEdges.resize(m);
    for (int i = 0; i < m; i++)
        _allEdges[i] = i;

    // Tables for mapping nodes and edges between the two graphs
    // -1 means no match has been assigned yet
    _h2gNodes.clear();
    _h2gNodes.resize(_h->numNodes(), -1);
    _g2hNodes.clear();
    _g2hNodes.resize(n, -1);

    // Keeps track of the number of search edges mapped to a particular
    // node, so we can know if we need to reset its mapping when removing
    // edges from a search trail.
    _numSearchEdgesForNode.clear();
    _numSearchEdgesForNode.resize(n, 0);

    // the edge of the motif that we are trying to map to so far
    vector<Edge> dfs_edges(_h->numEdges(), Edge(-1, -1, -1, -1));


    // Stores all edges found that match our query.
    // Stack used to backtrack when a particular search ends up a dead-end.
    while (_sg_edgeStack.empty() == false) // Make sure it's empty to start
    {
        _sg_edgeStack.pop();
        //_h_edgeStack.pop();
    }

    
    // vars to keep track of which motif edge we are currently trying to map to, init to 0
    _region_i = 0;
    _matching_h_i = 0;

    // The edge from H we are trying to match in G
    int h_i = _matching_regions_order[_region_i][_matching_h_i];
    // The current edge from G we are testing out (yes, should start at -1)
    int g_i = start_edge_idx;

    // upper bounds to see if we've run out of edges in the current level of DFS search tree
    _upper_id = _g->numEdges();
    _upper_time = INT_MAX;

    // Loop until we can account for all subgraphs matching our edges
    while (true)
    {

        // If we've run out of edges, we need to pop the last edge used,
        // and start the search back at the edge after that one
        while ((g_i >= _upper_id) || _sg_edgeStack.empty() || (_sg_edgeStack.empty() == false && _g->edges()[g_i].time() > _upper_time))
        {
            // If the edge stack is empty, then we have no options left
            // and need to give up.
            if (_sg_edgeStack.empty())
            {
                if (g_i >= end_edge_idx)
                    return results;
                else
                    break;
            }

            // Pop the stack
            int last_g_i = _sg_edgeStack.top();
            _sg_edgeStack.pop();

            // Get edge object
            const Edge &g_edge = _g->edges()[last_g_i];
            if (debugOutput)
                cout << "Giving up on edge " << last_g_i << ": " << g_edge.source() << ", " << g_edge.dest() << endl;

            // Decrement the number of edges for the nodes in g_i
            _numSearchEdgesForNode[g_edge.source()]--;
            _numSearchEdgesForNode[g_edge.dest()]--;

            // If any of them reach zero, then we need to remove the
            // node mapping for that node, since none of our edges are
            // currently using it (making it free to be re-assigned).
            if (_numSearchEdgesForNode[g_edge.source()] == 0)
            {
                int old_h_u = _g2hNodes[g_edge.source()];
                _h2gNodes[old_h_u] = -1;
                _g2hNodes[g_edge.source()] = -1;
            }
            if (_numSearchEdgesForNode[g_edge.dest()] == 0)
            {
                int old_h_v = _g2hNodes[g_edge.dest()];
                _h2gNodes[old_h_v] = -1;
                _g2hNodes[g_edge.dest()] = -1;
            }

            // Decrement region_i and search_h_i, so that we can find a new one
            decreNextEdge();
            // Make sure we start the search immediately after the failed edge
            g_i = last_g_i + 1;
        }

        // Get query edge
        const Edge &h_edge = _h->edges()[_matching_regions_order[_region_i][_matching_h_i]];
        int h_u = h_edge.source();
        int h_v = h_edge.dest();

        g_i = this->findNextMatch(_matching_regions_order[_region_i][_matching_h_i], g_i);
        // cout << "g_i: " << g_i << ", h_i: [" << _region_i << "][" << _matching_h_i << "]" << endl;

        if (g_i < _upper_id)
        {

            // Test to see if whole graph is found
            if ((_region_i == _matching_regions_order.size() - 1) && (_matching_h_i == _matching_regions_order[_region_i].size() - 1))
            {

                // GraphMatch match = convert_flip(_sg_edgeStack, g_i);
                // for(int i=0; i < match.edges().size(); i++)
                //     cout << match.edges()[i] << ",";
                // cout << endl;
                // Add new subgraph to the results
                // results.push_back(match);
                const Edge &g_edge = _g->edges()[g_i];
                int g_u = g_edge.source();
                int g_v = g_edge.dest();

                dfs_edges[h_edge.index()] = g_edge;

                // Map the nodes from each graph, adding the last edge
                vector<int> h2gNodes = _h2gNodes; // a copy of _h2gNodes
                h2gNodes[h_u] = g_u;
                h2gNodes[h_v] = g_v;

                // after the first few levels (spanning tree chosen) of DFS, we can start to count the number of subgraphs that extended
                // from the the selected spanning tree by using pointer technique
                results += deriveMotifCounts(dfs_edges, spanning_tree, sp_tree_range_edges, h2gNodes);

                g_i++;

                // Don't increment h_i (or perform mappings), because we want
                // to find if there are other alternative subgraphs for that edge.
                // cout << "Found subgraph #" << results.size() << endl;
                // Test if we've reached our limit, and stop if we have.
                if (results >= limit)
                    return results;
            }
            // Otherwise, add the edge and mappings to the subgraph search
            // and continue on to find next edges.
            else
            {
                // Get matched edge
                const Edge &g_edge = _g->edges()[g_i];
                int g_u = g_edge.source();
                int g_v = g_edge.dest();

                // Set the regionFirstEdge and regionLastEdge id and time, if needed
                if (_matching_h_i == 0)
                {
                    MatchRegion *cur_region = new MatchRegion();
                    cur_region->regionFirstEdgeTime = g_edge.time();
                    cur_region->regionFirstEdgeid = g_edge.index();
                    _matching_regions_info.push_back(cur_region);
                }
                if (_matching_h_i == _matching_regions_order[_region_i].size() - 1)
                {
                    MatchRegion *cur_region = _matching_regions_info[_matching_regions_info.size() - 1];
                    cur_region->regionLastEdgeTime = g_edge.time();
                    cur_region->regionLastEdgeid = g_edge.index();
                }

                // Map the nodes from each graph
                _h2gNodes[h_u] = g_u;
                _h2gNodes[h_v] = g_v;
                _g2hNodes[g_u] = h_u;
                _g2hNodes[g_v] = h_v;

                // Increment number of search edges for each node in our G edge
                _numSearchEdgesForNode[g_u]++;
                _numSearchEdgesForNode[g_v]++;

                // Add it to the stack
                _sg_edgeStack.push(g_i);
                // Add it to the dfs_edges
                dfs_edges[h_edge.index()] = g_edge;

                // Increment to next edge to find
                g_i = increNextEdge(g_i);
            }
        }
        // Increment the edge to test
        // g_i++;
    }
    return results;
}

/** Get vector of matching order of regions given the matching order. In each region, the
 * matching order is consecutively increasing. For example, given matching order as {1, 2, 3, 0},
 * return value is {{1, 2, 3}, {0}}.
 */
vector<vector<int>> GraphSearch::getMatchingRegions(const vector<int> &matching_order)
{
    vector<vector<int>> regions;
    vector<int> region;
    if (matching_order.size() == 0)
        return regions;
    for (int i = 0; i < matching_order.size(); i++)
    {
        if ((region.size() != 0) && (matching_order[i] != matching_order[i - 1] + 1))
        {
            regions.push_back(region);
            region.clear();
        }
        region.push_back(matching_order[i]);
    }
    regions.push_back(region);
    return regions;
}

int GraphSearch::increNextEdge(int g_i)
{
    if (_matching_h_i + 1 < _matching_regions_order[_region_i].size())
    {
        // only need to update once for each region
        if (_matching_h_i == 0)
        {
            updateRegionUpperTimeOnCur();
        }
        _matching_h_i++;
        // next edge is in the same region, increment the g_i to start
        return g_i + 1;
    }
    else
    {
        _region_i++;
        _matching_h_i = 0;
        if (_region_i < _matching_regions_order.size())
        {
            return updateRegionUpperBoundOnPrev();
        }
        else
        {
            return _g->edges().size();
        }
    }
}

void GraphSearch::decreNextEdge()
{
    if (_matching_h_i > 0)
    {
        _matching_h_i--;
    }
    else
    {
        _region_i--;
        if (_region_i >= 0)
        {
            _matching_h_i = _matching_regions_order[_region_i].size() - 1;
            updateRegionUpperBoundOnPrev();
        }
    }
    if (_matching_h_i == 0)
    {
        // clear current level's matching region info
        delete(_matching_regions_info.back());
        _matching_regions_info.pop_back();
    }
}

/** Update the upper bound for the current region _region_i given the previous region info
 * So that when g_i < upper_id or g.edges()[g_i].time() <= upper_time,
 * we know that we have run out of input edges g_i matching h_i
 * Return the start idx of the next level edge to search
 */
int GraphSearch::updateRegionUpperBoundOnPrev()
{
    int lower_id = 0;
    // reset
    _upper_id = _g->numEdges();
    _upper_time = INT_MAX;
    if (_region_i == 0)
    {
        updateRegionUpperTimeOnCur();
    }
    for (int i = 0; i < _region_i; i++)
    {
        vector<int> &prev_region = _matching_regions_order[i];
        if (_matching_regions_order[_region_i][0] > prev_region[prev_region.size() - 1])
        {
            _upper_time = min(_upper_time, _matching_regions_info[i]->regionFirstEdgeTime + _delta);
            lower_id = max(lower_id, _matching_regions_info[i]->regionLastEdgeid + 1);
        }
        // _matching_regions_order[_region_i][-1] < prev_region[0]
        else
        {
            _upper_id = min(_upper_id, _matching_regions_info[i]->regionFirstEdgeid);
            lower_id = max(lower_id, findRootEdgeStartLowerTime(_matching_regions_info[i]->regionLastEdgeTime - _delta));
        }
    }
    return lower_id;
}

/** Update the upper time for current region next _matching_h_i based on this region's
 * first edge time.
 * For example, region = {edge1, edge2, edge3}. Edge3's time <= edge1's time + delta
 */
void GraphSearch::updateRegionUpperTimeOnCur()
{
    _upper_time = _matching_regions_info[_region_i]->regionFirstEdgeTime + _delta;
}

int GraphSearch::findRootEdgeStartLowerTime(time_t lower_time)
{
    const vector<int> &edgeIndexes = _allEdges;
    // perform binary search
    int left = 0, right = edgeIndexes.size() - 1;
    while (true)
    {
        if (right <= left)
        {
            return left;
        }
        int i = (right + left) / 2;
        int ei = edgeIndexes[i];
        if (_g->edges()[ei].time() == lower_time)
            right = i;
        if (_g->edges()[ei].time() >= lower_time && i == left)
            return i;
        if (_g->edges()[ei].time() < lower_time)
            left = i + 1;
        else
        {
            if (_g->edges()[edgeIndexes[i - 1]].time() < lower_time)
            {
                return i;
            }
            right = i - 1;
        }
    }
}

/*Counting increasing 3-stars.*/
long long int GraphSearch::count_stars(int delta, int id_begin, int id_end)
{
    long long int star_count = 0;

    for(int i = id_begin; i < id_end; i++) /*Consider a temporal edge u->v*/
    {
        Edge e = _g->edges()[i];
        int e_src = e.source();
        int e_dst = e.dest();
        if(e_src == e_dst){
            continue;
        }
        const vector<int>& e_dst_out_edges = _g->nodes()[e_dst].outEdges(); /*A vector consisting of the id's of the edges directed away from v*/
        
        // e_dst_out_edges.index() > e.index()
        vector<int>::const_iterator e_dst_out_edges_start_it = upper_bound(           /*The starting point of edges directed away from v, with timestamp at least t(u->v)*/
            e_dst_out_edges.begin(), e_dst_out_edges.end(), e.index()
        );
        // e_dst_out_edges.time() <= e.time() + _delta
        vector<int>::const_iterator e_dst_out_edges_end_it = upper_bound(           /*The ending point of edges directed away from v, with timestamp at most t(u->v) + delta*/
            e_dst_out_edges.begin(), e_dst_out_edges.end(), e.time() + _delta,
            [&](const time_t &a, const int &b) { return a < _g->edges()[b].time(); }
        );

        long long int number_of_candidate_edges = 0;
        long long int number_of_pairs_on_same_edges = 0;
        map<int, int>edge_count_per_vertex;

        /*Creates a map, that for each dst vertex w of the outgoing edges v->w, stores the number of edges from v to w.*/
        for(vector<int>::const_iterator edge_iterator = e_dst_out_edges_start_it; edge_iterator != e_dst_out_edges_end_it; edge_iterator++)
        {
            if(_g->edges()[*edge_iterator].dest() == e_src || (_g->edges()[*edge_iterator].dest() == _g->edges()[*edge_iterator].source()))
                continue;
            else{
                edge_count_per_vertex[_g->edges()[*edge_iterator].dest()] += 1;
                number_of_candidate_edges += 1;
            }
        }

        /*Total number of pairs obtained from edges between same pair of nodes (v, w).*/
        for(auto itr = edge_count_per_vertex.begin(); itr != edge_count_per_vertex.end(); itr++)
            number_of_pairs_on_same_edges += (itr->second)*((itr->second) - 1)/2;

        /*Since you need to form pairs out of edges within the same list, no need to execute 'merge' subroutine. */
        long long int this_star_count = (((number_of_candidate_edges)*(number_of_candidate_edges - 1)/2) - number_of_pairs_on_same_edges);
        star_count += this_star_count;
    }
    return star_count;
}


long long int GraphSearch::count_stars_multi_thread(
    int delta, int num_of_threads, int partition_per_thread)
{
    long long int sum_star_count = 0;
    vector<long long int> star_counts(num_of_threads * partition_per_thread, 0);

    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < num_of_threads * partition_per_thread; i++)
    {
        int id_begin = _g->numEdges() / (num_of_threads * partition_per_thread) * i;
        int id_end =  min(_g->numEdges() / (num_of_threads * partition_per_thread) * (i+1), _g->numEdges());
        star_counts[i] = count_stars(delta, id_begin, id_end);
    }

    for (int i = 0; i < num_of_threads * partition_per_thread; i++)
    {
        sum_star_count += star_counts[i];
    }
    
    return sum_star_count;
}


bool GraphSearch::sample_directed_3path(
    int edge_id, vector<vector<vector<int>::const_iterator>>& sampling_neigh_edges_it,
    vector<Edge>& sampled_edges, vector<long long int>& failure)
{
    unsigned int seed = time(NULL) ^ omp_get_thread_num();
    sampled_edges.clear();
    // sample edge (u,v) with probability (in_d_u) * (out_d_v)
    // int edge_id = sample_weight(gen);
    // if(edge_id == 453617 || edge_id == 577917)
    //     cout << "edge_id:" << edge_id << endl;
    Edge u_v_edge = _g->edges()[edge_id];
    int u = u_v_edge.source();
    int v = u_v_edge.dest();
    // if (u == v) return 0; // already checked by preprocess_sampling_weights()
    sampled_edges.push_back(u_v_edge);
    // sample u'->u edge from all in edges of u within delta timestamp constraints
    int num_u_in_edges = distance(sampling_neigh_edges_it[edge_id][0], sampling_neigh_edges_it[edge_id][1]);
    if (num_u_in_edges == 0)
        return false;
    int u_in_edges_randidx = rand_r(&seed) % num_u_in_edges;
    int u_in_edge_id = *(sampling_neigh_edges_it[edge_id][0] + u_in_edges_randidx);
    // check if u'->u edge meets requirement
    Edge u_in_edge = _g->edges()[u_in_edge_id];
    int u_prime = u_in_edge.source();
    // u' and v must be different, u'->u edge time must < u->v time
    if ((u_prime == v) || (u_prime == u))
    {
        failure[0]++;
        return false;
    }
    sampled_edges.push_back(u_in_edge);
    // sample v->v' edge from all out edges of v within delta timestamp constraints
    int num_v_out_edges = distance(sampling_neigh_edges_it[edge_id][2], sampling_neigh_edges_it[edge_id][3]);
    if (num_v_out_edges == 0)
        return false;
    int v_out_edge_randidx = rand_r(&seed) % num_v_out_edges;
    int v_out_edge_id = *(sampling_neigh_edges_it[edge_id][2] + v_out_edge_randidx);
    // v' and u must be different
    Edge v_out_edge = _g->edges()[v_out_edge_id];
    int v_prime = v_out_edge.dest();
    if ((v_prime == u) || (v_prime == v) || (v_prime == u_prime))
    {
        failure[0]++;
        return false;
    }
    // delta constraints
    if (v_out_edge.time() > u_in_edge.time() + _delta)
    {
        failure[1]++;
        return false;
    }
    sampled_edges.push_back(v_out_edge);
    return true;
}

vector<long long int> GraphSearch::check_motif(vector<Edge> &sampled_edges)
{
    vector<long long int> motifs_cnts(8, 0);
    Edge u_v_edge = sampled_edges[0];
    Edge u_in_edge = sampled_edges[1];
    Edge v_out_edge = sampled_edges[2];
    int u = u_v_edge.source();
    int v = u_v_edge.dest();
    int u_prime = u_in_edge.source();
    int v_prime = v_out_edge.dest();

    // M1: 3-path
    motifs_cnts[1] = 1;
    int v_prime_u_edges_num = 0;

    vector<int>::iterator v_prime_u_edges_left_it;
    vector<int>::iterator v_prime_u_edges_right_it;
    vector<int>::iterator v_prime_u_prime_edges_left_it;
    vector<int>::iterator v_prime_u_prime_edges_right_it;
    vector<int>::iterator u_prime_v_edges_left_it;
    vector<int>::iterator u_prime_v_edges_right_it;
    vector<int>::iterator u_u_prime_edges_left_it;
    vector<int>::iterator u_u_prime_edges_right_it;
    bool v_prime_u_exist = false;
    bool v_prime_u_prime_exist = false;
    bool u_prime_v_exist = false;
    bool u_u_prime_exist = false;

    // M2: tailed triangle v'->u
    if (_g->nodeEdges().find(v_prime) != _g->nodeEdges().end())
    {
        if(_g->nodeEdges()[v_prime].find(u) != _g->nodeEdges()[v_prime].end())
        {
            v_prime_u_exist = true;
            vector<int> &v_prime_u_edges_ids = _g->nodeEdges()[v_prime][u];
            // find idx > v_out_edge.index()
            v_prime_u_edges_left_it = upper_bound(
                v_prime_u_edges_ids.begin(), v_prime_u_edges_ids.end(), v_out_edge.index()
            );
            // find timestamp <= u_in_edge.time() + delta
            v_prime_u_edges_right_it = upper_bound(
                v_prime_u_edges_ids.begin(), v_prime_u_edges_ids.end(), u_in_edge.time() + _delta,
                [&](const time_t &a, const int &b) { return a < _g->edges()[b].time(); }
            );
            v_prime_u_edges_num = distance(v_prime_u_edges_left_it, v_prime_u_edges_right_it);
            if (v_prime_u_edges_num > 0) {
                motifs_cnts[2] = v_prime_u_edges_num;
            }
        }
    }
    // M3: four cycle_1 v'->u'
    if (_g->nodeEdges().find(v_prime) != _g->nodeEdges().end())
    {
        if(_g->nodeEdges()[v_prime].find(u_prime) != _g->nodeEdges()[v_prime].end())
        {
            v_prime_u_prime_exist = true;
            vector<int> &v_prime_u_prime_edges_ids = _g->nodeEdges()[v_prime][u_prime];
            // find idx > v_out_edge.index()
            v_prime_u_prime_edges_left_it = upper_bound(
                v_prime_u_prime_edges_ids.begin(), v_prime_u_prime_edges_ids.end(), v_out_edge.index()
            );
            // find timestamp <= u_v_edge.time() + _delta
            v_prime_u_prime_edges_right_it = upper_bound(
                v_prime_u_prime_edges_ids.begin(), v_prime_u_prime_edges_ids.end(), u_in_edge.time() + _delta,
                [&](const time_t &a, const int &b) { return a < _g->edges()[b].time(); }
            );
            int v_prime_u_prime_edges_num = distance(v_prime_u_prime_edges_left_it, v_prime_u_prime_edges_right_it);
            if (v_prime_u_prime_edges_num > 0){
                motifs_cnts[3] = v_prime_u_prime_edges_num;
            }
        }
    }

    // M4: 4 cyecle_2 v'->u', u->u'
    // find number of pair (e3, e4) that e3.time() <= e4.time()
    // e3 in (v_prime_u_prime_edges_left_it, v_prime_u_prime_edges_right_it)
    // e4 in (u_u_prime_edges_left_it, u_u_prime_edges_right_it)
    if (_g->nodeEdges().find(u) != _g->nodeEdges().end())
    {
        if(_g->nodeEdges()[u].find(u_prime) != _g->nodeEdges()[u].end())
        {
            u_u_prime_exist = true;
            vector<int> &u_u_prime_edges_ids = _g->nodeEdges()[u][u_prime];
            // find idx > v_out_edge.index()
            u_u_prime_edges_left_it = upper_bound(
                u_u_prime_edges_ids.begin(), u_u_prime_edges_ids.end(), v_out_edge.index()
            );
            // find timestamp <= u_in_edge.time() + delta
            u_u_prime_edges_right_it = upper_bound(
                u_u_prime_edges_ids.begin(), u_u_prime_edges_ids.end(), u_in_edge.time() + _delta,
                [&](const time_t &a, const int &b) { return a < _g->edges()[b].time(); }
            );
        }
    }
    int num_M4 = 0;
    if (v_prime_u_prime_exist && u_u_prime_exist)
    {
        vector<int>::iterator v_prime_u_prime_edges_it = v_prime_u_prime_edges_left_it;
        vector<int>::iterator u_u_prime_edges_it = u_u_prime_edges_left_it;
        while(
            (v_prime_u_prime_edges_it != v_prime_u_prime_edges_right_it) &&
            (u_u_prime_edges_it != u_u_prime_edges_right_it)
        )
        {
            // keep moving u_u_prime_edges_it to right if v'u'.time() > uu'.time()
            if (_g->edges()[*v_prime_u_prime_edges_it].time() > _g->edges()[*u_u_prime_edges_it].time())
            {
                u_u_prime_edges_it++;
            }
            // it is gauranteened that v'u'.time() <= uu'.time()
            // update count and increment v_prime_u_edges_it
            else
            {
                num_M4 += distance(u_u_prime_edges_it, u_u_prime_edges_right_it);
                v_prime_u_prime_edges_it++;
            }
        }
    }
    motifs_cnts[4] = num_M4;
    
    // M5: chordal-4-cycle_1 v'->u, v'->u'
    // find number of pair(e3, e4) that e3.time() <= e4.time() 
    // e3 in (v_prime_u_edges_left_it, v_prime_u_edges_right_it)
    // e4 in (v_prime_u_prime_edges_left_it, v_prime_u_prime_edges_right_it)
    int num_chordal_4_cycle = 0;
    if(v_prime_u_exist && v_prime_u_prime_exist)
    {
        vector<int>::iterator v_prime_u_edges_it = v_prime_u_edges_left_it;
        vector<int>::iterator v_prime_u_prime_edges_it = v_prime_u_prime_edges_left_it;
        while(
            (v_prime_u_edges_it != v_prime_u_edges_right_it) &&
            (v_prime_u_prime_edges_it != v_prime_u_prime_edges_right_it) 
        )
        {
            // keep moving v_prime_u_prime_edges_it to right if v'u.time() > v'u'.time()
            if(_g->edges()[*v_prime_u_edges_it].time() > _g->edges()[*v_prime_u_prime_edges_it].time())
            {
                v_prime_u_prime_edges_it++;
            }
            // it is gauranteened that v'u.time() <= v'u'.time()
            // update count and increment v_prime_u_edges_it
            else
            {
                num_chordal_4_cycle += distance(v_prime_u_prime_edges_it, v_prime_u_prime_edges_right_it);
                v_prime_u_edges_it++;
            }
        }
    }
    motifs_cnts[5] = num_chordal_4_cycle;

    // M6: chordal-4-cycle_2 v'->u, v'->u', u->u'
    int num_M6 = 0;
    if (v_prime_u_exist && v_prime_u_prime_exist && u_u_prime_exist)
    {
        vector<int>::iterator v_prime_u_edges_it = v_prime_u_edges_left_it;
        // vector<int>::iterator v_prime_u_prime_edges_it = v_prime_u_prime_edges_left_it;
        vector<int>::iterator u_u_prime_edges_it = u_u_prime_edges_left_it;
        for (vector<int>::iterator v_prime_u_prime_edges_it = v_prime_u_prime_edges_left_it; v_prime_u_prime_edges_it != v_prime_u_prime_edges_right_it; v_prime_u_prime_edges_it++)
        {
            while(v_prime_u_edges_it != v_prime_u_edges_right_it && _g->edges()[*v_prime_u_edges_it].time() <= _g->edges()[*v_prime_u_prime_edges_it].time())
            {
                v_prime_u_edges_it++;
            }
            while(u_u_prime_edges_it != u_u_prime_edges_right_it && _g->edges()[*u_u_prime_edges_it].time() < _g->edges()[*v_prime_u_prime_edges_it].time())
            {
                u_u_prime_edges_it++;
            }
            int num_v_prime_u = distance(v_prime_u_edges_left_it, v_prime_u_edges_it);
            int num_u_u_prime = distance(u_u_prime_edges_it, u_u_prime_edges_right_it);
            num_M6 += num_v_prime_u * num_u_u_prime;
        }
    }
    motifs_cnts[6] = num_M6;

    // M7: 4-clique
    int num_4_clique = 0;
    // u'->v
    if (_g->nodeEdges().find(u_prime) != _g->nodeEdges().end())
    {
        if(_g->nodeEdges()[u_prime].find(v) != _g->nodeEdges()[u_prime].end())
        {
            u_prime_v_exist = true;
            vector<int> &u_prime_v_edges_ids = _g->nodeEdges()[u_prime][v];
            // find idx > v_out_edge.index()
            u_prime_v_edges_left_it = upper_bound(
                u_prime_v_edges_ids.begin(), u_prime_v_edges_ids.end(), v_out_edge.index()
            );
            // find timestamp <= u_v_edge.time() + _delta
            u_prime_v_edges_right_it = upper_bound(
                u_prime_v_edges_ids.begin(), u_prime_v_edges_ids.end(), u_in_edge.time() + _delta,
                [&](const time_t &a, const int &b) { return a < _g->edges()[b].time(); }
            );
        }
    }
    
    if(v_prime_u_exist && v_prime_u_prime_exist && u_prime_v_exist)
    {
        vector<int>::iterator v_prime_u_edges_it = v_prime_u_edges_left_it;
        vector<int>::iterator v_prime_u_prime_edges_it = v_prime_u_prime_edges_left_it;
        vector<int>::iterator u_prime_v_edges_it = u_prime_v_edges_left_it;
        for (vector<int>::iterator v_prime_u_prime_edges_it = v_prime_u_prime_edges_left_it; v_prime_u_prime_edges_it != v_prime_u_prime_edges_right_it; v_prime_u_prime_edges_it++)
        {
            while(v_prime_u_edges_it != v_prime_u_edges_right_it && _g->edges()[*v_prime_u_edges_it].time() <= _g->edges()[*v_prime_u_prime_edges_it].time())
            {
                v_prime_u_edges_it++;
            }
            while(u_prime_v_edges_it != u_prime_v_edges_right_it && _g->edges()[*u_prime_v_edges_it].time() < _g->edges()[*v_prime_u_prime_edges_it].time())
            {
                u_prime_v_edges_it++;
            }
            int num_v_prime_u = distance(v_prime_u_edges_left_it, v_prime_u_edges_it);
            int num_u_prime_v = distance(u_prime_v_edges_it, u_prime_v_edges_right_it);
            // cout << "num_v_prime_u * num_u_prime_v: " << num_v_prime_u * num_u_prime_v<< endl;
            num_4_clique +=  num_v_prime_u * num_u_prime_v;
        }
    }
    motifs_cnts[7] = num_4_clique;
    return motifs_cnts;
}


vector<float> GraphSearch::estimate_motif(const vector<long long int> &motifs_cnts, int num_sample, long long int W)
{
    vector<float> estimated_cnts(motifs_cnts.size(), 0);
    estimated_cnts[0] = (float) motifs_cnts[0]; // exact 3-star
    float_t W_div_k = (float_t) W / num_sample;
    for(size_t i=1; i < motifs_cnts.size(); i++)
    {
        estimated_cnts[i] = (float) motifs_cnts[i]  * W_div_k;
    }
    return estimated_cnts;
}


bool GraphSearch::out_of_range(time_t target, time_t left, time_t right)
{
    return (target < left) || (target > right);
}

long long int GraphSearch::preprocess_sampling_weights(
    vector<int>& sampling_weights, vector<vector<vector<int>::const_iterator>>& sampling_neigh_edges_it)
{
    long long int W = 0;

    #pragma omp parallel for reduction(+:W)
    for(int i=0; i < _g->numEdges(); i++)
    {
        Edge e = _g->edges()[i];
        int e_src = e.source();
        int e_dst = e.dest();
        if(e_src == e_dst){
            sampling_weights[i] = 0;
            continue;
        }
        const vector<int>& e_src_in_edges = _g->nodes()[e_src].inEdges();
        const vector<int>& e_dst_out_edges = _g->nodes()[e_dst].outEdges();
        // e_src_in_edges.index() < e.index()
        vector<int>::const_iterator e_src_in_edges_right_it = lower_bound(
            e_src_in_edges.begin(), e_src_in_edges.end(), e.index()
        );
        // e_src_in_edges.time() >= e.time() - _delta
        vector<int>::const_iterator e_src_in_edges_left_it = lower_bound(
            e_src_in_edges.begin(), e_src_in_edges.end(), e.time() - _delta,
            [&](const int &a, const time_t &b) { return _g->edges()[a].time() < b; }
        );
        int num_e_src_in_edges = distance(e_src_in_edges_left_it, e_src_in_edges_right_it);
        
        // e_dst_out_edges.index() > e.index()
        vector<int>::const_iterator e_dst_out_edges_left_it = upper_bound(
            e_dst_out_edges.begin(), e_dst_out_edges.end(), e.index()
        );
        // e_dst_out_edges.time() <= e.time() - _delta
        vector<int>::const_iterator e_dst_out_edges_right_it = upper_bound(
            e_dst_out_edges.begin(), e_dst_out_edges.end(), e.time() + _delta,
            [&](const time_t &a, const int &b) { return a < _g->edges()[b].time(); }
        );
        int num_e_dst_out_edges = distance(e_dst_out_edges_left_it, e_dst_out_edges_right_it);
        
        sampling_weights[i] = num_e_src_in_edges * num_e_dst_out_edges;
        sampling_neigh_edges_it[i] = {
            e_src_in_edges_left_it, e_src_in_edges_right_it,
            e_dst_out_edges_left_it, e_dst_out_edges_right_it
        };
        W += sampling_weights[i];
    }
    return W;
}

vector<Dependency> GraphSearch::analyze_spanning_tree(vector<vector<int>>& spanning_tree)
{
    // for edges in motif that are not in the spanning tree, just use empty Depdency
    vector<Dependency> dep_edges;
    // map from node_id to vector of edges in spanning tree that connect to the node
    unordered_map<int, vector<int>> node2edges;
    // if this edge has been used as dependent edge before
    // every edge should be at most used as dependency once
    // TODO: such assignment may not be the optimal one
    // think motif 110, two cases
    // case 1: e1 depends on e0, e4; e2 depends on e3, e4
    // case 2: e1 depends on e0; e2 depends on e1, e3, e4
    // case 2 has better sampling weights because if enforce more constraints (e1 and e4)
    unordered_set<int> visited;
    for(int i = 0; i < _h->numEdges(); i++)
    {
        dep_edges.push_back(Dependency(i));
    }
    for(int j = 0; j < spanning_tree[0].size(); j++)
    {
        int m_edge_id = spanning_tree[0][j];
        int src = _h->edges()[m_edge_id].source();
        int dst = _h->edges()[m_edge_id].dest();
        // dep_edges[m_edge_id] = Dependency(m_edge_id);
        node2edges[src].push_back(m_edge_id);
        node2edges[dst].push_back(m_edge_id);
    }
    for(int i = 1; i < spanning_tree.size(); i++)
    {
        for (int j = 0; j < spanning_tree[i].size(); j++)
        {
            int m_edge_id = spanning_tree[i][j];
            int src = _h->edges()[m_edge_id].source();
            int dst = _h->edges()[m_edge_id].dest();
            Dependency& dep = dep_edges[m_edge_id];
            for(int src_dst = 0; src_dst < 2; src_dst++)
            {
                int connected_node;
                if (src_dst == 0) // m_edge and dep_edge connect to m_edge.src
                {
                    connected_node = src;
                }
                else // m_edge and dep_edge connect to m_edge.dst
                {
                    connected_node = dst;
                }
                for (int k = 0; k < node2edges[connected_node].size(); k++)
                {
                    int dep_edge_id = node2edges[connected_node][k];
                    if (visited.find(dep_edge_id) != visited.end())
                        continue;
                    dep.add_dep_edge(
                        src_dst,
                        _h->edges()[dep_edge_id],
                        connected_node == _h->edges()[dep_edge_id].source(),
                        dep_edge_id > m_edge_id,
                        i == 1
                        );
                    visited.insert(dep_edge_id);
                }
            }
            
            node2edges[src].push_back(m_edge_id);
            node2edges[dst].push_back(m_edge_id);
        }
    }
    return dep_edges;
}

map<pair<int, int>, vector<int>> GraphSearch::analyze_exatra_edges(const Graph &h, const vector<int>& flattened_spanning_tree)
{
    map<pair<int, int>, vector<int>> hash;
    vector<int> extra_edges; // extra edges in h that are not in the spanning tree
    int sp_idx = 0;
    int h_idx = 0;
    while(h_idx < h.numEdges())
    {
        if (sp_idx < flattened_spanning_tree.size() && flattened_spanning_tree[sp_idx] == h_idx)
        {
            sp_idx++;
            h_idx++;
        }
        else
        {
            extra_edges.push_back(h_idx);
            h_idx++;
        }
    }
    // key0: largest spanning tree edge id that < extra edge id (-1 if not exist)
    // key1: smallest spanning tree edge id that > extra edge id (-1 if not exist)
    // pair<key0, key1> -> vector of extra edge id in that range
    for(int i = 0; i < extra_edges.size(); i++)
    {
        int extra_edge_id = extra_edges[i];
        int left_spanning_tree_edge_id = -1;
        int right_spanning_tree_edge_id = -1;
        for(int j = 0; j < flattened_spanning_tree.size(); j++)
        {
            if (flattened_spanning_tree[j] < extra_edge_id)
            {
                left_spanning_tree_edge_id = flattened_spanning_tree[j];
            }
            else if (flattened_spanning_tree[j] > extra_edge_id)
            {
                right_spanning_tree_edge_id = flattened_spanning_tree[j];
                break;
            }
        }
        hash[make_pair(left_spanning_tree_edge_id, right_spanning_tree_edge_id)].push_back(extra_edge_id);
    }

    return hash;
}


long long int GraphSearch::preprocess(
    vector<vector<int>> &spanning_tree, vector<Dependency> &dep_edges,
    vector<vector<long long int>>& e_sampling_weights)
{
    long long int W = 0;
    // e_sampling_weights.clear();
    for(int s = 1; s < spanning_tree.size(); s++)
    {
        // motif edges in s'th level of spanning_tree
        vector<int> &m_edges_in_level = spanning_tree[s];
        for (auto m_edge: m_edges_in_level)
            e_sampling_weights[m_edge].resize(_g->numEdges(), 0); // resize the vector
    }
    // no need to preprocess for the first level spanning tree because they just need uniform sampling
    for(int s = 1; s < spanning_tree.size(); s++)
    {
        // motif edges in s'th level of spanning_tree
        // vector<int> &m_edges_in_level = spanning_tree[s];
        #pragma omp parallel for reduction(+:W)
        for(int i = 0; i < _g->numEdges(); i++)
        {
            Edge e = _g->edges()[i];
            int u = e.source();
            int v = e.dest();
            if (u == v)
                continue;
            for(int j = 0; j < spanning_tree[s].size(); j++)
            {
                int m_edge_id = spanning_tree[s][j];
                Dependency dep = dep_edges[m_edge_id];
                // weights of the dependent edges that connect to e.source() and e.dest()
                vector<long long int> W_src_dst(2, 1);
                for(int src_dst = 0; src_dst < 2; src_dst++)
                {
                    int connected_node;
                    if (src_dst == 0) // e and dep_edge connect to e.src
                    {
                        connected_node = u;
                    }
                    else // e and dep_edge connect to e.dst
                    {
                        connected_node = v;
                    }
                    for (int k = 0; k < dep.dep_edges[src_dst].size(); k++) // dep_edge connect to edge_id.source()
                    {
                        // only consider the first dependency in the same orbit
                        // in two ranges, < and > current edge
                        // use select n from m to compute the weight
                        int less_than_cur_edge_num = 0;
                        for (int l = 0; l < dep.dep_edges[src_dst][k].size(); l++)
                        {
                            int dep_edge_id = dep.dep_edges[src_dst][k][l].index();
                            if (dep_edge_id > m_edge_id)
                                break;
                            less_than_cur_edge_num++;
                        }
                        vector<int> two_range_dep_edges_num{
                            less_than_cur_edge_num, (int) dep.dep_edges[src_dst][k].size() - less_than_cur_edge_num};
                        for(int l = 0; l < 2; l++)
                        {

                            int select_n =two_range_dep_edges_num[l];
                            if (select_n == 0)
                                continue; // no dependent edge in this range
                            int anchor_idx;
                            if (l == 0)
                                anchor_idx = 0;
                            else
                                anchor_idx = less_than_cur_edge_num;
                            Edge dep_edge = dep.dep_edges[src_dst][k][anchor_idx];
                            int dep_edge_id = dep_edge.index();
                            bool dep_edge_in_out = dep.dep_edges_in_out[src_dst][k][anchor_idx];
                            bool dep_edge_ordering = dep.dep_edges_ording[src_dst][k][anchor_idx];
                            const vector<int>* search_edges_ptr;
                            if (!dep_edge_in_out) // incoming edges
                            {
                                search_edges_ptr = &(_g->nodes()[connected_node].inEdges());
                            }
                            else // outgoing edges
                            {
                                search_edges_ptr = &(_g->nodes()[connected_node].outEdges());
                            }
                            const vector<int>& search_edges = *search_edges_ptr;
                            vector<int>::const_iterator search_edges_left_it;
                            vector<int>::const_iterator search_edges_right_it;
                            if (!dep_edge_ordering)
                            {
                                // dep_edge_id.time < m_edge_id.time
                                // binary search in [m_edge_id.time - _delta, m_edge_id.time]
                                // use lower_bound
                                // search_edges.index() < e.index()
                                search_edges_right_it = lower_bound(
                                    search_edges.begin(), search_edges.end(), e.index()
                                );
                                // search_edges.time() >= e.time() - _delta
                                search_edges_left_it = lower_bound(
                                    search_edges.begin(), search_edges.end(), e.time() - _delta,
                                    [&](const int &a, const time_t &b) { return _g->edges()[a].time() < b; }
                                );
                            }
                            else
                            {
                                // dep_edge_id.time > m_edge_id.time
                                // binary search in [m_edge_id.time, m_edge_id.time + _delta]
                                // use upper_bound
                                // search_edges.index() > e.index()
                                search_edges_left_it = upper_bound(
                                    search_edges.begin(), search_edges.end(), e.index()
                                );
                                // search_edges.time() <= e.time() + _delta
                                search_edges_right_it = upper_bound(
                                    search_edges.begin(), search_edges.end(), e.time() + _delta,
                                    [&](const time_t &a, const int &b) { return a < _g->edges()[b].time(); }
                                );
                            }
                            long long int tmp_src_dst_weight = 0;
                            if (dep_edges[dep_edge_id].num_dep_edge == 0) // the dependent edge is the leaf edge
                            {
                                int num_edges = distance(search_edges_left_it, search_edges_right_it);
                                tmp_src_dst_weight = (long long int) num_edges;
                            }
                            else // use the dependent edge's sampling weights to compute the sampling weights of edge_id
                            {
                                
                                for(vector<int>::const_iterator it = search_edges_left_it; it != search_edges_right_it; it++)
                                {
                                    tmp_src_dst_weight += e_sampling_weights[dep_edge_id][*it];
                                }
                                
                            }
                            if (select_n > 1)
                            {
                                if (tmp_src_dst_weight >= select_n)
                                {
                                    // select n from m
                                    long long int dividend = 1;
                                    long long int divisor = 1;
                                    while(select_n > 0)
                                    {
                                        dividend *= tmp_src_dst_weight;
                                        divisor *= select_n;
                                        tmp_src_dst_weight--;
                                        select_n--;
                                    }
                                    tmp_src_dst_weight = dividend / divisor;
                                }
                                else
                                {
                                    tmp_src_dst_weight = 0;
                                }
                            }
                            W_src_dst[src_dst] *= tmp_src_dst_weight;
                        }
                    }

                }
                e_sampling_weights[m_edge_id][i] = W_src_dst[0] * W_src_dst[1];
                if (s == spanning_tree.size() - 1)
                    W += e_sampling_weights[m_edge_id][i];
            }
        }
    }
    return W;
}


void GraphSearch::sampleSpanningTree(
    int iter,
    mt19937 &eng,
    discrete_distribution<> &e_center_weight_distr,
    vector<vector<long long int>>& e_sampling_weights,
    vector<vector<int>> &spanning_tree,
    vector<Dependency> &dep_edges,
    vector<Edge> &sampled_edges
)
{
    unsigned int seed = time(NULL) ^ omp_get_thread_num() ^ static_cast<unsigned int>(iter);
    int spanning_tree_num_edges = _h->numNodes() - 1;
    vector<bool> visited(spanning_tree_num_edges, false);

    // sample ordering is the reverse order of preprocessing (the spanning tree ordering)
    // sample center edge first using e_sampling_weights
    int e_center_idx_in_h = spanning_tree[spanning_tree.size() - 1][0];
    int e_center_idx_in_g = e_center_weight_distr(eng);
    // e_center_idx_in_g = 10444; // for debug, e_center is e1
    Edge e_center = _g->edges()[e_center_idx_in_g];
    sampled_edges[e_center_idx_in_h] = e_center;
    visited[e_center_idx_in_h] = true;

    // sample the dependent edges in graph of the previous level in motif
    for(int l = spanning_tree.size() - 1; l > 0; l--)
    {
        vector<int> &m_edges_in_level = spanning_tree[l];
        for(int j = 0; j < m_edges_in_level.size(); j++)
        {
            int m_edge_id = m_edges_in_level[j];
            // dependency of m_edge_id in motif
            Dependency& dep = dep_edges[m_edge_id];
            // the sampled edge in graph that serve as the m_edge_id's role in motif
            Edge sampled_edge = sampled_edges[m_edge_id];
            for(int src_dst = 0; src_dst < 2; src_dst++)
            {
                int connected_node;
                if (src_dst == 0) // e and dep_edge connect to e.src
                {
                    connected_node = sampled_edge.source();
                }
                else // e and dep_edge connect to e.dst
                {
                    connected_node = sampled_edge.dest();
                }
                for (int k = 0; k < dep.dep_edges[src_dst].size(); k++) // dep_edge connect to edge_id.source()
                {
                    // only consider the first dependency in the same orbit
                    // in two ranges, < and > current edge
                    // use select n from m to compute the weight
                    int less_than_cur_edge_num = 0;
                    for (int l = 0; l < dep.dep_edges[src_dst][k].size(); l++)
                    {
                        int dep_edge_id = dep.dep_edges[src_dst][k][l].index();
                        if (dep_edge_id > m_edge_id)
                            break;
                        less_than_cur_edge_num++;
                    }
                    vector<int> two_range_dep_edges_num{
                        less_than_cur_edge_num, (int)dep.dep_edges[src_dst][k].size() - less_than_cur_edge_num};
                    for(int l = 0; l < 2; l++)
                    {

                        int select_n =two_range_dep_edges_num[l];
                        if (select_n == 0)
                            continue; // no dependent edge in this range
                        int anchor_idx;
                        if (l == 0)
                            anchor_idx = 0;
                        else
                            anchor_idx = less_than_cur_edge_num;
                        Edge m_dep_edge = dep.dep_edges[src_dst][k][anchor_idx];
                        int m_dep_edge_id = m_dep_edge.index();
                        if (visited[m_dep_edge_id]) // if sampled before, then skip
                            continue;
                        bool m_dep_edge_in_out = dep.dep_edges_in_out[src_dst][k][anchor_idx];
                        bool m_dep_edge_ordering = dep.dep_edges_ording[src_dst][k][anchor_idx];
                        const vector<int>* search_edges_ptr;
                        if (!m_dep_edge_in_out) // incoming edges
                        {
                            search_edges_ptr = &(_g->nodes()[connected_node].inEdges());
                        }
                        else // outgoing edges
                        {
                            search_edges_ptr = &(_g->nodes()[connected_node].outEdges());
                        }
                        const vector<int>& search_edges = *search_edges_ptr;
                        vector<int>::const_iterator search_edges_left_it;
                        vector<int>::const_iterator search_edges_right_it;
                        if (!m_dep_edge_ordering)
                        {
                            // binary search in [sampled_edge.time - _delta, sampled_edge.time]
                            // use lower_bound
                            // search_edges.index() < sampled_edge.index()
                            search_edges_right_it = lower_bound(
                                search_edges.begin(), search_edges.end(), sampled_edge.index()
                            );
                            // search_edges.time() >= sampled_edge.time() - _delta
                            search_edges_left_it = lower_bound(
                                search_edges.begin(), search_edges.end(), sampled_edge.time() - _delta,
                                [&](const int &a, const time_t &b) { return _g->edges()[a].time() < b; }
                            );
                        }
                        else
                        {
                            // binary search in [sampled_edge.time, sampled_edge.time + _delta]
                            // use upper_bound
                            // search_edges.index() > sampled_edge.index()
                            search_edges_left_it = upper_bound(
                                search_edges.begin(), search_edges.end(), sampled_edge.index()
                            );
                            // search_edges.time() <= sampled_edge.time() + _delta
                            search_edges_right_it = upper_bound(
                                search_edges.begin(), search_edges.end(), sampled_edge.time() + _delta,
                                [&](const time_t &a, const int &b) { return a < _g->edges()[b].time(); }
                            );
                        }
                        int num_edges = distance(search_edges_left_it, search_edges_right_it);
                        if (num_edges == 0)
                            continue;
                        if (dep_edges[m_dep_edge_id].num_dep_edge == 0)
                        {
                            vector<int> selected;
                            if (select_n > 1)
                            {
                                // uniform sampling the first level edges
                                vector<int> search_edges_in_range(search_edges_left_it, search_edges_right_it);
                                // random_unique(search_edges_in_range.begin(), search_edges_in_range.end(), select_n, seed);
                                // sort(search_edges_in_range.begin(), search_edges_in_range.begin() + select_n);
                                selected = random_select_n(search_edges_in_range, select_n, eng);
                                sort(selected.begin(), selected.end());
                            }
                            else // uniform sample from left_it to right_it
                            {
                                int idx = rand_r(&seed) % num_edges;
                                selected = {*(search_edges_left_it+idx)};
                            }
                            for(int i = 0; i < select_n; i++)
                            {
                                int m_dep_edge_id_ith = dep.dep_edges[src_dst][k][i+anchor_idx].index();
                                // int g_dep_edge_id = search_edges_in_range[i];
                                int g_dep_edge_id = selected[i];
                                Edge g_dep_edge = _g->edges()[g_dep_edge_id];
                                sampled_edges[m_dep_edge_id_ith] = g_dep_edge;
                                visited[m_dep_edge_id_ith] = true;
                            }
                        }
                        else
                        {
                            // use the dependent edge's sampling weights to compute the sampling weights of edge_id
                            vector<long long int> cur_e_weights(num_edges, 0);
                            for(int i = 0; i < cur_e_weights.size(); i++)
                            {
                                cur_e_weights[i] = e_sampling_weights[m_dep_edge_id][*(search_edges_left_it+i)];
                            }
                            // select_n must be 1 (otherwise does not work ...)
                            discrete_distribution<> cur_e_weights_distr(cur_e_weights.begin(), cur_e_weights.end());
                            int g_dep_edge_id = *(search_edges_left_it + cur_e_weights_distr(eng));
                            Edge g_dep_edge = _g->edges()[g_dep_edge_id];
                            sampled_edges[m_dep_edge_id] = g_dep_edge;
                            visited[m_dep_edge_id] = true;
                        }
                    }
                }   
            }
        }
    }
    return;
}


bool GraphSearch::checkSpanningTree(
    vector<Edge> &sampled_edges, vector<int> &flattened_spanning_tree, vector<int> &h2gNodes)
{
    h2gNodes.resize(_h->numNodes(), -1);
    time_t upper_time;
    // Tables for mapping nodes between the two graphs
    unordered_map<int, int> g2hNodes;
    int prevEdgeId;
    for(int i = 0; i < flattened_spanning_tree.size(); i++)
    {
        int m_edge_id = flattened_spanning_tree[i];
        Edge g_e = sampled_edges[m_edge_id]; // e in graph g
        Edge h_e = _h->edges()[m_edge_id]; // e in motif h
        if (i == 0)
        {
            upper_time = g_e.time() + _delta;
        }
        // time ordering and delta constraint not satisfied
        else if(g_e.index() <= prevEdgeId || g_e.time() > upper_time)
        {
            return false;
        }
        // if different nodes in h is mapped to same node in g, invalid
        if(g2hNodes.find(g_e.source()) != g2hNodes.end() && g2hNodes[g_e.source()] != h_e.source())
            return false;
        if (g2hNodes.find(g_e.source()) == g2hNodes.end())
            g2hNodes[g_e.source()] = h_e.source();
        if(g2hNodes.find(g_e.dest()) != g2hNodes.end() && g2hNodes[g_e.dest()] != h_e.dest())
            return false;
        if (g2hNodes.find(g_e.dest()) == g2hNodes.end())
            g2hNodes[g_e.dest()] = h_e.dest();
        // if different nodes in g is mapped to same node in h, invalid
        if (h2gNodes[h_e.source()] != -1 && h2gNodes[h_e.source()] != g_e.source())
            return false;
        if (h2gNodes[h_e.source()] == -1)
            h2gNodes[h_e.source()] = g_e.source();
        if (h2gNodes[h_e.dest()] != -1 && h2gNodes[h_e.dest()] != g_e.dest())
            return false;
        if (h2gNodes[h_e.dest()] == -1)
            h2gNodes[h_e.dest()] = g_e.dest();
        prevEdgeId = g_e.index();
    }
    return true;
}

vector<int> GraphSearch::precompute2PairCounts(
    vector<vector<int>::const_iterator>& search_edges_left_its,
    vector<vector<int>::const_iterator>& search_edges_right_its)
{
    // pair_counts[i] is the number of pairs (l0_j, l1_i) where l0_j <= l1_i
    int list_num = search_edges_left_its.size();
    int l1_size = distance(search_edges_left_its[0], search_edges_right_its[0]);
    int l2_size = distance(search_edges_left_its[1], search_edges_right_its[1]);
    vector<int> pair_counts(l2_size, 0);
    vector<int>::const_iterator l1_it = search_edges_left_its[0];
    for(int i = 0; i < l2_size; i++)
    {
        int l2_val = *(search_edges_left_its[1] + i);
        // keep moving l1_it to right until l1_it >= l2_val
        while(l1_it != search_edges_right_its[0] && *l1_it < l2_val)
        {
            l1_it++;
        }
        pair_counts[i] = distance(search_edges_left_its[0], l1_it);
    }
    return pair_counts;
}

vector<int> GraphSearch::precompute2PairCounts(
    vector<vector<int>::const_iterator>& search_edges_left_its,
    vector<vector<int>::const_iterator>& search_edges_right_its,
    vector<int>& pair_counts
)
{
    int l1_size = distance(search_edges_left_its[0], search_edges_right_its[0]);
    int l2_size = distance(search_edges_left_its[1], search_edges_right_its[1]);
    vector<int> new_pair_cnts(l2_size); // pair_counts considering pair_counts and the l2 relationship
    vector<int>::const_iterator l1_it = search_edges_left_its[0];
    int sum = 0;
    for(int i = 0; i < l2_size; i++)
    {
        int l2_val = *(search_edges_left_its[1] + i);
        // find index in l1 where l1[idx] < l2_val
        while (l1_it < search_edges_right_its[0] && *l1_it < l2_val)
        {
            int l1_idx = distance(search_edges_left_its[0], l1_it);
            sum += pair_counts[l1_idx];
            l1_it++;
        }
        new_pair_cnts[i] = sum;
    }
    return new_pair_cnts;
}

long long int GraphSearch::countSortedPairs(
    vector<vector<int>::const_iterator> search_edges_left_its,
    vector<vector<int>::const_iterator> search_edges_right_its
)
{
    /*
    For sorted lists l0, l1, l2, ..., ln-1 
    with begin and end iterator as search_edges_left_it and search_edges_right_it,
    count the number of pairs (t0, t1, t2, ..., tn-1) where ti is an element in list li
    and t0 < t1 < t2 < ... < tn-1
    */
    long long int cnt = 0;
    int list_num = search_edges_left_its.size();
    // in case left iterator is already greater than right iterator
    for(int i = 0; i < list_num; i++)
    {
        if (distance(search_edges_left_its[i], search_edges_right_its[i]) <= 0)
            return 0;
    }
    // base case: only one list
    // distance of two binary searches to count the number of motif instances
    if (list_num == 1)
    {
        cnt = distance(search_edges_left_its[0], search_edges_right_its[0]);
    }
    else if(list_num == 2)
    {
        // two-pointer technique to count the number of motif instances
        // move iteration from left to right to maintain the relative temporal ordering of extra edges
        vector<vector<int>::const_iterator> search_edges_its(search_edges_left_its.size());
        for(int i = 0; i < search_edges_left_its.size(); i++)
        {
            // iterate from left to right
            search_edges_its[i] = search_edges_left_its[i];
        }
        while(
            (search_edges_its[0] != search_edges_right_its[0]) &&
            (search_edges_its[1] != search_edges_right_its[1])
        )
        {
            // keep moving search_edges_its[1] to right until it > search_edges_its[0]
            if (*search_edges_its[1] <= *search_edges_its[0])
            {
                search_edges_its[1]++;
            }
            else
            {
                // all other extra edges must > previous edge index
                cnt += distance(search_edges_its[1], search_edges_right_its[1]);
                // update search_edges_its[0]
                search_edges_its[0]++;
            }
        }
    }
    else if(list_num == 3)
    {
        // use 3-pointer technique to count the number of motif instances
        vector<vector<int>::const_iterator> search_edges_its(search_edges_left_its.size());
        for(int i = 0; i < search_edges_left_its.size(); i++)
        {
            // iterate from left to right
            search_edges_its[i] = search_edges_left_its[i];
        }
        // anchor on the second extra edge
        for(; search_edges_its[1] != search_edges_right_its[1]; search_edges_its[1]++)
        {
            // keep moving search_edges_its[0] to right until it >= search_edges_its[1]
            while(search_edges_its[0] != search_edges_right_its[0] && *search_edges_its[0] < *search_edges_its[1])
            {
                search_edges_its[0]++;
            }
            // keep moving search_edges_its[2] to right until it > search_edges_its[1]
            while(search_edges_its[2] != search_edges_right_its[2] && *search_edges_its[2] <= *search_edges_its[1])
            {
                search_edges_its[2]++;
            }
            int num_edegs_0 = distance(search_edges_left_its[0], search_edges_its[0]);
            int num_edges_2 = distance(search_edges_its[2], search_edges_right_its[2]);
            cnt += num_edegs_0 * num_edges_2;
        }
    }
    else
    {
        // for list_num > 3, use recursion
        // At 1st step, precompute the number of pairs in the first two lists l0, l1
        // and store the vector as pair_cnt0;
        // then use pair_cnt1, l1 and l2 to count number of pairs in (t1, t2, t3) in 
        // l0, l1, l2 where t1 <= t2 <= t3, store the counts in pair_cnt2
        // repeat this until we reach the last two lists, then we could use 3-pointer technique
        // to count the number of pairs of ln-3, ln-2, ln-1 (with pair_cnt_{n-3})
        vector<int> pair_cnts;
        for(int i = 1; i < list_num - 2; i++)
        {
            vector<vector<int>::const_iterator> list_left_its{search_edges_left_its[i-1], search_edges_left_its[i]};
            vector<vector<int>::const_iterator> list_right_its{search_edges_right_its[i-1], search_edges_right_its[i]};
            if (i == 1)
            {
                pair_cnts = precompute2PairCounts(list_left_its, list_right_its);
            }
            else
            {
                pair_cnts = precompute2PairCounts(list_left_its, list_right_its, pair_cnts);
            }
        }
        vector<int>::const_iterator l1_it = search_edges_left_its[list_num-3];
        vector<int>::const_iterator l3_it = search_edges_left_its[list_num-1];
        // binary search to find the number of pairs of ln-3 (t1), ln-2 (t2), ln-1 (t3)
        int cntLessThanT2 = 0; // prefix sum
        for(auto it=search_edges_left_its[list_num-2]; it != search_edges_right_its[list_num-2]; it++)
        {
            int l2_val = *it;
            // keep moving search_edges_its[list_num-3] to right until it >= search_edges_its[1]
            while(l1_it != search_edges_right_its[list_num-3] && *l1_it < l2_val)
            {
                int l1_idx = distance(search_edges_left_its[list_num-3], l1_it);
                cntLessThanT2 += pair_cnts[l1_idx];
                l1_it++;
            }
            // keep moving search_edges_its[list_num-1] to right until it > search_edges_its[1]
            while(l3_it != search_edges_right_its[list_num-1] && *l3_it <= l2_val)
            {
                l3_it++;
            }
            int cntGreaterThanT2 = distance(l3_it, search_edges_right_its[list_num-1]);
            cnt += cntLessThanT2 * cntGreaterThanT2;
        }
    }
    return cnt;
}


long long int GraphSearch::deriveMotifCounts(
    vector<Edge>& sampled_edges,
    const vector<int>& flattened_spanning_tree,
    const map<pair<int, int>, vector<int>>& sp_tree_range_edges,
    vector<int> &h2gNodes
)
{
    long long int motif_cnt = 1;
    time_t firstEdgeTime = sampled_edges[flattened_spanning_tree[0]].time();
    time_t upperTime = sampled_edges[flattened_spanning_tree[flattened_spanning_tree.size() - 1]].time();
    // for each extra edge, count the number of motif instances
    for(auto &[k, vs]: sp_tree_range_edges)
    {
        int left_spanning_tree_edge_id = k.first;
        int right_spanning_tree_edge_id = k.second;
        // every edge in extra_edges have the same upper and lower bound with respect to spanning tree edges,
        // they are sorted by index
        const vector<int>& extra_edges = vs;
        vector<vector<int>::const_iterator> search_edges_left_its;
        vector<vector<int>::const_iterator> search_edges_right_its;
        for(int h_extra_edge_id: extra_edges)
        {
            Edge h_extra_edge = _h->edges()[h_extra_edge_id];
            int h_extra_edge_src = h_extra_edge.source();
            int h_extra_edge_dst = h_extra_edge.dest();
            int g_extra_edge_src = h2gNodes[h_extra_edge_src];
            int g_extra_edge_dst = h2gNodes[h_extra_edge_dst];

            bool edge_exist = false;
            const vector<int>* search_edges_ptr;
            // the range of edges in graph g that may map to the extra edge in motif h
            if(_g->nodeEdges().find(g_extra_edge_src) != _g->nodeEdges().end())
            {
                if(_g->nodeEdges()[g_extra_edge_src].find(g_extra_edge_dst) != _g->nodeEdges()[g_extra_edge_src].end())
                {
                    edge_exist = true;
                    search_edges_ptr = &(_g->nodeEdges()[g_extra_edge_src][g_extra_edge_dst]);
                }
            }
            if (!edge_exist) // search edges is empty, no motif instance
                return 0;
            
            const vector<int>& search_edges = *search_edges_ptr;
            vector<int>::const_iterator search_edges_left_it;
            vector<int>::const_iterator search_edges_right_it;
            if (left_spanning_tree_edge_id == -1)
                // this edge is the smallest edge regarding to the spanning tree,
                // constraints by the upper bound time - delta
                // search_edges.time >= upper_time - delta
                search_edges_left_it = lower_bound(
                    search_edges.begin(), search_edges.end(), upperTime - _delta,
                    [&](const int &a, const time_t &b) { return _g->edges()[a].time() < b; }
                );
            else
                // search_edges.index > sampled_edges[left_spanning_tree_edge_id].index()
                search_edges_left_it = upper_bound(
                    search_edges.begin(), search_edges.end(), sampled_edges[left_spanning_tree_edge_id].index()
                );
            if (right_spanning_tree_edge_id == -1)
                // this edge is the largest edge regarding to the spanning tree,
                // constraints by the first edge time + delta
                // search_edges.time <= firstEdgeTime + delta
                search_edges_right_it = upper_bound(
                    search_edges.begin(), search_edges.end(), firstEdgeTime + _delta,
                    [&](const time_t &a, const int &b) { return a < _g->edges()[b].time(); }
                );
            else
                // search_edges.index < sampled_edges[right_spanning_tree_edge_id].index()
                search_edges_right_it = lower_bound(
                    search_edges.begin(), search_edges.end(), sampled_edges[right_spanning_tree_edge_id].index()
                );

            search_edges_left_its.push_back(search_edges_left_it);
            search_edges_right_its.push_back(search_edges_right_it);
        }
        long long int cur_range_cnts = countSortedPairs(search_edges_left_its, search_edges_right_its);
        motif_cnt *= cur_range_cnts;
        if (motif_cnt == 0) // early terminate if no motif instance
            return 0;
    }
    return motif_cnt;
}

long long int GraphSearch::sampleAndCheckMotifSpanningTree(
    long long int max_trial, vector<vector<long long int>>& e_sampling_weights,
    vector<vector<int>> &spanning_tree, vector<int> &flattened_spanning_tree,
    vector<Dependency> &dep_edges, map<pair<int, int>, vector<int>> sp_tree_range_edges,
    long long int &valid_sp_cnt, long long int& nz_sp_cnt
)
{
    random_device rd;
    mt19937 eng(rd() ^ omp_get_thread_num());

    valid_sp_cnt = 0;
    long long int motifs_cnt = 0;
    int spanning_tree_num_edges = _h->numNodes() - 1;
    int e_center_idx = spanning_tree[spanning_tree.size() - 1][0];
    discrete_distribution<> e_center_weight_distr(
        e_sampling_weights[e_center_idx].begin(), e_sampling_weights[e_center_idx].end());
    
    for(long long int trial = 0; trial < max_trial; trial++)
    {
        vector<int> h2gNodes(_h->numNodes(), -1); // map from node_id in h to node_id in g
        vector<Edge> sampled_edges(_h->numEdges(), Edge(-1, -1, -1, -1));
        sampleSpanningTree(
            trial, eng, e_center_weight_distr, e_sampling_weights, spanning_tree, dep_edges, sampled_edges);
        bool is_valid = checkSpanningTree(sampled_edges, flattened_spanning_tree, h2gNodes);
        // cout << "is_valid: " << is_valid << endl;
        if (is_valid)
        {
            valid_sp_cnt++;
            long long int derive_cnt = deriveMotifCounts(sampled_edges, flattened_spanning_tree, sp_tree_range_edges, h2gNodes);
            motifs_cnt += derive_cnt;
            if (derive_cnt > 0)
            {
                nz_sp_cnt++;
            }
        }
    }
    return motifs_cnt;
}

vector<float> GraphSearch::SpanningTreeSample(const Graph &g, const Graph &h,
                                                int num_of_threads, int partition_per_thread,
                                                int delta, long long int max_trial, 
                                                vector<vector<int>>& spanning_tree)

{
    // Store class data structures
    _g = &g;
    _h = &h;
    _delta = delta;

    bool debugOutput = false;

    // Stores the matching subgraphs as list of edge indices
    long long int results = 0;

    int n = _g->numNodes();
    int m = _g->numEdges();

    // the edge index of the spanning tree in motif h, sorted by index
    vector<int> flattened_spanning_tree;
    for(int i = 0; i < spanning_tree.size(); i++)
    {
        for(int j = 0; j < spanning_tree[i].size(); j++)
        {
            flattened_spanning_tree.push_back(spanning_tree[i][j]);
        }
    }
    sort(flattened_spanning_tree.begin(), flattened_spanning_tree.end());

    // vector<vector<int>>& spanning_tree represents the order of preprocessing weights
    // we will first compute weights for spanning_tree[1], then spanning_tree[2], ...
    // Note that spanning_tree[0] are the leaves of the spanning tree that just need uniform sampling
    vector<Dependency> dep_edges = analyze_spanning_tree(spanning_tree);
    int m_spanning_tree = _h->numNodes() - 1;

    for(int i =0; i < dep_edges.size(); i++)
    {
        dep_edges[i].print_info();
    }

    // store the lower and upperbound respect to spanning tree edges of the extra edges that are not in the spanning tree
    // SPTreeRangeEdges: (SPtree edge_id 1, SPtree edge_id 2, vector<int> h_edge_ids of extra edges)
    // the edges in the vector<int> have the index range of (SPtree edge_id 1, SPtree edge_id 2)
    // edges with same range are counted by n-pointer technique (linear time complexity)
    map<pair<int, int>, vector<int>> sp_tree_range_edges = analyze_exatra_edges(h, flattened_spanning_tree); 
    cout << "sp_tree_range_edges: " << endl;
    for(auto &[k, vs]: sp_tree_range_edges)
    {
        cout << "[ " << k.first << "," << k.second << " ] : ";
        for(auto v: vs)
        {
            cout << v << ", ";
        }
        cout << endl;
    }

    vector<vector<long long int>> e_sampling_weights(_h->numEdges());
    Timer t;
    t.Start();
    long long int W = preprocess(spanning_tree, dep_edges, e_sampling_weights);
    t.Stop();
    // cout << "preprocess time: " << t.Seconds() << endl;

    // cout << "W: " << W << endl;

    long long int motif_cnt = 0;
    long long int valid_sp_cnt = 0;
    long long int nz_sp_cnt = 0;
    t.Start();
    if (num_of_threads > 1)
    {
        // prepare omp
        // trial number for each thread and partition
        long long int single_trial_num = (max_trial + num_of_threads * partition_per_thread - 1) / (num_of_threads * partition_per_thread);
        // the last trial may has less trial_num than others
        long long int last_trial_num = max_trial - single_trial_num * (num_of_threads * partition_per_thread - 1);
        vector<long long int> trial_motif_cnts(num_of_threads * partition_per_thread);
        vector<long long int> valid_sp_cnts(num_of_threads * partition_per_thread);
        vector<long long int> nz_sp_cnts(num_of_threads * partition_per_thread);
        // trial_motifs_cnts[i][j] , i the parition i, j is the jth motif count

        #pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < num_of_threads * partition_per_thread; i++)
        {
            int trial_num = single_trial_num;
            if (i == num_of_threads * partition_per_thread - 1)
                trial_num = last_trial_num;
            trial_motif_cnts[i] = sampleAndCheckMotifSpanningTree(
                trial_num, e_sampling_weights, spanning_tree, flattened_spanning_tree, dep_edges, sp_tree_range_edges,
                valid_sp_cnts[i], nz_sp_cnts[i]);
        }
        // # pragma omp parallel for reduction(+:motif_cnt)
        for(int i = 0; i < num_of_threads * partition_per_thread; i++)
        {
            valid_sp_cnt += valid_sp_cnts[i];
            nz_sp_cnt += nz_sp_cnts[i];
            motif_cnt += trial_motif_cnts[i];
        }

    }
    else
    {
        motif_cnt = sampleAndCheckMotifSpanningTree(
            max_trial, e_sampling_weights, spanning_tree, flattened_spanning_tree, dep_edges, sp_tree_range_edges, valid_sp_cnt, nz_sp_cnt);
    }
    t.Stop();
    // cout << "sample time: " << t.Seconds() << endl;
    // cout << "valid_sp_cnt: " << valid_sp_cnt << endl;
    // cout << "nz_sp_cnt: " << nz_sp_cnt << endl;
    // cout << "motif_cnt: " << motif_cnt << endl;

    float_t W_div_k = (float_t) W / max_trial;
    float estimated_cnt = (float) motif_cnt  * W_div_k;

    return {estimated_cnt};
}