/*
 * File:   GraphSearch.h
 * Author: D3M430
 *
 * Created on January 13, 2017, 3:45 PM
 */

#ifndef GRAPHSEARCH_H
#define GRAPHSEARCH_H

#include "Graph.h"
#include "MatchCriteria.h"
#include <limits.h>
#include <vector>
#include <stack>
#include <unordered_set>
#include <tuple>
#include <random>
#include <functional>
#include <memory>
#include "GraphMatch.h"

class Dependency
{
public:
    // int vertex_id;
    int edge_id;
    int num_dep_edge;
    std::vector<std::vector<std::vector<Edge>>> dep_edges; // edge is dependent on dep_edge
    // dep_edge is in or out of vertex
    std::vector<std::vector<std::vector<bool>>> dep_edges_in_out; // 0: in, 1: out
    // dep_edge.time < or > edge.time
    std::vector<std::vector<std::vector<bool>>> dep_edges_ording; // 0: smaller, 1: larger
    // Dependency(int vertex_id, int edge_id) : vertex_id(vertex_id), edge_id(edge_id) {};
    Dependency(int edge_id = -1) : edge_id(edge_id), num_dep_edge(0) {
        dep_edges.resize(2);
        dep_edges_in_out.resize(2);
        dep_edges_ording.resize(2);
    };

    void add_dep_edge(int src_dst, Edge dep_edge, bool in_out, bool edge_ording, bool check_same = false)
    {
        int check_vertex;
        if (check_same)
        {
            check_vertex = in_out ? dep_edge.source() : dep_edge.dest();
        }
        else
        {
            check_vertex = in_out ? dep_edge.dest() : dep_edge.source();
        }
        for(int i = 0; i < dep_edges[src_dst].size(); i++)
        {
            int other_vertex;
            if (check_same)
            {
                other_vertex = in_out ? dep_edges[src_dst][i][0].source() : dep_edges[src_dst][i][0].dest();
            }
            else
            {
                other_vertex = in_out ? dep_edges[src_dst][i][0].dest() : dep_edges[src_dst][i][0].source();
            }
            if (check_vertex == other_vertex)
            {
                dep_edges[src_dst][i].push_back(dep_edge);
                dep_edges_in_out[src_dst][i].push_back(in_out);
                dep_edges_ording[src_dst][i].push_back(edge_ording);
                num_dep_edge++;
                return;
            }
        }
        dep_edges[src_dst].push_back({dep_edge});
        dep_edges_in_out[src_dst].push_back({in_out});
        dep_edges_ording[src_dst].push_back({edge_ording});
        num_dep_edge++;
    }

    void print_info() {
        std::cout << "edge " << edge_id << " has " << num_dep_edge << " dependencies." << std::endl;
        if (num_dep_edge == 0)
            return;
        for(int j = 0; j < dep_edges.size(); j++) // src dep or dst dep
        {
            if (j == 0)
                std::cout << "src deps: ";
            else
                std::cout << "dst deps: ";
            for(int k = 0; k < dep_edges[j].size(); k++)
            {
                std::cout << "{";
                for(int l = 0; l < dep_edges[j][k].size(); l++)
                {
                    bool in_out = dep_edges_in_out[j][k][l];
                    bool ordering = dep_edges_ording[j][k][l];
                    std::string in_out_str = in_out ? "out" : "in";
                    std::string ordering_str = ordering ? ">" : "<";
                    std::cout << "(" << dep_edges[j][k][l].index()<< ", " \
                         << in_out_str << ", " << ordering_str << "), ";
                }
                std::cout << "}, ";
                
            }
            std::cout << std::endl;
        }
    }
};

/**
 * Main class for performing subgraph searches.
 */
class GraphSearch
{
public:
    /**
     * Performs a subgraph search, in which the ORDER of the edges between the
     * query graph and original graph must match.  If a -> b comes before b -> c
     * in the query graph, it must also come before it in the original graph.
     * ALL matching subgraphs are returned.  The returned data structure is a list
     * of all the subgraphs (which in this case are indices to the edges in g)
     * that match h.  If empty, then no match was found.
     * @param g  The directed graph to search on.
     * @param h  The directed query graph we are trying to match.
     * @param limit  The max number of subgraphs to find.
     * @param delta  The max time duration allowed between edge matches.
     * @return  List of subgraphs that match h.
     */
    long long int findOrderedSubgraphs(const Graph &g, const Graph &h, long long int limit = LONG_LONG_MAX, int delta = INT_MAX);

    /**
     * Performs a subgraph search, in which the ORDER of the edges between the
     * query graph and original graph must match.  If a -> b comes before b -> c
     * in the query graph, it must also come before it in the original graph.
     * ALL matching subgraphs are returned.  The returned data structure is a list
     * of all the subgraphs (which in this case are indices to the edges in g)
     * that match h.  If empty, then no match was found.
     * @param g  The directed graph to search on.
     * @param h  The directed query graph we are trying to match.
     * @param criteria  Polymorphic class specifies whether or not two given edges match the query criteria.
     * @param limit  The max number of subgraphs to find.
     * @param delta  The max time duration allowed between edge matches.
     * @return  List of subgraphs that match h.
     */
    long long int findOrderedSubgraphs(const Graph &g, const Graph &h, const MatchCriteria &criteria, long long int limit = LONG_LONG_MAX, int delta = INT_MAX);
    long long int findOrderedSubgraphsMultiThread(
        const Graph &g, const Graph &h, const MatchCriteria &criteria, long long int limit = LONG_LONG_MAX, int delta = INT_MAX, int num_of_threads=1, int partition_per_thread=1);
    void findOrderedSubgraphs(
        long long int *results, int start_edge_idx, int end_edge_idx,
        const Graph &g, const Graph &h, const MatchCriteria &criteria, long long int limit = LONG_LONG_MAX, int delta = INT_MAX);
    
    long long int count3Star(
        const Graph &g, int num_of_threads, int partition_per_thread, int delta);


    // three path sampling for M1 - M7
    std::vector<float> threePathSample(
        const Graph &g, int num_of_threads=1, int partition_per_thread=1,
        int delta = INT_MAX, long long int max_trial = 2e5);

    std::vector<long long int> threePathSampleAndCheckMotif(
        long long int max_trial, std::discrete_distribution<>& edge_in_mult_out_weights,
        std::vector<std::vector<std::vector<int>::const_iterator>>& sampling_neigh_edges_it);

    /**
     * Performs a subgraph search, in which the edge order does NOT matter,
     * and date/time of events is usually ignored.
     * The returned data structure is a list of all the subgraphs (which in this
     * case are indices to the edges in g) that match h.  If empty, then no match was found.
     * @param g  The directed graph to search on.
     * @param h  The directed query graph we are trying to match.
     * @param limit  The max number of subgraphs to find.
     * @return  List of subgraphs that match h.
     */
    std::vector<GraphMatch> findAllSubgraphs(const Graph &g, const Graph &h, long long int limit = LONG_LONG_MAX);

    /**
     * Performs a subgraph search, in which the edge order does NOT matter,
     * and date/time of events is usually ignored.
     * The returned data structure is a list of all the subgraphs (which in this
     * case are indices to the edges in g) that match h.  If empty, then no match was found.
     * @param g  The directed graph to search on.
     * @param h  The directed query graph we are trying to match.
     * @param criteria  Polymorphic class specifies whether or not two given edges match the query criteria.
     * @param limit  The max number of subgraphs to find.
     * @return  List of subgraphs that match h.
     */
    std::vector<GraphMatch> findAllSubgraphs(const Graph &g, const Graph &h, const MatchCriteria &criteria, long long int limit = LONG_LONG_MAX);

    bool sample_directed_3path(
        int edge_id, std::vector<std::vector<std::vector<int>::const_iterator>>& sampling_neigh_edges_it,
        std::vector<Edge>& sampled_edges, std::vector<long long int>& failure);
    std::vector<long long int> check_motif(std::vector<Edge> &sampled_edges);
    std::vector<float> estimate_motif(const std::vector<long long int> &motifs_cnts, int num_sample, long long int W);
    bool out_of_range(time_t target, time_t left, time_t right);
    long long int preprocess_sampling_weights(
        std::vector<int>& sampling_weights, std::vector<std::vector<std::vector<int>::const_iterator>>& sampling_neigh_edges_it);

    std::vector<float> sixNodePathSample(
        const Graph &g, int spanning_tree_no=112, int num_of_threads=1, int partition_per_thread=1,
        int delta = INT_MAX, long long int max_trial = 2e5);

    long long int sixNode108PreprocessSamplingWeights(
        std::vector<long long int>& e1_sampling_weights,
        std::vector<long long int>& e2_sampling_weights,
        std::vector<long long int>& e3_sampling_weights);

    long long int sixNode109PreprocessSamplingWeights(
        std::vector<long long int>& e1_sampling_weights,
        std::vector<long long int>& e2_sampling_weights,
        std::vector<long long int>& e3_sampling_weights);

    long long int sixNode110PreprocessSamplingWeights(
        std::vector<long long int>& e1_sampling_weights,
        std::vector<long long int>& e2_sampling_weights,
        std::vector<long long int>& e3_sampling_weights);

    long long int sixNode111PreprocessSamplingWeights(
        std::vector<long long int>& e1_sampling_weights,
        std::vector<long long int>& e2_sampling_weights,
        std::vector<long long int>& e3_sampling_weights);

    long long int sixNode112PreprocessSamplingWeights(
        std::vector<long long int>& e1_sampling_weights,
        std::vector<long long int>& e2_sampling_weights,
        std::vector<long long int>& e3_sampling_weights);
    
    std::vector<long long int> sixNodeSampleAndCheckMotif(
        long long int max_trial,
        int spanning_tree_no,
        std::unique_ptr<std::discrete_distribution<>>& e_center_weight_distr,
        std::vector<long long int>* e_left_sampling_weights,
        std::vector<long long int>* e_right_sampling_weights
        );

    bool sample108(
        int iter,
        std::mt19937& eng,
        std::unique_ptr<std::discrete_distribution<>>& e3_weights_distr,
        std::vector<long long int>* ei_sampling_weights_ptr,
        std::vector<long long int>* ej_sampling_weights_ptr,
        std::vector<Edge>& sampled_edges
        );

    bool sample109(
        int iter,
        std::mt19937& eng,
        std::unique_ptr<std::discrete_distribution<>>& e3_weights_distr,
        std::vector<long long int>* ei_sampling_weights_ptr,
        std::vector<long long int>* ej_sampling_weights_ptr,
        std::vector<Edge>& sampled_edges
        );

    bool sample110(
        int iter,
        std::mt19937& eng,
        std::unique_ptr<std::discrete_distribution<>>& e1_weights_distr,
        std::vector<long long int>* ei_sampling_weights_ptr,
        std::vector<long long int>* e2_sampling_weights_ptr,
        std::vector<Edge>& sampled_edges
        );

    bool sample111(
        int iter,
        std::mt19937& eng,
        std::unique_ptr<std::discrete_distribution<>>& e2_weights_distr,
        std::vector<long long int>* ei_sampling_weights_ptr,
        std::vector<long long int>* e3_sampling_weights_ptr,
        std::vector<Edge>& sampled_edges
        );

    bool sample112(
        int iter,
        std::mt19937& eng,
        std::unique_ptr<std::discrete_distribution<>>& e2_weights_distr,
        std::vector<long long int>* e1_sampling_weights_ptr,
        std::vector<long long int>* e3_sampling_weights_ptr,
        std::vector<Edge>& sampled_edges
        );

    std::vector<long long int> check_motif108(std::vector<Edge> &sampled_edges);
    std::vector<long long int> check_motif109(std::vector<Edge> &sampled_edges);
    std::vector<long long int> check_motif110(std::vector<Edge> &sampled_edges);
    std::vector<long long int> check_motif111(std::vector<Edge> &sampled_edges);
    std::vector<long long int> check_motif112(std::vector<Edge> &sampled_edges);

    std::vector<float> estimate_motif_general(
        const std::vector<long long int> &motifs_cnts, long long int num_sample, long long int W);

    std::map<int,std::function<long long int(std::vector<long long int>&, std::vector<long long int>&, std::vector<long long int>&)>> preprocess_funcs {
            { 108, std::bind(&GraphSearch::sixNode108PreprocessSamplingWeights, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3)},
            { 109, std::bind(&GraphSearch::sixNode109PreprocessSamplingWeights, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3)},
            { 110, std::bind(&GraphSearch::sixNode110PreprocessSamplingWeights, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3)},
            { 111, std::bind(&GraphSearch::sixNode111PreprocessSamplingWeights, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3)},
            { 112, std::bind(&GraphSearch::sixNode112PreprocessSamplingWeights, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3)},
        };

    std::map<int,std::function<bool(
        int, std::mt19937&, std::unique_ptr<std::discrete_distribution<>>&,
        std::vector<long long int>*, std::vector<long long int>*, std::vector<Edge>&)>> sample_funcs {
            { 108, std::bind(&GraphSearch::sample108, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6)},
            { 109, std::bind(&GraphSearch::sample109, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6)},
            { 110, std::bind(&GraphSearch::sample110, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6)},
            { 111, std::bind(&GraphSearch::sample111, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6)},
            { 112, std::bind(&GraphSearch::sample112, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6)},
        };

    std::map<int,std::function<std::vector<long long int>(std::vector<Edge> &)>> check_motif_funcs {
            { 108, std::bind(&GraphSearch::check_motif108, this, std::placeholders::_1)},
            { 109, std::bind(&GraphSearch::check_motif109, this, std::placeholders::_1)},
            { 110, std::bind(&GraphSearch::check_motif110, this, std::placeholders::_1)},
            { 111, std::bind(&GraphSearch::check_motif111, this, std::placeholders::_1)},
            { 112, std::bind(&GraphSearch::check_motif112, this, std::placeholders::_1)},
        };

    std::vector<Dependency> analyze_spanning_tree(std::vector<std::vector<int>>& spanning_tree);

    long long int preprocess(
        std::vector<std::vector<int>> &spanning_tree, std::vector<Dependency> &dep_edges,
        std::vector<std::vector<long long int>>& e_sampling_weights);

    std::vector<Edge> sampleSpanningTree(
        int iter,
        std::mt19937 &eng,
        std::discrete_distribution<> &e_center_weight_distr,
        std::vector<std::vector<long long int>>& e_sampling_weights,
        std::vector<std::vector<int>> &spanning_tree,
        std::vector<Dependency> &dep_edges
        );

    bool checkSpanningTree(
        std::vector<Edge> &sampled_edges,
        std::vector<int> &flattened_spanning_tree
        );

    std::vector<long long int> sampleAndCheckMotifSpanningTree(
        long long int max_trial,
        std::vector<std::vector<long long int>>& e_sampling_weights,
        std::vector<std::vector<int>> &spanning_tree,
        std::vector<Dependency> &dep_edges
    );

    std::vector<float> SpanningTreeSample(const Graph &g, const Graph &h,
                                                int num_of_threads, int partition_per_thread,
                                                int delta, long long int max_trial, 
                                                std::vector<std::vector<int>>& spanning_tree);
    

private:
    /** Creates map of which nodes in G can map to the nodes we are searching for from H */
    std::vector<std::unordered_set<int>> mapPossibleNodes();

    /** Performs recursive unordered graph search, stopping at first matching subgraph */
    bool search(int &numAssigned, std::vector<std::unordered_set<int>> &h2gPossible);

    /** Performs recursive unordered graph search, storing each matching subgraph in results */
    bool search(int &numAssigned, std::vector<std::unordered_set<int>> &h2gPossible, std::vector<GraphMatch> &results);

    /** Returns true if the number assigned all match up with the appropriate edges */
    bool matchesSoFar(int numAssigned);

    /** Picks an efficient list of edge indexes to search from before searching through them
     * for a edge that matches query edge h_i. The return value is the index of
     * the matching edge in G.  If no edge is found, it will return the size of edges in G. */
    int findNextMatch(int h_i, int startIndex);

    /** Searches through edge indexes listed in edgesToSearch (starting at the startIndex)
     * for a edge that matches query edge h_i. The return value is the index of
     * the matching edge in G.  If no edge is found, it will return the size of edges in G. */
    int findNextMatch(int h_i, const std::vector<int> &edgesToSearch, int startIndex);

    /**
     * Performs binary search to find best starting place.
     * @param g_i  Edge index we want to find (or greater).
     * @param edgeIndexes  List to search through.
     * @return   Index into the list where there is an edge >= g_i.
     */
    int findStart(int g_i, const std::vector<int> &edgeIndexes);

    /** Converts the given stack into a vector, without modifying it */
    std::vector<int> convert(std::stack<int> s);

    /** Converts the given stack and final edge to a GraphMatch object,
     without modifying the stack. */
    GraphMatch convert(const std::stack<int> &s, int g_lastEdge);

    // count 3-star
    long long int count_stars(int delta, int id_begin, int id_end);
    long long int count_stars_multi_thread(int delta, int num_of_threads, int partition_per_thread);

    // helper function of find duplicates
    bool containDuplicates(std::vector<int> list);

    // choose num_random of the list between begin and end
    template<class BidiIter > BidiIter random_unique(BidiIter begin, BidiIter end, size_t num_random, unsigned int seed);

    std::vector<int> random_select_n(std::vector<int>& list, size_t num_random, unsigned int seed);

    // Private data members
    const Graph *_g, *_h;
    const MatchCriteria *_criteria;
    int _delta;
    time_t _firstEdgeTime;
    int _firstEdgeid;
    std::vector<int> _h2gNodes, _g2hNodes;
    std::vector<int> _numSearchEdgesForNode;
    std::stack<int> _sg_edgeStack; //, _h_edgeStack;
    std::vector<int> _allEdges;

    time_t _upper_time; // upper bound of the g.edges()[g_i].time() when we run out of edges
};

#endif /* GRAPHSEARCH_H */
