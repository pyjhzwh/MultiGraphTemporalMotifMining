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
#include "GraphMatch.h"

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

    std::vector<long long int> sampleAndCheckMotif(
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

    std::vector<float> sixNode112PathSample(
        const Graph &g, int num_of_threads=1, int partition_per_thread=1,
        int delta = INT_MAX, long long int max_trial = 2e5);

    long long int sixNode112PreprocessSamplingWeights(
        std::vector<int>& e1_sampling_weights,
        std::vector<int>& e3_sampling_weights,
        std::vector<int>& u_sampling_weights,
        std::vector<int>& v_sampling_weights,
        std::vector<int>& e2_sampling_weights);
    
    std::vector<long long int> sixNode112SampleAndCheckMotif(
        long long int max_trial,
        std::discrete_distribution<>& e1_weights_distr,
        std::vector<int>& e1_sampling_weights,
        std::vector<int>& e3_sampling_weights
        );

    bool sample112(
        int iter,
        std::mt19937& eng,
        std::discrete_distribution<>& e2_weights_distr,
        std::vector<int>& e1_sampling_weights,
        std::vector<int>& e3_sampling_weights,
        std::vector<Edge>& sampled_edges
        );

    std::vector<long long int> check_motif112(std::vector<Edge> &sampled_edges);
    std::vector<float> estimate_motif_general(const std::vector<long long int> &motifs_cnts, int num_sample, long long int W);

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
