
#include <iostream>
#include "CmdArgs.h"
#include "DataGraph.h"
#include "FileIO.h"
#include "GraphFilter.h"
#include "GraphSearch.h"
#include "MatchCriteria_DataGraph.h"
#include "timer.h"
#include "Tests.h"
#include "omp.h"

using namespace std;

int main(int argc, char **argv)
{
    try
    {
        // Parse command line arguments
        CmdArgs args(argc, argv);
        if (!args.success())
            return -1;

        // Largest graph to display on console (for testing purposes)
        const int MAX_NUM_EDGES_FOR_DISP = 50;
        const int PARTITION_PER_THREAD = 4;

        if (args.verbose())
        {
            cout << "Loading data graph from " << args.graphFname() << endl;
        }
        DataGraph g = FileIO::loadGenericGDF(args.graphFname());
        if (args.verbose())
        {
            cout << g.nodes().size() << " nodes, " << g.edges().size() << " edges" << endl;
            if (g.numEdges() < MAX_NUM_EDGES_FOR_DISP)
                g.disp();
            cout << endl;
        }
        long long int limit = LONG_LONG_MAX; // No limit

        // Keeps track of the subgraph counts for each query and each delta value
        vector<vector<int>> queryDeltaCounts;

        int algo = args.algorithm();
        int spanning_tree_no = args.spanning_tree_no();
        time_t delta = args.delta();
        if (args.verbose())
        {
            cout << "Using delta value = " << delta << endl;

            cout << "Searching for query graph in larger data graph" << endl;
        }
        
        // vector<int> cnts(args.queryFnames().size(), 0);
        // vector<int> runtimes(args.queryFnames().size(), 0);
        // Try each of the requested query graphs
        for (int i = 0; i < args.queryFnames().size(); i++)
        {
            const string &queryFname = args.queryFnames()[i];
            if (args.verbose())
                cout << "Loading query graph from " << queryFname << endl;
            DataGraph h;
            h = FileIO::loadGenericGDF(queryFname);
            if (args.verbose())
            {
                h.disp();
                cout << h.nodes().size() << " nodes, " << h.edges().size() << " edges" << endl;
                if (h.numEdges() < MAX_NUM_EDGES_FOR_DISP)
                    h.disp();
                cout << endl;
            }

            if (h.nodeAttributesDef() != g.nodeAttributesDef())
                throw "Node attribute definitions don't match between the query graph and data graph.";
            if (h.edgeAttributesDef() != h.edgeAttributesDef())
                throw "Edge attribute definitionsgit p don't match between the query graph and data graph.";

            MatchCriteria_DataGraph criteria;
            if (args.verbose())
                cout << "Filtering data graph to improve query performance." << endl;
            // DataGraph g2;
            // g2.setNodeAttributesDef(g.nodeAttributesDef());
            // g2.setEdgeAttributesDef(g.edgeAttributesDef());
            // GraphFilter::filter(g, h, criteria, g2);
            g.forceUpdateOrderedEdges();
            h.forceUpdateOrderedEdges();
            int num_of_threads = args.num_of_threads;
            if (num_of_threads * PARTITION_PER_THREAD > g.edges().size())
            {
                num_of_threads = 1;
            }
            cout << "Running in " << num_of_threads << " threads" << endl;
            omp_set_num_threads(num_of_threads);

            // Try each of the requested delta time restrictions
            
            GraphSearch search;

            if(algo == 0)
            {
                cout << "========== Backtracking baseline ==========" << endl;
                long long int results;
                Timer t;
                double avg_time;
                t.Start();
                if (num_of_threads > 1)
                    results = search.findOrderedSubgraphsMultiThread(g, h, criteria, limit, delta, num_of_threads, PARTITION_PER_THREAD);
                else
                    results = search.findOrderedSubgraphs(g, h, criteria, limit, delta);
                t.Stop();
                avg_time += t.Seconds();
                cout << "count for " << queryFname << " : " << results << endl;
                std::cout << "BT runtime(s): " << avg_time << " s." << std::endl;
            }
            else if (algo == 1)
            {
                cout << "========== Spanning tree sampling ==========" << endl;
                const long long int max_trial = args.max_trial();
                cout << "max_trial: " << max_trial << endl;
                vector<float> results;
                Timer t;
                double avg_time = 0;
                vector<vector<int>> spanning_tree = h.spanning_tree();//{{0,4}, {1}, {3}};
                if (spanning_tree.size() == 0) // if no spanning tree is provided, use heuristics to generate one
                {
                    cout << "No spanning tree provided, using heuristics to generate one" << endl;
                    spanning_tree = search.analyzeSPTreeSample(h);
                }
                // print spanning tree
                cout << "spanning tree: ";
                for (int i = 0; i < spanning_tree.size(); i++)
                {
                    cout << "{";
                    for (int j = 0; j < spanning_tree[i].size(); j++)
                    {
                        cout << spanning_tree[i][j] << ", ";
                    }
                    cout << "},";
                }
                cout << endl;
                t.Start();
                // vector<vector<int>> spanning_tree = search.analyzeSPTreeSample(h); //{{4, 9} ,{0}};
                results = search.SpanningTreeSample(
                    g, h, num_of_threads, PARTITION_PER_THREAD, delta, max_trial, spanning_tree);
                t.Stop();
                avg_time += t.Seconds();
                cout << "count for " << queryFname << " : " << results[0] << endl;
                std::cout << "Sampling runtime(s): " << avg_time << " s." << std::endl;
            }
            else if (algo == 2)
            {
                cout << "========== DFS with pointer-technique (multi-graph motifs) ==========" << endl;
                long long int results;
                Timer t;
                double avg_time = 0;
                t.Start();
                vector<int> spanning_tree = search.analyzeSPTreeBT(h);
                if (num_of_threads > 1)
                    results = search.findOrderedSubgraphsSpanningTreeMultiThread(g, h, criteria, spanning_tree, limit, delta, num_of_threads, PARTITION_PER_THREAD);
                else
                    results = search.findOrderedSubgraphsSpanningTree(g, h, criteria, spanning_tree, limit, delta);
                t.Stop();
                avg_time += t.Seconds();
                cout << "count for " << queryFname << " : " << results << endl;
                std::cout << "BT+derivecnts runtime(s): " << avg_time << " s." << std::endl;
            }
        }
        
        if (args.verbose())
            cout << "Done!\n"
                 << endl;

    }
    catch (exception &e)
    {
        cout << "An error occurred: " << e.what() << endl;
    }
    catch (const char *msg)
    {
        cout << "An error occurred: " << msg << endl;
    }
    catch (...)
    {
        cout << "An unknown exception occurred." << endl;
    }
    return 0;
}
