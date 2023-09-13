
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

        // Keeps track of the subgraph counts for each query and each delta value
        vector<vector<int>> queryDeltaCounts;

        int algo = args.algorithm();
        time_t delta = args.delta();
        if (args.verbose())
        {
            cout << "Using delta value = " << delta << endl;

            cout << "Searching for query graph in larger data graph" << endl;
        }
        if(algo == 0)
        {
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
                DataGraph g2;
                g2.setNodeAttributesDef(g.nodeAttributesDef());
                g2.setEdgeAttributesDef(g.edgeAttributesDef());
                GraphFilter::filter(g, h, criteria, g2);
                int num_of_threads = args.num_of_threads;
                if (num_of_threads * PARTITION_PER_THREAD > g2.edges().size())
                {
                    num_of_threads = 1;
                }
                cout << "Running in " << num_of_threads << " threads" << endl;
                omp_set_num_threads(num_of_threads);

                // Try each of the requested delta time restrictions
                    
                long long int limit = LONG_LONG_MAX; // No limit
                GraphSearch search;

                cout << "========== Backtracking baseline ==========" << endl;
                long long int results;
                Timer t;
                double avg_time;
                t.Start();
                if (num_of_threads > 1)
                    results = search.findOrderedSubgraphsMultiThread(g2, h, criteria, limit, delta, num_of_threads, PARTITION_PER_THREAD);
                else
                    results = search.findOrderedSubgraphs(g2, h, criteria, limit, delta);
                t.Stop();
                avg_time += t.Millisecs();
                cout << "count for " << queryFname << " : " << results << endl;
                std::cout << "total time search.findOrderedSubgraphs (ms) is: " << avg_time << " ms." << std::endl;
            }
        }
        else
        {
            GraphSearch search;
            int num_of_threads = args.num_of_threads;
            if (num_of_threads * PARTITION_PER_THREAD > g.edges().size())
            {
                num_of_threads = 1;
            }
            cout << "Running in " << num_of_threads << " threads" << endl;
            omp_set_num_threads(num_of_threads);
            if(algo == 1)
            {
                cout << "========== TEACUPS 3 path sampling (M1-M7) ==========" << endl;
                const long long int max_trial = args.max_trial();
                cout << "max_trial: " << max_trial << endl;
                vector<float> results;
                Timer t;
                double avg_time;
                t.Start();
                results = search.threePathSample(g, num_of_threads, PARTITION_PER_THREAD, delta, max_trial);
                t.Stop();
                avg_time += t.Millisecs();
                for (int i = 1; i < results.size(); i++)
                    cout << "M" << i << ": " << results[i] << ", ";
                cout << endl;
                std::cout << "total time searchPB.findOrderedSubgraphs (ms) is: " << avg_time << " ms." << std::endl;
            }
            else if(algo == 2)
            {
                cout << "========== TEACUPS 3 star (M0) ==========" << endl;
                long long int results;
                Timer t;
                double avg_time;
                t.Start();
                results = search.count3Star(g, num_of_threads, PARTITION_PER_THREAD, delta);
                t.Stop();
                avg_time += t.Millisecs();
                cout << "M0: " << results << endl;
                std::cout << "total time search.count3Star (ms) is: " << avg_time << " ms." << std::endl;
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
