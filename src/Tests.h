#include "CertGraph.h"
#include "GraphMatch.h"

#ifndef TESTS__H
#define TESTS__H

/**
 * Functions for testing results from searches, etc.
 */
class Tests
{
public:
    static int countDuplicateEdges(const CertGraph &g);

    /** Sort Graph Match gEdges by the order of hEdges*/
    static std::vector<std::vector<int>> sortGraphMatch(std::vector<GraphMatch> &gms);

    static bool assertSameGraphMatch(std::vector<GraphMatch> gm1, std::vector<GraphMatch>gm2, Graph g);
};

#endif

