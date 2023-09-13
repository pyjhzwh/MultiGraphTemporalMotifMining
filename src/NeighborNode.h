#ifndef NEIGHBORNODE_H
#define	NEIGHBORNODE_H

#include <vector>
#include "Node.h"
#include <algorithm>

/**
 * Basic node in our graph.  Contains a list of out going and incoming
 * edges in the order they occurred.
 */
class NeighborNode : public Node
{
public:
    void addOutNeighbors(int v)
    {
        if(std::find(_outNeighbors.begin(), _outNeighbors.end() ,v) == _outNeighbors.end())
            _outNeighbors.push_back(v);
    };
    void addInNeighbors(int v)
    {
        if(std::find(_inNeighbors.begin(), _inNeighbors.end() ,v) == _inNeighbors.end())
            _inNeighbors.push_back(v);
    };
    /** outgoing neighbors*/
    std::vector<int> &outNeighbors() { return _outNeighbors; }
    /** incoming neighbors*/
    std::vector<int> &inNeighbors() { return _inNeighbors; }
    /** sort outgoing neighbors by id*/
    void sortOutNeighbors() { sort(_outNeighbors.begin(), _outNeighbors.end()); }
    /** sort incoming neighbors by id*/
    void sortInNeighbors() { sort(_outNeighbors.begin(), _outNeighbors.end()); }
    std::vector<int> &Neighbors(bool direction) { return direction? inNeighbors() : outNeighbors();}
    int numOutNeighbors() const { return _outNeighbors.size(); }
    int numInNeighbors() const { return _inNeighbors.size(); }
    int numNeighbors(bool direction) { return direction? numInNeighbors() : numOutNeighbors(); }

    
    /** outgoing neighbors*/
    void addOutEdges(int e) {_outEdges.push_back(e); }
    void addInEdges(int e) {_inEdges.push_back(e); }
    std::vector<int> &edgesinDir(bool direction) { return direction? _inEdges: _outEdges;}
    const std::vector<int> &edgesinDir(bool direction) const { return direction? inEdges(): outEdges();}
    int numOutEdges() const { return _outEdges.size(); }
    int numInEdges() const { return _inEdges.size(); }
    int numEdges(bool direction) { return direction? numInEdges() : numOutEdges(); }
private:
    std::vector<int> _edges, _outEdges, _inEdges;
    std::vector<int> _outNeighbors, _inNeighbors;
};

#endif	/* NEIGHBORNODE_H */

