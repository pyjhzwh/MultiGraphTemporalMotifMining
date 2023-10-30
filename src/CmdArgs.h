#ifndef CMDARGS_H
#define CMDARGS_H

#include <string>
#include <vector>
#include <time.h>

/**
 * Our custom command line argument parsing class.
 */
class CmdArgs
{
public:
    CmdArgs(int argc, char **argv);
    //const std::string &queryFname() const { return _queryFname; }    
    //time_t delta() const { return _delta; }
    const std::string &graphFname() const { return _graphFname; }
    const std::vector<std::string> &queryFnames() const { return _queryFnames; }
    const std::string &outFname() const { return _outFname; }
    //const std::vector<std::string> &outFnames() const { return _outFnames; }
    const time_t &delta() const { return _delta; }    
    bool success() const { return _success; }
    bool unordered() const { return _unordered; }
    bool verbose() const { return _verbose; }
    bool need_h() const { return _need_h; }
    const long long int max_trial() const { return _max_trial; }
    int algorithm() const{ return _algo; }
    int spanning_tree_no() const { return _spanning_tree_no; }
    void dispHelp() const;
     /** Parses the time in the string as seconds. Makes putting in long durations
     * less painful on the command line.
     * @param str  String containing time. 
     * If it's just a number (or ends with 's'), treat it as seconds.
     * If it ends with 'm', treat it as minutes.
     * If it ends with 'h', treat it as hours.
     * If it ends with 'd', treat it as days.
     * If it ends with 'w', treat it as weeks.*/
    int parseDuration(std::string str) const;
    /**
     * Creates a default output file name based on the input filenames and the delta used.
     * @param gFname  The data graph we are searching against.
     * @param hFname  The query graph we are looking for.
     * @param delta  Time in seconds that the query graph must take place over.
     * @return New unique filename based on the input parameters.
     */
    std::string createOutFname(const std::string &gFname, const std::string &hFname, time_t delta);
    int num_of_threads;

    // parse order "1,2,3,0" to vector {1, 2, 3, 0}
    std::vector<int> parseOrder(std::string str) const;

private:
    std::string _graphFname, _outFname; // _queryFname
    std::vector<std::string> _queryFnames; // _outFnames;
    time_t _delta;
    bool _success, _unordered, _verbose;
    int _algo, _spanning_tree_no;
    long long int _max_trial;
    bool _need_h; // if query file is needed (true if algo=0, false otherwise)
};

#endif
