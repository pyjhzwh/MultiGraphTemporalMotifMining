# Accurate and Fast Estimation of Temporal Motifs using Path Sampling

This is the code for the paper "Accurate and Fast Estimation of Temporal Motifs using Path Sampling".

## Build

Change the c++ compiler in line 13 in `Makefile` to that works for you. Then run
```
make -j
```
It will compile to an executable file `graph_search`.

## Prepare the inputs
`cd dataset`, modify the DATASET name in line 3 and/or the download path in line 7 of `download.sh`. Then run `bash download.sh`. It will download the required dataset, converted to gdf format in `examples/non-attribute/` directory.

## Run `graph_search`

```
./graph_search -g {input graph path} -q {motif path} -delta {delta value} -algo {algorithm} -thread {num-of-threads} 
```

* -g: path to input graph (required)
* -q: path to motif graph (required if algo=0)
* -delta: 1M(minute), 1H (hour), 1D (day), 1W(week), ... or just the number like 86400, the default unit is seconds
* -algo: 0, 1, 2 or 3; 0 means running the backtracking exact count baseline [BT](https://ieeexplore.ieee.org/abstract/document/8622100); 1 means running hard-coded 6-node motif counting algorithm to get counts; 2 means running spanning tree sampling algorithm; 3 means using pointer-technique to quickly get exact count of motifs that are multi-graphs
* -thread: thread number, default value is 1
* -max_trial: the number of samples to take in TEACUPS 3-path-sampling algorithm, only takes effect for algo=1; default value is "1e3". Note if you use scientific notation, make sure surround the number by "".


A sample output of `./graph_search -g {path-to-wikitalk-temporal} -delta 4W -algo 1 -thread 32 -max_trial "2e7"` is the following
```
Reading node data
Reading edge data
Reading node data
Reading edge data
Filtering graph
Running in 32 threads
========== TEACUPS 3 path sampling (M1-M7) ==========
max_trial: 20000000
M1: 1.04415e+10, M2: 1.59627e+09, M3: 3.05108e+08, M4: 1.1191e+09, M5: 7.381e+08, M6: 3.4886e+09, M7: 1.07794e+09, 
total time searchPB.findOrderedSubgraphs (ms) is: 2228.94 ms.
```

## Reproduce results

<!-- To reproduce the estimated counts and runtime of TEACUPS 3-path sampling for M1-M7, and counts and runtime of efficient counting M0 -->
To reproduce the runtime and counts for BT and TEACUPS, run `python reproduce.py`.
`RUN_TEACUPS` and `RUN_BT` is boolean variable to decide whether to run TEACUPS and BT, respectively. Be aware to set `RUN_BT = True` because BT takes days or weeks.
The path to the directory of input graphs and motifs can be set in `GRAPH_DIR` and `MOTIF_DIR`, respectively.

Part of the sample output of `python reproduce.py` is as follows:
```
graph=wiki-talk-temporal, delta=4W
                    |     M0    |    M1    |    M2    |    M3    |    M4    |    M5    |    M6    |    M7    
TEACUPS cnt         |   7.09e+11|  1.04e+10|  1.58e+09|  2.93e+08|  1.01e+09|  7.22e+08|  3.53e+09|  1.23e+09
TEACUPS runtime (s) |       1.15|                                2.46
graph=wiki-talk-temporal, delta=8W
                    |     M0    |    M1    |    M2    |    M3    |    M4    |    M5    |    M6    |    M7    
TEACUPS cnt         |   2.16e+12|  3.39e+10|  7.78e+09|  1.51e+09|  9.51e+09|  6.75e+09|  6.23e+10|  2.14e+10
TEACUPS runtime (s) |       2.02|                                2.25
graph=wiki-talk-temporal, delta=16W
                    |     M0    |    M1    |    M2    |    M3    |    M4    |    M5    |    M6    |    M7    
TEACUPS cnt         |   6.04e+12|  1.01e+11|  3.43e+10|  6.79e+09|  7.74e+10|  5.41e+10|  9.91e+11|  4.40e+11
TEACUPS runtime (s) |       4.16|                                2.25
```
