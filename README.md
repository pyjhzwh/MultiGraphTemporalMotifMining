# Accurate and Fast Estimation of Temporal Motifs using Path Sampling

This is the code for the paper "Accurate and Fast Estimation of Temporal Motifs using Path Sampling" (ICDM'24).
The camera-ready paper for ICDM'24 is coming soon.
For author's version, please refer to: https://arxiv.org/abs/2409.08975

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
* -q: path to motif graph (required)
* -delta: 1M(minute), 1H (hour), 1D (day), 1W(week), ... or just the number like 86400, the default unit is seconds
* -algo: 0, 1; 0 means running the backtracking exact count baseline [BT](https://ieeexplore.ieee.org/abstract/document/8622100); 1 means running the sampling-based algorithm;
* -thread: thread number, default value is 1
* -max_trial: the number of samples to take in STEMMING 3-path-sampling algorithm, only takes effect for algo=1; default value is "1e3". Note if you use scientific notation, make sure surround the number by "".


A sample output of `./graph_search -g {path-to-wikitalk-temporal} -q {path-to-M3-1} -delta 1W -algo 1 -thread 32 -max_trial "1e8"` is the following
```
GRAPH=wiki-talk-temporal MOTIF=M3-1 delta=1W
Reading node data
Reading edge data
Reading node data
Reading edge data
Filtering graph
Running in 32 threads
========== Spanning tree sampling ==========
max_trial: 100000000
Choose spanning tree: 0, 4, 
Choose spanning tree: {0, }, {4, }, 
edge 0 has 0 dependencies.
edge 1 has 0 dependencies.
edge 2 has 0 dependencies.
edge 3 has 0 dependencies.
edge 4 has 1 dependencies.
src deps: {(0, in, <), }, 
dst deps: 
edge 5 has 0 dependencies.
edge 6 has 0 dependencies.
edge 7 has 0 dependencies.
edge 8 has 0 dependencies.
sp_tree_range_edges: 
[ 0,4 ] : 1, 2, 3, 
[ 4,-1 ] : 5, 6, 7, 8, 
count for ./examples/non-attributed/M3-1.gdf : 2.55813e+08
Sampling runtime(s): 11.0727 s.
```

## Reproduce results

<!-- To reproduce the estimated counts and runtime of STEMMING 3-path sampling for M1-M7, and counts and runtime of efficient counting M0 -->
To reproduce the runtime and counts for BT and STEMMING, run `python reproduce.py`.
`RUN_STEMMING` and `RUN_BT` are boolean variables to decide whether to run STEMMING, BT+derivetns and BT, respectively. Be aware to set `RUN_BT = True` because BT takes days or weeks.
The path to the directory of input graphs and motifs can be set in `GRAPH_DIR` and `MOTIF_DIR`, respectively.

Part of the sample output of `python reproduce.py` is as follows:
```
graph=wiki-talk-temporal, delta=4W
                     |    M3-0   |   M3-1   |   M4-0   |   M4-1   |   M4-2   |   M4-3   |   M4-4   |   M4-5   
STEMMING cnt         |   6.29e+11|  8.60e+11|  3.00e+08|  8.68e+12|  2.87e+10|  1.34e+11|  3.52e+11|  1.16e+09
STEMMING runtime (s) |      10.40|     10.31|      4.42|     11.88|     13.01|     12.49|     12.45|     12.85
```
