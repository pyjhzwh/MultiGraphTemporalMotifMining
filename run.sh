#!/bin/bash
make -j4
MOTIF="73"
GRAPH="wiki-talk-temporal"
THREAD=32
MAX_TRIAL=1e8
SPANNING_TREE=108
for GRAPH in "CollegeMsg" "wiki-talk-temporal" "sx-stackoverflow" #"bitcoin-temporal" "temporal-reddit-reply"
do
    echo ""
    for MOTIF in "73"
    do
    # echo "GRAPH="$GRAPH "MOTIF="$MOTIF
    for delta in "1D" "1W" "4W" "8W" #"1M" "1H" "1D" "1W" "4W" "8W" "16W"
        do
            echo "GRAPH="$GRAPH "MOTIF="$MOTIF "delta="$delta
            ./graph_search -g ../TEACUPS/examples/non-attributed/${GRAPH}.gdf -q examples/non-attributed/${MOTIF}.gdf -delta $delta -algo 0 -max_trial ${MAX_TRIAL} -thread ${THREAD} -spanning_tree ${SPANNING_TREE} | tail -2
        done
    done
done
