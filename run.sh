#!/bin/bash
make -j4
MOTIF="73"
GRAPH="wiki-talk-temporal"
THREAD=32
MAX_TRIAL=1e8
for GRAPH in "wiki-talk-temporal" #"sx-stackoverflow" #"bitcoin-temporal" "temporal-reddit-reply"
do
    echo ""
    for MOTIF in "M3-1" #"M3-2" "M4-6" "M4-7" "M4-8" "M4-9"
    do
    # echo "GRAPH="$GRAPH "MOTIF="$MOTIF
    for delta in "1W" #"1M" "1H" "1D" "1W" "4W" "8W" "16W"
        do
            echo "GRAPH="$GRAPH "MOTIF="$MOTIF "delta="$delta
            # ./graph_search -g ../STEMMING/examples/non-attributed/${GRAPH}.gdf -q examples/non-attributed/${MOTIF}.gdf -delta $delta -algo 0 -max_trial ${MAX_TRIAL} -thread ${THREAD} -spanning_tree ${SPANNING_TREE} | tail -2
            # ./graph_search -g ../SpanningTreeSampling/examples/non-attributed/${GRAPH}.gdf -q ./examples/non-attributed/${MOTIF}.gdf -delta $delta -algo 0 -max_trial ${MAX_TRIAL} -thread ${THREAD} | tail -2
            ./graph_search -g ../SpanningTreeSampling/examples/non-attributed/${GRAPH}.gdf -q ./examples/non-attributed/${MOTIF}.gdf -delta $delta -algo 1 -max_trial ${MAX_TRIAL} -thread ${THREAD} #| tail -2
            # ./graph_search -g ../SpanningTreeSampling/examples/non-attributed/${GRAPH}.gdf -q ./examples/non-attributed/${MOTIF}.gdf -delta $delta -algo 2 -max_trial ${MAX_TRIAL} -thread ${THREAD} | tail -2
        done
    done
done
