#!/bin/bash
make -j4
MOTIF="73"
GRAPH="wiki-talk-temporal"
THREAD=32
MAX_TRIAL=1e9
for GRAPH in "sx-stackoverflow" #"wiki-talk-temporal" "sx-stackoverflow" #"bitcoin-temporal" "temporal-reddit-reply"
do
    echo ""
    for MOTIF in "M4-9"
    do
    # echo "GRAPH="$GRAPH "MOTIF="$MOTIF
    for delta in "16W" #"1M" "1H" "1D" "1W" "4W" "8W" "16W"
        do
            echo "GRAPH="$GRAPH "MOTIF="$MOTIF "delta="$delta
            # ./graph_search -g ../TEACUPS/examples/non-attributed/${GRAPH}.gdf -q examples/non-attributed/${MOTIF}.gdf -delta $delta -algo 0 -max_trial ${MAX_TRIAL} -thread ${THREAD} -spanning_tree ${SPANNING_TREE} | tail -2
            ./graph_search -g ./examples/non-attributed/${GRAPH}.gdf -q examples/non-attributed/${MOTIF}.gdf -delta $delta -algo 2 -max_trial ${MAX_TRIAL} -thread ${THREAD} #| tail -5
        done
    done
done
