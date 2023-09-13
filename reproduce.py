import subprocess
import shlex
import numpy as np
import json

GRAPH_DIR = "examples/non-attributed" # directory of input graphs
MOTIF_DIR = "examples/non-attributed" # directory of motifs
# input graphs
graphs = ["wiki-talk-temporal", "sx-stackoverflow", "bitcoin-temporal", "temporal-reddit-reply"]
# delta constraints for each input graph
deltas = [["4W", "8W", "16W"], ["4W", "8W", "16W"], ["1D", "1W"], ["1D", "1W"]]
# M0-M7
motifs = [f"M{i}" for i in range(8)]
# number of samples to take in TEACUPS for each input graph
max_trials = ["2e7", "1e8", "1e8", "1e8"]
threads = 32

# whether to run TEACUPS or not, default is true
RUN_TEACUPS = True
# whether to run backtracking baseline or not, default is false because it's very slow (take days to weeks)
RUN_BT = False 


for i, graph in enumerate(graphs):
    for delta in deltas[i]:
        print(f"graph={graph}, delta={delta}")
        if RUN_TEACUPS:
            # for motif in motifs[]:
            command = f"./graph_search -g {GRAPH_DIR}/{graph}.gdf -delta {delta} -algo 1 -thread {threads} -max_trial {max_trials[i]}"
            # run it sequantially
            pipe = subprocess.run(shlex.split(command), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            trialn = 2 # matches and runtime for baseline and sampling
            pipe2 = subprocess.run(["tail", f"-{trialn}",], input=pipe.stdout, stdout=subprocess.PIPE)
            pipe2_lines = pipe2.stdout.splitlines()

            line = pipe2_lines[0].decode().strip()
            cnts = line.split(",")[:-1]
            cnts = [float(cnt.split(":")[-1]) for cnt in cnts]
            
            line = pipe2_lines[1].decode().strip()
            M1_M7_runtime = float((line.split(":")[-1]).split("ms")[0]) / 1e3
            
            command = f"./graph_search -g {GRAPH_DIR}/{graph}.gdf -delta {delta} -algo 2 -thread {threads}"
            pipe = subprocess.run(shlex.split(command), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            trialn = 2 # matches and runtime for baseline and sampling
            pipe2 = subprocess.run(["tail", f"-{trialn}",], input=pipe.stdout, stdout=subprocess.PIPE)
            pipe2_lines = pipe2.stdout.splitlines()
            
            line = pipe2_lines[0].decode().strip()
            M0_cnt = float(line.split(":")[-1])
            cnts = [M0_cnt, *cnts]
            
            line = pipe2_lines[1].decode().strip()
            M0_runtime = float((line.split(":")[-1]).split("ms")[0]) / 1e3
            
            print("                    |", f"{'|'.join(format(motif, '^10s') for motif in motifs)}")
            print("TEACUPS cnt         |", f"{'|'.join(format(cnt, '10.2e') for cnt in cnts)}")
            print("TEACUPS runtime (s) |", f"{M0_runtime:10.2f}|", "        "*3, f"{M1_M7_runtime:10.2f}")
        
        if RUN_BT:
            BT_cnts = []
            BT_runtime = []
            for motif in motifs:
                command = f"./graph_search -g {GRAPH_DIR}/{graph}.gdf -q {MOTIF_DIR}/{motif}.gdf -delta {delta} -algo 0 -thread {threads}"
                pipe = subprocess.run(shlex.split(command), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

                trialn = 2 # matches and runtime for baseline and sampling
                pipe2 = subprocess.run(["tail", f"-{trialn}",], input=pipe.stdout, stdout=subprocess.PIPE)
                pipe2_lines = pipe2.stdout.splitlines()
                
                line = pipe2_lines[0].decode().strip()
                Mi_cnt = float(line.split(":")[-1])
                BT_cnts.append(Mi_cnt)
                
                line = pipe2_lines[1].decode().strip()
                Mi_runtime = float((line.split(":")[-1]).split("ms")[0]) / 1e3
                BT_runtime.append(Mi_runtime)
                
            print("BT      cnt         |", f"{'|'.join(format(cnt, '10.2e') for cnt in BT_cnts)}")
            print("BT      runtime (s) |", f"{'|'.join(format(runtime, '10.2f') for runtime in BT_runtime)}")