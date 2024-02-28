import subprocess
import shlex
import numpy as np
import json

GRAPH_DIR = "examples/non-attributed" # directory of input graphs
MOTIF_DIR = "examples/non-attributed" # directory of motifs
# input graphs
graphs = ["wiki-talk-temporal", "sx-stackoverflow", "bitcoin-temporal", "temporal-reddit-reply"]
# delta constraints for each input graph
deltas = [["4W", "8W"], ["8W", "16W"], ["1D"], ["1H", "1D"]]
motifs = ["M3-1", "M3-2", "M4-1", "M4-2", "M4-3", "M4-4"]
# number of samples to take in STEMMING for each input graph
max_trials = "1e8"
threads = 32

# whether to run STEMMING or not, default is true
RUN_STEMMING = True
# whether to run backtracking + derivecnts or not, default is false
RUN_STEMMING_BT = True
# whether to run backtracking baseline or not, default is false because it's very slow (take days to weeks)
RUN_BT = True 


for i, graph in enumerate(graphs):
    for delta in deltas[i]:
        print(f"graph={graph}, delta={delta}")
        if RUN_STEMMING:
            cnts = []
            runtime = []
            for motif in motifs:
                if (graph == "sx-stackoverflow" and motif == "M4-1" and delta == "8W") or \
                    (graph == "temporal-reddit-reply" and motif == "M4-1" and delta == "1D") or \
                    (graph == "sx-stackoverflow" and motif == "M4-4" and delta == "16W"):
                    max_trials = "1e9"
                elif graph == "temporal-reddit-reply" and motif == "M3-2" and delta == "1D":
                    max_trials = "5e8"
                else:
                    max_trials = "1e8"
                command = f"./graph_search -g {GRAPH_DIR}/{graph}.gdf -q {MOTIF_DIR}/{motif}.gdf -delta {delta} -algo 1 -thread {threads} -max_trial {max_trials}"
                # run it sequantially
                pipe = subprocess.run(shlex.split(command), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

                trialn = 2 # matches and runtime for baseline and sampling
                pipe2 = subprocess.run(["tail", f"-{trialn}",], input=pipe.stdout, stdout=subprocess.PIPE)
                pipe2_lines = pipe2.stdout.splitlines()

                line = pipe2_lines[0].decode().strip()
                Mi_cnt = float(line.split(":")[-1])
                cnts.append(Mi_cnt)
                
                line = pipe2_lines[1].decode().strip()
                Mi_runtime = float((line.split(":")[-1]).split("s")[0])
                runtime.append(Mi_runtime)
            
            print("                     |", f"{'|'.join(format(motif, '^10s') for motif in motifs)}")
            print("STEMMING cnt         |", f"{'|'.join(format(cnt, '10.2e') for cnt in cnts)}")
            print("STEMMING runtime (s) |", f"{'|'.join(format(time, '10.2f') for time in runtime)}")
            
        if RUN_STEMMING_BT: # backtracking + derivecnts
            cnts = []
            runtime = []
            for motif in motifs:
                # for motif in motifs[]:
                command = f"./graph_search -g {GRAPH_DIR}/{graph}.gdf -q {MOTIF_DIR}/{motif}.gdf -delta {delta} -algo 2 -thread {threads}"
                # run it sequantially
                pipe = subprocess.run(shlex.split(command), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

                trialn = 2 # matches and runtime for baseline and sampling
                pipe2 = subprocess.run(["tail", f"-{trialn}",], input=pipe.stdout, stdout=subprocess.PIPE)
                pipe2_lines = pipe2.stdout.splitlines()

                line = pipe2_lines[0].decode().strip()
                Mi_cnt = float(line.split(":")[-1])
                cnts.append(Mi_cnt)
                
                line = pipe2_lines[1].decode().strip()
                Mi_runtime = float((line.split(":")[-1]).split("s")[0])
                runtime.append(Mi_runtime)
                
            # print("                                   |", f"{'|'.join(format(motif, '^10s') for motif in motifs)}")
            print("BT+derivecnts cnt         |", f"{'|'.join(format(cnt, '10.2e') for cnt in cnts)}")
            print("BT+derivecnts runtime (s) |", f"{'|'.join(format(time, '10.2f') for time in runtime)}")
    
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
                Mi_runtime = float((line.split(":")[-1]).split("s")[0])
                BT_runtime.append(Mi_runtime)
                
            print("BT      cnt         |", f"{'|'.join(format(cnt, '10.2e') for cnt in BT_cnts)}")
            print("BT      runtime (s) |", f"{'|'.join(format(runtime, '10.2f') for runtime in BT_runtime)}")