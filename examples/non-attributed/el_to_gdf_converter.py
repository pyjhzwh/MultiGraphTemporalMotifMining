import argparse
import sys
import re
import numpy
import os
import math
import glob

parser = argparse.ArgumentParser()
parser.add_argument('-i', required=True)
args = parser.parse_args()

tmp = (args.i).split(".")
output_filename = tmp[0] + ".gdf"
print("output_filename:", output_filename)
out_file = open(output_filename, "w")

src_list = []
dst_list = []
timestamp_list = []


with open(args.i) as fp:
    print("Reading from file", args.i)
    for ln in fp:
        if ',' in ln:
            ln_split = ln.split(',')
        if ' ' in ln:
            ln_split = ln.split(' ')
        # print("ln_split:", ln_split)
        src = ln_split[0]
        dst = ln_split[1]
        if len(ln_split) == 3:
            timestamp = ln_split[2]
        elif len(ln_split) == 4:
            timestamp = ln_split[3]
        else:
            sys.exit("[ERROR] Unknown format!!")
        src_list.append(int(src))
        dst_list.append(int(dst))
        timestamp_list.append(int(timestamp))

if len(src_list) != len(dst_list) or len(src_list) != len(timestamp_list):
    sys.exit("[ERROR] Inconsistent list lengths!!")

max_node_id = max(max(src_list), max(dst_list))

out_file.write("nodedef> name INT\n")

for i in range(0, max_node_id + 1):
    out_file.write(str(i) + "\n")

out_file.write("edgedef> node1 INT, node2 INT, time INT\n")
for i in range(len(src_list)):
    out_file.write(str(src_list[i]) + "," + str(dst_list[i]) + "," + str(timestamp_list[i]) + "\n")

out_file.close()