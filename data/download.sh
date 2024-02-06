#!/bin/bash

DATASET="CollegeMsg"
echo $DATASET

cd "${0%/*}"
wget https://snap.stanford.edu/data/$DATASET.txt.gz
gunzip $DATASET.txt.gz
mv $DATASET.txt $DATASET.el
mv $DATASET.el ../examples/non-attributed/
cd ../examples/non-attributed/
python el_to_gdf_converter.py -i $DATASET.el