## Downloading the datasets

- Get dataset from [https://snap.stanford.edu/temporal-motifs/data.html](https://snap.stanford.edu/temporal-motifs/data.html)
- Rename it to *.el and move to examples/non_attributed/
- Run `python el_to_gdf_converter.py -i .el`


Example: 
```
wget https://snap.stanford.edu/data/CollegeMsg.txt.gz
gunzip CollegeMsg.txt.gz
mv CollegeMsg.txt CollegeMsg.el
mv CollegeMsg.el ../examples/non-attributed/
cd ../examples/non-attributed/
python el_to_gdf_converter.py -i CollegeMsg.el
cd ../..
make -j
./graph_search -g examples/non-attributed/CollegeMsg.gdf -q examples/non-attributed/ex1_q.gdf -delta 10000 > run.out
```


These are verified correctness links: 

- email-Eu-core: [https://snap.stanford.edu/data/email-Eu-core-temporal.txt.gz](https://snap.stanford.edu/data/email-Eu-core-temporal.txt.gz)
- CollegeMsg: [https://snap.stanford.edu/data/CollegeMsg.txt.gz](https://snap.stanford.edu/data/CollegeMsg.txt.gz)