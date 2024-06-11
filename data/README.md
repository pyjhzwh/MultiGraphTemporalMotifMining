## Downloading the datasets

- Get dataset from [https://snap.stanford.edu/temporal-motifs/data.html](https://snap.stanford.edu/temporal-motifs/data.html) and [https://www.cs.cornell.edu/~arb/data/](https://www.cs.cornell.edu/~arb/data/)
- Rename it to *.el and move to examples/non_attributed/
- Run `python el_to_gdf_converter.py -i .el`


Example: 
```
wget https://snap.stanford.edu/data/wiki-talk-temporal.txt.gz
gunzip wiki-talk-temporal.txt.gz
mv wiki-talk-temporal.txt wiki-talk-temporal.el
mv wiki-talk-temporal.el ../examples/non-attributed/
cd ../examples/non-attributed/
python el_to_gdf_converter.py -i wiki-talk-temporal.el
cd ../..
make -j
./graph_search -g examples/non-attributed/wiki-talk-temporal.gdf -q examples/non-attributed/ex1_q.gdf -delta 10000 > run.out
```


These are verified correctness links: 
- wiki-talk (WT): [https://snap.stanford.edu/data/wiki-talk-temporal.txt.gz](https://snap.stanford.edu/data/wiki-talk-temporal.txt.gz)
- stackoverflow (SO) [https://snap.stanford.edu/data/sx-stackoverflow.txt.gz](https://snap.stanford.edu/data/sx-stackoverflow.txt.gz)
- bitcoin (BI) [https://drive.google.com/open?id=1XuLGdDNshvASMyILQb41GKqpAgQZ_sET](https://drive.google.com/open?id=1XuLGdDNshvASMyILQb41GKqpAgQZ_sET)
- reddit-reply (RE) [https://drive.google.com/open?id=1JCiDApvzOml1adq1tKUiWe8KsBmqarbM](https://drive.google.com/open?id=1JCiDApvzOml1adq1tKUiWe8KsBmqarbM)
