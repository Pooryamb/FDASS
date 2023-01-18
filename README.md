# FDASS
Functional Domain Annotation by Structural Similarity



By running the DownloadFilesAndMakeDirs.sh script, all the data needed for running the notebooks will be downloaded. Besides, it will download the databases of PfamSeeds (Both FoldSeek and MMseqs).

In the rawinput directory, there are four directories each for an organism used in the study (Tb, Ec, Mj, and Sc). In each directory, Pfams v35.0 predictions are available named Pfam{Org}.txt. The output of FoldSeek and MMseqs when the organism’s proteome was searched against Pfam Instances is found in aln_{Org}_pf_e3.tsv and aln_{Org}_pf_seq_e3.tsv respectively. 
The following commands have generated the FoldSeek output:

```
foldseek search {OrgDB} {PfamDB} aln_{Org}_pf tmpFolder -a --max-seqs 10000000 --cov-mode 1 -c 0.8 -e 3
foldseek convertalis {OrgDB} {PfamDB} aln_{Org}_pf aln_{Org}_pf_e3.tsv --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,qlen,tlen,evalue,bits,alntmscore,lddt"
```

And the following commands generate the MMseqs output:
```
mmseqs search {OrgDB} {PfamDB} aln_{Org}_pf_seq tmpFolder -a --max-seqs 10000000 --cov-mode 1 -c 0.8 -e 3
mmseqs convertalis {OrgDB} {PfamDB} aln_{Org}_pf aln_{Org}_pf_seq_e3.tsv --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,qlen,tlen,evalue,bits"
```
By running AliLabeler.py and AliLabelerSeq.py, the output of FoldSeek and MMseqs will be labeled, which will later be used for training and benchmarking.

