# FDASS
Functional Domain Annotation by Structural Similarity


We aim to annotate the functional domains of T. brucei based on the
structural similarity of its proteins to annotated domains. In this 
regard, we gathered a database of structure of Pfam instances. It 
was done by cutting the structure of Proteins used as Pfam instance
on their Pfam borders. Structure database of Pfam instances can be 
downloaded from shorturl.at/noXZ0.


We searched proteome of T. brucei and three other organisms against 
the Pfam instances by foldseek using following commands:


```
foldseek search {OrgDB} {PfamDB} aln_{Org}_pf tmpFolder -a --max-seqs 10000000 --cov-mode 1 -c 0.8 -e 3
foldseek convertalis {OrgDB} {PfamDB} aln_{Org}_pf aln_{Org}_pf_e3.tsv --format-output \
"query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,qlen,tlen,evalue,bits,alntmscore,lddt"
```


In the notebooks, we have shown that based on alignment characteristics of 
query-pfam instance we can train a model to predict if the aligned part of 
the query has the same pfam as the instance it has aligned to. We have used
the trained model to predict new domains in T. brucei.


To reproduce the data, first run DownloadFilesAndMakeDirs.sh to download all
the data needed for the next steps. AliLabeler.py and AliLabelerSeq.py label
the alignments. TrainingAndBenchmarking.ipynb trains a model and benchmarks 
the trained model.

Preprint available at: 

https://www.biorxiv.org/content/10.1101/2023.01.18.524644v1
