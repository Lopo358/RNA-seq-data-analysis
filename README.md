# RNA-seq data analysis
Perform differential gene expression analysis starting from .gtf files

## StringTie

StringTie is is a tool for transcript assembly and qunatification for RNA-seq. StringTie outputs .gtf files that contian assembled transcripts (https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual#output).

## prepDE.py extracts count data from the .gtf files

The prepDE.py3 program can be obtained here:
https://github.com/gpertea/stringtie/blob/master/prepDE.py3

I was only able to run the program on one file at a time with directories arranged like this:

```
+---prepDE.py3
|
+---sample1
|    \sample1
|        \sample1.gtf
+---sample2
    \sample2
        \sample2.gtf
```

In the command line, navigate to the top level sample1 directory and execute this line:
```
python ../prepDE.py3 -v
```
This will result in the gene counts for each sample in gene_count_matrix.csv, which can be loaded into R for differential gene expression analysis.
