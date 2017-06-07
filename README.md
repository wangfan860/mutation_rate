## correlation of gene features and the mutation density (cook-book)

Most following steps were running at bash (ubuntu14.04) or ipython (2.7) environment. To use the repo code:
```
import mutation_rate.mutation_density as md
```

####uncompress the attachment.
```
tar -xvf test.assignment.txt.tar.gz
```
####remove decoy/virus sequences and sort.
```
grep -v ^GL test.assignment.txt | sort -k1,1V -k2,2n > test.assignment.no.decoy.tsv
```
####split patients.
```
awk '{print $NF}' test.assignment.no.decoy.tsv | sort -u | while read line; do grep -w $line test.assignment.no.decoy.tsv >> $line.tsv; done
```
####prepare .bed file from individuals for ucsc genome browser, and modify chrMT to chrM for loading.
```
md.creat_bed('input')
```
####calculate mutation density on 1M base window.
```
md.create_density_csv('ucsc.hg19.fasta.fai', 1000000, 'test.assignment.no.decoy.tsv', 'test.assignment.no.decoy.1M_window.mutation_density.csv')
```
####get gene features from gencode.
```
for i in CDS UTR lincRNA; do grep -w $i gencode.v19.annotation.gtf > $i.gtf; done
grep 'gene_type \"protein_coding\"' gencode.v19.annotation.gtf >> protein_coding.gtf
```
####calculate gene features on 1M base window. e.g:
```
for i in range(22):
    n = str(i + 1)
    print 'processing chr%s.utr.csv' % n
    md.create_gc_csv('chr' + n, 'ucsc.hg19.fasta.fai', 1000000, 'utr.gtf', 'chr'+n+'.utr.csv')
    print 'chr%s is done' % n
```
####install R 3.2.1, bioconductor, SomatiCA to get 10000 base window GC content. need to modify the output format to run next step. (remove first column)
```
library(SomatiCA)
data(GCcontent)
write.csv(GCcontent, file="10000_gccontent.csv")
```
####calculate GC content on 1M base window. e.g:
```
for i in range(22):
    n = str(i + 1)
    print 'processing chr%s.gc.csv' % n
    md.create_gc_csv('chr' + n, 'ucsc.hg19.fasta.fai', 1000000, '10000_gccontent.1.csv', 'chr'+n+'.gc.csv')
    print 'chr%s is done' % n
```
P.S: The first attempt was to calculate GC content by UCSC data. However, the code seemed extremely inefficient to handle such a large file. (UCSC GC data is on 5 bases windows, 6.6G in size.)
