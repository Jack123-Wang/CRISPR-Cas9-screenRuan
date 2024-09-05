

```
#mageck count
mageck count -l 'Brunello_library_mageck.csv' -n demo --sample-label lib1,lib2,CTRL1,CTRL2  --fastq lib-1.fastq lib-2.fastq Input-B2B-180924.top20.fasta Input-B2B-181105.top20.fasta

#mageck test
mageck test -k demo.count.txt -t lib1 -c CTRL1 -n demo.lib1
mageck test -k demo.count.txt -t lib2 -c CTRL2 -n demo.lib2
```

