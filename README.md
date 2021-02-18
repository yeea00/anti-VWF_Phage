## To run it

1. download blastn 2.11.0+ from https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.11.0/

2. after unzipping the blast file, add the program to your PATH. If it is Linux, and the file was unzipped under /usr/local, the command will belike
```
	$ export PATH=$PATH:/usr/local/ncbi-blast-2.11.0+/bin
```

3. A command example
```
	$ ./BLAST_frag.py -v -q VWF_phage.fa -k 500 -o count.txt -p 97 -w 11 -e 20 *_R1.fq.gz *_R2.fq.gz
```
