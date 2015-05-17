# Tutorial_10
HUMAnN: The HMP Unified Metabolic Analysis Network

#Resouces
https://bitbucket.org/biobakery/humann

https://bitbucket.org/biobakery/biobakery/wiki/humann

#Overview

1.Place one or more translated BLAST results using KEGG gene identifiers in the "input" directory (optionally gzipped or bzipped).

2.Edit the "SConstruct" file; in particular, make sure that the input processors include one configured for your BLAST file name(s) and format(s).

3.Run the "scons" command, optionally parallelizing multiple analyses using the "-j" flag. Results will be placed in the "output" directory.

#Step 1: PREREQUISITES  

1.A network connection. At least for the first run.

2.A bunch of RAM. At least 8-10GB

3.Python >= 2.7

4."-outfmt 6" blastx result. Alternatively, input processors are also provided for accelerated BLAST implementations such as mapx, mblastx, or usearch.

#Step 2: Input

1.One or more of the following:

(1) Tabular translated BLAST (blastx) output files files matching sequence read IDs to gene IDs.

(2)Mapping (bowtie, bwa, etc.) output in BAM format.

(3)Tab-delimited text containing one or more pre-quantified gene abundances.

Place (or symlink) each file with a .txt, .txt.gz, .txt.bz, .bam, .csv, .tsv, or .pcl extension as appropriate in the "input" directory before running HUMAnN. The pipeline includes processors for three tab-delimited text formats by default (below) and BAM binary format and can easily by modified to accept more.

Example blast command:
```bash 
blastx -outfmt 6 -db 28_kegg_genomes < mock_even_lc.fasta.gz | gzip -c > mock_even_lc.txt.gz
```
Now kegg has gone commercial, we can use COG, NOG database insdead. Replace 28_kegg_genomes by your own database.

You can build your COG, NOG database by:

```bash 
wget ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/prot2003-2014.fa.gz

gzip -d prot2003-2014.fa.gz

makeblastdb -in *.fa -dbtype prot -out COG -parse_seqids
```

The above example is for COG database. For NOG database, you can find more detail from here
```bash 
http://eggnogdb.embl.de/download/eggnog_4.1/
```

Note that eggNOG data base contain NOG database so we can use eggNOG insdead of NOG.











