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





