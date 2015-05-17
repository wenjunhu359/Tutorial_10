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

#Step 2: Put your inputs into the input folder /humann-0.99/input, edit the SConstruct.

One or more of the following:

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
(blast module is needed)

The above example is for COG database. For NOG database, you can find it from here
```bash 
http://eggnogdb.embl.de/download/eggnog_4.1/
```

Note that eggNOG database contains NOG database.

And then we should modify the SConstruct file to specify the exact format of our input data:
```bash
blastx -outfmt 6
----------------
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

mapx
----
template-name frame read-name template-start template-end template-length read-start read-end read-length template-protein read-protein alignment identical %identical positive %positive mismatches raw-score bit-score e-score

mblastx
-------
Query_Id Reference_Id E_Value_Bit_Score Identity_Percentage_Length_of_HSP Number_of_Positives_Frame_#s Alignment_Start_in_Query Alignment_End_in_Query Alignment_Start_in_Reference Alignment_End_in_Reference

BAM (bowtie/bwa/etc.)
--------------------
qname flag rname pos mapq cigar rnext pnext tlen seq qual

usearch
-------
usearch can be treated identically to "blastx -outfmt 6" if the "--blast6out" flag is provided.
```

#Step 3: Using alternatives to KEGG (optional)

The process of using HUMAnN with a database other than KEGG (e.g. COG, NOG, etc.) requires:
```bash
1.A FASTA file of nucleotide sequences against which the meta'ome is searched, each labeled with a gene ID (for KEGG, the genes.pep file distributed with KEGG).
2.A file of nucleotide sequence lengths, each labeled with a gene ID (for KEGG, data/genels).
3.A file of gene-to-OG mappings (for KEGG, data/koc).
4.A file of OG-to-pathway mappings (for KEGG, data/keggc).
```
Note that the newest version of HuManN use both KEGG and Metacyc as default. After we use COG or eggNOG to blast, the HuManN would check both. We don't have to modified the code if you willing to use KEGG and MetaCyc.

There is no clue that which one is better. The KEGG vs. MetaCyc will be presented in the latter sildes.

#Step 4: Cd to the humann-0.99 folder, run ./$ scons

#Step 5: Check the output

Three or more files per input including:
```bash
The relative abundances of each gene in the input metagenome. By default tagged as type "01". Two columns of tab-delimited text: geneid abundance.
The coverages of each pathway, expressed as a fraction between 0 and 1 inclusive. By default tagged as type "04a". Two columns of tab-delimited text: pathid coverage.
The relative abundances of each pathway. By default tagged as type "04b". Two columns of tab-delimited text: pathid abundance.
Optionally, a table of individual gene abundances appropriate for loading into METAREP. By default tagged as type "99". Five columns of tab-delimited text: geneid abundance e-score %identical identical. The abundance is relative and calculated as in HUMAnN's gene family abundances; e-score, percent identity, and identity length are averaged over all reads mapping to each gene in the input translated BLAST results.
```
Three or more merged files are also produced, all tab-delimited text:
```bash
A table containing the relative abundances of all genes in any input metagenomes. By default named as "01-*.txt".
A table containing the coverages of all pathways in any input files. By default named as "04a-*.txt".
A table containing the relative abundances of all pathways in any input files. By default named as "04b-*.txt".
A table containing the relative abundances of all individual genes in any input files. By default named as "99-*.txt".
```
#Step 6 Further analysis

GraPhlAn:

The following HUMAnN output files may be used as inputs for visualization with GraPhlAn:
04b-*-(mpm or mpt)-*-graphlan_rings.txt (Annotation file)
04b-*-(mpm or mpt)-*-graphlan_tree.txt (Tree)

Go to http://huttenhower.sph.harvard.edu/galaxy/

Click the "GraPhlAn" link on the left and follow the instructions

MaAsLin:

MaAsLin is a multivariate statistical framework that finds associations between clinical metadata and potentially high-dimensional experimental data. 

Modifications:
```bash
1. Select a file from the HUMAnN output folder (named 04b-*-mpt-*.txt or 04b-*-mpm-*.txt)
2. Open the file in Microsoft Excel or a text editor.
3. Remove the first column.
```

Tutorial of MaAsLin: https://bitbucket.org/biobakery/biobakery/wiki/maaslin

LEfSe:

LDA Effect Size is an algorithm for high-dimensional biomarker discovery and explanation that identifies genomic features (genes, pathways, or taxa) characterizing the differences between two or more biological conditions (or classes).

Modifications:
```bash
1.Select a file from the HUMAnN output folder (named 04b-*-mpt-*.txt or 04b-*-mpm-*.txt)
2.Open the file in Microsoft Excel or a text editor.
3.Remove the first column.
4.Remove every metadata row (anything including and above InverseSimpson) except the class (and optional subclass), and the top row: ID/NAME.
5.Please ensure only 1-2 metadata rows remain apart from the Name/ID row at the top.
```

Go to http://huttenhower.sph.harvard.edu/galaxy/

Click the "LEfSe" link on the left and follow the instructions

















