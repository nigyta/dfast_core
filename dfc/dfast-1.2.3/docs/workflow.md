# DFAST workflow


## Input
DFAST accepts a FASTA-formatted genome sequence file as an input.  
* `GENOME_FASTA`
  Genomic FASTA file. Multi-FASTA files are also acceptable.
  Normally, this is given by a command-line option `--genome` or `-g`. You can specify this in the configuration file.

The following options specified by `GENOME_CONFIG` in the configuration file are available.  
- `use_original_name` (default: False)  
  If set to True, the first word in the Fasta header line will be used as a sequence name.  
  Otherwise, sequences will be renamed as 'Sequence001, Sequence002, ...'.
- `sort_sequence` (default: True)  
  By default, sequences will be sorted so that longer sequences come first.
  If set to False, the original order in the FASTA file will be retained.
- `minimum_length` (default: 200)  
  Sequences shorter than this value will be eliminated. 
  GenBank requires a minimum length of 200 base pairs for the submission.
- `complete` (default: False)  
  If set to True, the input genome will be treated as a complete genome. 
  Also, `sort_sequence` is set to False, and `minimum_length` is set to 0.  
  In the current version, changing this flag does not affect the annotation result.
  It is only used when creating submission files to DDBJ and GenBank,
  thus you may leave this as default unless you need submission files.
  
You can override the values specified in the configuratin file by providing command line options: 
`--use_original_name`, `--sort_sequence`, `--minimum_length`, `--complete`.  


## Structural annotation

DFAST first detects biological features such as CDS, rRNA, and tRNA. 
DFAST calls external programs for this, and they are executed in parallel.  
Currently, programs listed below are incorporated. 
The programs marked with * asterisk are used in the default workflow and 
their executables are bundled in the DFAST distribution.  

The workflow is defined by the `STRUCTURAL_ANNOTATION` attribute in the configuration file,
in which options for each setting is specified as a dictionary. 
You can switch enabled/disabled and set options passed to each prorgam.

1. MGA* (MetaGeneAnnotator)  
CDS prediction tool.  
http://metagene.nig.ac.jp
2. Prodigal  
CDS prediction tool.  
http://prodigal.ornl.gov
3. Barrnap*  
rRNA prediction tool.  
http://www.vicbioinformatics.com/software.barrnap.shtml
4. RNAmmer  
rRNA prediction tool. RNAmmer requires the hmmscan program version 2.3, which is not included in the DFAST distibution. 
To install RNAmmer, follow the [instruction](https://blog.karinlag.no/2013/10/rnammer-install/) by original authors.  
http://www.cbs.dtu.dk/services/RNAmmer/
5. Aragorn*  
tRNA prediction tool.  
http://mbio-serv2.mbioekol.lu.se/ARAGORN/Downloads/
6. tRNAscan-SE  
tRNA prediction tool.  
http://eddylab.org/software.html
7. CRT* (Crisper Recognition Tool)  
CRISPR prediction tool.
CRT requires Java Runtime Environment (JRE). Type "java" in your console to check if Java is installed.  
Depending on your system, you may need to configure options to launch Java.
If not working properly, set `java_options` in the configuration file to "-Xmx256m", for example. 
You can also export _JAVA_OPTIONS to environmental variables before running DFAST.  
e.g. `export _JAVA_OPTIONS="-Xmx256m -XX:ParallelGCThreads=1`  
http://www.room220.com/crt/
8. Gap annotator*  
Gap annotator is not an external program, but works as other programs do.
It finds gap regions (runs of Ns or ns) in the sequence, 
and annotate them as assembly_gap features.  
By default `len_cutoff` is set to 5, which means gaps shorter than this value will be ignored.
`linkage_evidence` and `gap_type` are assigned with the values specified in the configuration file.  
(default: "paired-ends" and "within scaffold")  
Please see this [page](https://www.ncbi.nlm.nih.gov/genbank/wgs_gapped/) to submit gapped sequences to the INSDC.

If you want to use programs not included in the DFAST distribution, 
you need to install them and put the executables or links in your PATH. 
You can also put them in $DFAST_APP_ROOT/bin/Darwin or $DFAST_APP_ROOT/bin/Linux,
which will be automatically added to the PATH variable when running DFAST.

## Feature optimization
After structural annotation, DFAST will cleanup unwanted features. 
This includes three methods below. 
You can alter default settings by specifying `FEATURE_ADJUSTMENT` in the configuration file.
1. Remove partial features  
Remove partial features predicted at the end of a sequence or abutting a gap.  
https://www.ncbi.nlm.nih.gov/genbank/genomesubmit_annotation/#partial_CDS  
2. Remove overlapping features  
Remove features that overlaps other features according to `feature_type_priority`.  
By default, it is defined as following:  
`["assembly_gap", "CRISPR", ("tmRNA", "tRNA", "rRNA"), "CDS"]`  
"assembly_gap" features have the highest priority, and CDS features has the lowest.
Thus, CDSs that overlap other types of features are removed. 
Features enclosed in parentheses have a same priority level.
3. Merge CDS (preliminary implementation)  
Current version is a preliminary implementation, and is disabled by default.  
This method merges CDS features predicted by different programs.
When two CDSs from different programs conflict (usually with different start-codon positions),
the one from the program with higher priority will be adopted.
## Functional Annotation
The main purpose of the functional annotation process is 
to infer protein function for each CDS feature by searching against various types of reference databases.

The workflow is defined by the `FUNCTIONAL_ANNOTATION` attribute in the configuration file.
`FUNCTIONAL_ANNOTATION` is a list of dictionaries (associative array or hash), 
and each dictionary represents options for each annotation process. 
Annotation process will be executed in the order described in the list.  

All annotation process have a `skipAnnotatedFeature` attribute. 
If set to True, CDS features whose function is already assigned in the earlier stage 
will be excluded in the process. This will enable a hierarchial approach 
similar to Prokka, and reduce running time.  

* DBsearch (Database search)  
This is a method for a conventional homology search against reference databases.  
By default, it uses GHOSTX as a sequence aligner, which is about 100 times more efficient than BLAST. 
You can choose an aligner from GHOSTX/BLASTP.  
You can specify thresholds for E-value (`evalue_cutoff`), 
query coverage % (`qcov_cutoff`), and subject coverage % (`scov_cutoff`).  
The reference database is specified by `database`. 
The DFAST-format reference file is a tab-separated table consisting of 7 columns,
which is detailed in the **Reference Database** section.  
You can prepare your own database following this format. 
Before running the pipeline, index files for BLASTP and GHOSTX must be created.
This can be done as following:
  ```
  python $DFAST_APP_ROOT/scripts/ref_util.py format your_database.ref
  ```
  DBsearch expects a small-sized curated database, as it loads all data in the reference file into memory.
  If you want to search against a large database use BlastSearch instead.

* BlastSearch  
This is meant for protein homology search against a large-sized reference database,
such as pre-formatted Blast databases like RefSeq Protein and SwissProt available at the NCBI FTP site.  
BlastSearch uses "blastp" and "blastdbcmd" executables, which are bundled in the DFAST distribution.  
DFAST can parse NCBI- and SwissProt-style Fasta descriptions.
If you want to use your custom-made database, 
format it using the `makeblastdb` command with `-parse_seqids` option.

* OrthoSearch (ortholog assignment based on Reciprocal-Best-Hit search)  
OrthoSearch identifies orthologous genes based on a simple Reciprocal-Best-Hit (RBH) approach.
Usually, it is executed in the first step of the functional annotation workflow.  
It first conducts all-against-all pairwise protein alignments between a query genome and each of reference genomes.  
When RBH relationship is found between a pair of genes, they are considered as orthologs.  
It also conducts self-to-self alignments within a query genome to find in-paralogs (co-orthologs), 
which arose from a duplication event that postdates the speciation event, 
and hence can be considered to share the same function.
When a gene in a query genomes possesses genes with higher similarity scores than its correspongind orthologous gene,
 those genes are considered as in-paralogs and the same protein function is assigned.  
OrthoSearch is effective in transferring annotation from closely related organisms and
in reducing the running time especially when using Blastp as an aligner.  

  To enable this, specify paths to the reference file in `references` in the configuration 
  or specify them with the `--references` option.
  The recommended format of the reference file is a GenBank flat file format containing annotated CDS features.
  We have a helper script to download a GenBank-formatted file from the Assembly Database of NCBI.
  ```
  python $DFAST_APP_ROOT/scripts/file_downloader.py --assembly GCF_000091005.1 GCA_000008865.1
  ```
  A FASTA-formatted file containing all protein sequences in a genome and a DFAST-format reference file are also acceptable as a reference. The file format is automatically recognized.


* CDDsearch  
CDDsearch uses RPS-Blast and rpsbproc (post-processing utility for RPS-Blast) 
to search against the [Conserved Domain Database](https://www.ncbi.nlm.nih.gov/cdd/),
which is a collection of domain models curated by NCBI. roS 
Precompiled databases are available at the [CDD FTP site](ftp://ftp.ncbi.nih.gov//pub/mmdb/cdd/little_endian),
which includes CDD, Pfam, Smart, COG, PRK, TIGRFAM, and Kog.
You can download CDD databases using the `prepare_database.py` script:
  ```
  python $DFAST_APP_ROOT/scripts/file_downloader.py --cdd Cdd Cog Pfam
  ```
  DFAST standard pipeline includes assignment of COG functional categories by CDDsearch.
  DFAST comes with the binaries for RPS-Blast and rpsbproc.

  To learn more, please refer to the READMEs for each program. [[CDD](ftp://ftp.ncbi.nih.gov//pub/mmdb/cdd/README), [rpsbproc](ftp://ftp.ncbi.nih.gov//pub/mmdb/cdd/rpsbproc/README)]

* HMMscan  
HMMscan searches for domain structures in query sequences using profile HMM databases. 
It uses hmmscan of the HMMer software package (ver. 3.1b).  
You can search against various kinds of publicly available HMM databases, 
such as [PFAM](ftp://ftp.ebi.ac.uk//pub/databases/Pfam/) and [TIGRFAM](ftp://ftp.tigr.org//pub/data/TIGRFAMs/TIGRFAMs_15.0_HMM.tar.gz).
After downloading the HMM profile data, it must be formatted for hmmscan by using hmmpress.
You can do it as following:
  ```
  python $DFAST_APP_ROOT/scripts/file_downloader.py --hmm Pfam TIGR
  ```
  DFAST standard pipeline includes HMMscan against TIGRFAM.



* PseudoGeneDetection  
PseudoGeneDetection finds internal stop codons or insertion/deletions that result in frameshifts.  
The process relies on its upstream annotation results, 
and thus should be placed as the final process of Functional Annotation.
It first lists candidates of pseudogenes, 
whose s_cov (subject coverage) in a protein alignment result is less than `scov_cutoff` (default: 85%).
Then their coding sequences are extended to flanking regions by the length specified by `extension`,
and re-aligned to their subject protein sequences using [LAST](http://last.cbrc.jp), which allows frameshift alignments.
When stop codons or frameshifts are found in the extended regions,
the query will be marked as possible pseudogenes.  
This also detects translation exceptions to selenocysteine/pyrrolysine.

  Kiełbasa SM, Wan R, Sato K, Horton P, Frith MC (2011)  
  Adaptive seeds tame genomic sequence comparison.  
  Genome Res 21: 487–493.

## Output
The output directory is defined by `WORK_DIR` in the config file. You can override this by the command line option, `--out` or `-o`.

The DFAST standard pipeline generates the following files:  
* genome.gff, genome.gbk  
  Genome sequence and annotation data in a standard GFF3 and GenBank format, respectively.
* genome.fna, cds.fna, rna.fna,   
  FASTA files for predicted CDSs, rRNAs/tRNAs, and translated protein sequences. 
* genome.fasta  
  Nucleotide FASTA file of the query genome 
* statistics.txt  
  Statistics for genome sequence and annotated features
* ddbj/\*\*\*.ann, ddbj/\*\*\*.fasta  
  Annotation file (\*\*\*.ann) and sequence file \*\*\*fasta) for the DDBJ Mass Submission System.
* genbank/\*\*\*.tbl, ddbj/\*\*\*.fsa  
  Annotation tbl file (\*\*\*.tbl) and sequence file (\*\*\*.fsa) for GenBank submission,
  which can be used as inputs for [tbl2asn](https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/) program.   

#### Output verbosity
DFAST cab set the verbosity level of output files in the range of 1 (simple) to 3 (rich).
In the default configuration, the verbosity level is set to 3 for GFF3 and GenBank files. For DDBJ/GenBank submission files, it is set to 1. You may increase this to 2 to include more information. 

#### For developpers
Optionally, DFAST can dump its native data object containing as much information as possible \[Advanced option for Python developpers\].
* genome.pickle

To generate a pickle file, run DFAST in the debug mode. Use `--debug` option or set `DEBUG` to True in the config file.

You can load and reuse this using Python's pickle module.
```
import pickle
import sys
sys.path.append("$DFAST_APP_ROOT")   # This enables importing modules from the dcore package.

genome = pickle.load(open('genome.pickle', 'rb'))
# genome is a master data object of DFAST that includes all sequence and annotation data.

# to itertate through all the sequences,
for seq_id, record in genome.seq_records.items():
      print(record.format('fasta')   # record is a BioPython's SeqRecord object.

# to walk around annotated features,
for feature_id, feature in genome.features.items():
      # feature is an extended SeqFeature object with additional attibutes
      # seq_id, annotations, ,primary_hit, and secondary_hits.
      print(feature_id, type(feature))
      print(feature.seq_id, feature.annotations, feature.primary_hit)

```

## Reference Databases

DFAST reference protein format is a tab-separated table with 8 columns.
Lines that start with "#" is considered comments.
The first line contains optional attributes of the database, which will be displayed in the log file when the database is loaded. The format should be "[attribute=value]".

```
# [dbname=DFAST-default] [version=1.0] [contributor=yt]
```

The second line and after contain protein sequence data in each line.

|# id |description |gene |EC_number |flag |organism |source_DB |sequence |
|:---|:---|:---|:---|:---|:---|:---|:---|
|O66428 |elongation factor G |fusA  | | |Aquifex aeolicus VF5 |UniProtKB |MAREVPIEKL... |
|WP_001015721.1	|aerobactin synthase IucC |	iucC |6.3.2.27 | | Escherichia coli O83:H1 str. NRG 857C |RefSeq |MNHKDWDLVN...|
|AAC73225.1 |pyruvate dehydrogenase |aceE |1.2.4.1 | | |INSD |MSERFPNDVD...|

* id  
  Sequecne ID (Mandatory. Must be unique.)
* description  
  Protein name (Mandatory)
* gene  
  Gene symbol.
* EC_number  
  EC number. For bifunctional proeins, use comma. e.g. 1.1.1.81, 1.1.1.79
* flag  
  This corresponds to the Flags field of the UniProtKB entry to indicate precursor or fragment. Current version of DFAST does not use this.
* organism  
  Organism name from which the sequence derived.
* source_DB  
  The database name from which the sequence is obtained. This will appear in the 'inference' qualifier of a CDS feature. e.g. /inference="DESCRIPTION:similar to AAsequence:UniProtKB:O66428"
* sequence (Mandatory)  
  Amino acid sequence for the protein

#### Utility script
reference_util.py is a utility script for format conversion and database indexing.
To show availbale commands,
```
python $DFAST_APP_ROOT/scripts/reference_util.py -h
```
To show help,
```
python $DFAST_APP_ROOT/scripts/reference_util.py fasta2dfast -h
```
