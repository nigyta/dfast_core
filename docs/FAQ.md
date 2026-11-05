# FAQ

## 1. CDS prediction: Prodigal (default) and MGA
Prodigal is the default CDS prediction tool from Ver 1.6 (MGA up to Ver 1.5), and `--use_prodigal` is no longer needed (it still selects Prodigal and logs a warning).

To use MGA (MetaGeneAnnotator) instead:
1. Make sure that 'mga' is in your `PATH`, e.g. install the Bioconda package `metagene_annotator` (linux-64 only).  
2. Set the option `--use_mga`.  
    When enabled, MGA will be used instead of the default prediction tool (Prodigal).

MGA supports the genetic code 11 only, so `--use_mga` stops with an error when `--gcode` is not 11. `--metagenome` runs MGA with `-m`.
Note that MGA is not bundled in the DFAST distribution. Please install it by yourself. 

Notes on Prodigal:
* DFAST stops when the Prodigal GFF has no `partial=` attribute. Some builds (e.g. the Ubuntu 24.04 apt package 2.6.3-6) write a broken GFF; use the Bioconda or Debian trixie build instead.
* Prodigal needs at least 20,000 bp of sequence to train. For a shorter input, DFAST automatically runs it with `-p meta` (genetic code 11 only; DFAST stops with other genetic codes).
* `--metagenome` runs Prodigal with `-p meta`, which supports the genetic code 11 only; DFAST stops when it is combined with another `--gcode`.
* Use `--gcode 4` (Mycoplasma) or `--gcode 25` (SR1/Gracilibacteria) for these genetic codes. Start codons are taken from the translation table.


## 2. How to use GeneMarkS2 to predict CDSs
1. Make sure that 'GeneMarkS2' is in your `PATH`, e.g. create a symbolic link for gms2.pl in '/usr/local/bin'.  
2. Set the option `--use_genemarks2 GENOME-TYPE`.  
	GENOME-TYPE should be 'auto', 'bact(eria)', or 'arch(aea)'.
    When enabled, GeneMarkS2 will be used instead of the default prediction tool (Prodigal).

GeneMarkS2 is invoked with the command like below:
`gms2.pl  --genome-type bacteria(or auto, archaea) --gcode 11 --format gff --seq input/genome.fna --output output_file.txt`  
Codon table for Mycoplasma can be used by specifying DFAST's commnand line option `--gcode 4`.
  
  
Note that GeneMarkS2 is not bundled in the DFAST distribution. Please install it by yourself. 


## 3. How to use tRNAscan-SE to predict tRNAs
1. Make sure that 'tRNAscan-SE' is in your `PATH`.  
2. Set the option `--use_trnascan [bact|arch]`.  
    `--use_trnascan bact` for bacterial genome, `--use_trnascan arch` for archaeal genome.  
    When enabled, tRNAscan-SE will be used instead of the default prediction tool (Aragorn).
    
Both tRNAscan-SE 1.3 and 2.0 can be used.  
Note that tRNAscan-SE is not bundled in the DFAST distribution. Please install it by yourself. 


## 4. How to use Diamond to align protein sequences
1. Make sure that the binary for Diamond is in your `PATH`.  
2. Set the option `--aligner diamond`.  

If an index file for Diamond (.dmnd) does not exist, DFAST attempts to build it. You can build it manually by `scripts/reference_util.py formatdb-dmnd`.
DFAST runs Diamond in the `--very-sensitive` mode. Diamond is not bundled in the DFAST distribution; install it, e.g. `conda install -c bioconda diamond`. 



## 5. What are the meanings of 'note' qualifiers in CDS features?
The note qualifier in GenBank or GFF format shows the result of an alignment to the reference sequence.
```
note="Q890K8 chromosomal replication initiator protein DnaA (Lactobacillus plantarum WCFS1)
[pid:65.2%, q_cov:99.8%, s_cov:99.8%, Eval:3.3e-166]"
```
`pid` and `Eval` represent percentage identity and E-value.  
`q_cov` and `s_cov` represents coverages against the query and subject (target) sequences, i.e. the percentages of the alignment length to the query and subject sequences, respectively.  
When either of `q_cov` or `s_cov` is below 70%, the alignment is marked as 'partial hit'.
When internal stop codons or frameshifts are found, they are also annotated in the note qualifier.  
  
`Q890K8` is an accession number of the reference sequence, and the following is its
description and, in parenthesis, the organism name from which it is derived.  
The accession numbers with 'WP_' come from the NCBI RefSeq database, and there are some from UniProt. This entry comes from UniProtKB [[Q890KB](https://www.uniprot.org/uniprot/Q890K8)].


## 6. How to change the location of the reference data.  
You can specify it with a `--dbroot` option. This may be helpful when running DFAST within a Docker container.  

Download the reference data:  
```
dfast_file_downloader.py --protein dfast --cdd Cog --hmm TIGR --dbroot /dfast_db
```
Invoke DFAST:
```
dfast --genome /data/your_genome.fa --out /data/your_result --dbroot /dfast_db
```
