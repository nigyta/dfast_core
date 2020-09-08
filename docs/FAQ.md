# FAQ

## 1. How to use Prodigal to predict CDSs
1. Make sure that 'prodigal' is in your `PATH`.  
2. Set the option `--use_prodigal`.  
    When enabled, Prodigal will be used instead of the default prediction tool (MGA).
    
Note that Prodigal is not bundled in the DFAST distribution. Please install it by yourself. 


## 2. How to use GeneMarkS2 to predict CDSs
1. Make sure that 'GeneMarkS2' is in your `PATH`, e.g. create a symbolic link for gms2.pl in '/usr/local/bin'.  
2. Set the option `--use_genemarks2 GENOME-TYPE`.  
	GENOME-TYPE should be 'auto', 'bact(eria)', or 'arch(aea)'.
    When enabled, GeneMarkS2 will be used instead of the default prediction tool (MGA).

GeneMarkS2 is invoked with the command like below:
`gms2.pl  --genome-type bacteria(or auto, archaea) --gcode 11 --format gff --seq input/genome.fna --output output_file.txt`  
Codon table for Mycoplasma can be used by specifying DFAST's commnand line option `--gcode 4`.
  
  
Note that GeneMarkS2 is not bundled in the DFAST distribution. Please install it by yourself. 


## 3. How to use RNAmmer to predict rRNAs
1. Make sure that 'RNAmmer' is in your `PATH`, e.g. create a symbolic link for the 'rnammer' executable in '/usr/local/bin'.  
2. Set the option `--use_rnammer [bact|arch]`.  
    `--use_rnammer bact` for bacterial genome, `--use_rnammer arch` for archaeal genome.  
    When enabled, RNAmmer will be used instead of the default prediction tool (Barnnap).

RNAmmer is invoked with the following command,  
`rnammer -S bac(or arc) -m tsu,lsu,ssu  -gff output_file.txt input/genome.fna`,  
meaning that all kind of rRNA genes will be predicted using the parameter for bacteria(or archaea).
Note that RNAmmer is not bundled in the DFAST distribution. Please install it by yourself. 

## 4. How to use tRNAscan-SE to predict tRNAs
1. Make sure that 'tRNAscan-SE' is in your `PATH`.  
2. Set the option `--use_trnascan [bact|arch]`.  
    `--use_trnascan bact` for bacterial genome, `--use_trnascan arch` for archaeal genome.  
    When enabled, tRNAscan-SE will be used instead of the default prediction tool (Aragorn).
    
Both tRNAscan-SE 1.3 and 2.0 can be used.  
Note that tRNAscan-SE is not bundled in the DFAST distribution. Please install it by yourself. 


## 5. How to use Diamond to align protein sequences
1. Make sure that the binary for Diamond is in your `PATH`.  
2. Set the option `--aligner diamond`.  

If an index file for Diamond (.dmnd) does not exist, DFAST attempts to build it. You can build it manually by `scripts/reference_util.py formatdb-dmnd`.
Note that Diamond is not bundled in the DFAST distribution. Please install it by yourself. 



## 6. What are the meanings of 'note' qualifiers in CDS features?
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
