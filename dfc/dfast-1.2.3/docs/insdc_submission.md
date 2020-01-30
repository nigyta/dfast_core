# INSDC submission
This document is intended for INSDC submitters and describes how to prepare submission files using DFAST. Please check the latest submission guideline at the INSDC partners.



For beginners not familiar with the data submission procedure, we recommend using our [DFAST web service](https://dfast.nig.ac.jp/), which offers you a graphical user interface to create submission files to DDBJ.

#### Submission type
The genome submission to the INSDC can be devided into two categories, __WGS__ mainly for a draft genome and __non-WGS__ mainly for a complete genome. DFAST can generate submission files for both types of categories.
###### WGS (draft genomes)
WGS genomes have three types of submission format depending on the presence of sequencing gap.
* Gapless contigs  
  A standard format of WGS submissoin
* Gapped sequences (scaffolds)  
  A newer style for gapped-sequence submission. Basically, it is same as above excepting that sequencing gaps (runs of Ns) must be annotated using assembly\_gap features.
* Gapless contigs and AGP file  
  An older style for gapped-sequence submissoin. AGP file describes the ordering and orientation of contigs within scaffold.   
  Note that DFAST does not support this type of submission as DFAST does not generate an AGP file.

###### non-WGS (Complete/finished genomes)
In a non-WGS genome, each chromosome is in a single sequence. There still can be gaps within the sequences, and  plasmids still can be in multiple pieces. Each sequence in the genome must be assigned to a chromosome or a plasmid.
DFAST requires additional parameters for non-WGS genomes. See later [section](#For a complete genome).
 
## Source modifier (Source qualifier)
In an INSDC flat file, each sequence has at least, and normally only, one __source feature__ that spans the entire sequence and describes the biological source of the sequence.  
The source must have "organism" and "mol\_type" as mandatory qualifiers, 
and the value of "mol\_type" is always "genomic DNA" for genome sequence data.
For microbial genome submission, it is recommended to use "strain" qualifier to describe the strain from which the sequence derived.
In addition, qualifiers such as "serovar", "type\_material", and "culture\_collection"
should be used if necessary.

These source modifiers can be specified by `GENOME_SOURCE_INFORMATION` in the configuration file,
which includes following attributes. Instead of specifying them in the configuration file, you can set them by providing command line options with the same name as source attributes, i.e. use `--organism` to override the attribute `organism` in the config.

#### Basic option
* organism  
  Scientific name for the organism from which the genome was obtained (mandatory for INSDC submission).
  The value specified in the config file can be overriden by the command line option `--organism`
* strain  
  Strain name (recommended for microbial genome submission)
  The value specified in the config file can be overriden by the command line option `--strain`

#### <a name="For a complete genome"></a>For a complete genome
The following 3 attributes are only valid for a comeplete (non-WGS) genome. Values corresponding to each sequence must be separated with commas (,) or semicolons (;), and the number of the values must match the number of sequences.
* `seq_names`  
  The sequence names to be used in a INSDC flat file. e.g. `"Chromosome, pXXX, pYYY"`
  The sequence name must not include any blanks. The plasmid name should start with lower case 'p' unless it is not known. In that case use 'unnamed', or ‘unnamed1’ & ‘unnamed2’ for distinct unnamed plasmids. 
  If you set `use_original_name = True`, the leading non-space letters in the query FASTA file will be taken as a sequence name, and the value specified by `seq_names` will be ignored.

* `seq_topologies`  
  Linear or circular. If not specified, the sequence will be treated as a linear sequence. Use `circular` or `c` to denote a circular sequence. For example, `"c,c,"` means the first and second replicons are circular, and the third one is linear (=`"circular,circular,linear"`).
* `seq_types`  
  Chromosome or plasmid assignment. If not specified, the sequence will be treated as a chromosome. Use `plasmid` or `p` to assign as a plasmid. For example, `",p,p"` means the first replicon is a chromosome, and the second and third replicons are plasmids (=`"chromosome,plasmid,plasmid"`).
  If a sequence is assigned as a plasmid, a "plasmid" qualifier with its plasmid name will appear in the source feature.

#### Advanced option
* `additional_modifiers` (optional)  
If you need additional source modifiers, you can specify them following the format like `"isolation_source=plant; culture_collection=JCM:XXXX; type_material=type strain of XXXX; note=ABC; note=DEF"`. The values specified by this attribute will be assigned to all the entries in the genome.

  __[Acceptable modifiers]__  
bio\_material,  collected\_by, collection\_date, country, culture\_collection,
ecotype, identified\_by, isolation\_source, host, note,
serotype, serovar, sub\_species, sub\_strain, type\_material, variety

  See also [Source Modifier List](https://www.ncbi.nlm.nih.gov/Sequin/modifiers.html) at the NCBI web site.




## Locus tag
Locus tags are identifiers that are systematically assigned to every gene in a genome.
A locus tag is represented by the combination of a prefix (locus\_tag\_prefix) and a tag value, 
separated by an under score "\_", eg A1C\_00001.  
Locus\_tag\_prefix can contain only alpha-numeric characters of at least 3 characters long,
and it should be registered at the INSDC in prior to genome submission.  
Usually, they are in a sequential order on the genome, 
and you can also leave gaps so that you can later add a new gene between gaps.
For more information, please read the [guideline](http://www.ddbj.nig.ac.jp/sub/locus_tag-e.html). 

You can adjust how locus tags should be assigned using the following parameters in `LOCUS_TAG_SETTINGS`.
* `locus_tag_prefix`: locus\_tag\_prefix
* `step`: The gap value between locus tags.
* `use_separate_tags`: If set to `True`, locus_tags are assigned separately by feature type. 
* `symbols`: Default setting = {"CDS": "", "rRNA": "r", "tRNA": "t", "tmRNA": "tm"},
  This has two functions. 
  First, genes whose feature type is included in the setting are assigned with locus_tag. Locus tags are not assigned to features not in the setting, such as assembly_gap and repeat_region. 
  Second, if `use_separate_tags` is set `True`, 
  the symbol specified by the setting will be added to tag values. 
  e.g. ABC_r0001 for rRNA, ABC_t0020 for tRNA.


## DDBJ submission
DFAST generates an annotation file (tab-separated) and a sequence file (FASTA-formatted),
which will be required to submit the genome through [DDBJ Mass Submission System](http://www.ddbj.nig.ac.jp/sub/mss_flow-e.html) (MSS). 


#### How to submit your genome
1. Before submission, you need to register BioProject and BioSample via the DDBJ submission portal [D-way](https://trace.ddbj.nig.ac.jp/D-way/). If necessary, raw sequence data should be deposited in Sequence Read Archive (SRA). A unique Locus Tag Prefix is required for submitting an annotated genome, which can be registered during the BioProject or BioSample registration process.

2. Run DFAST to create submission files. To create valid files, make sure source modifiers and locus tags are specified properly. Then, fill in the COMMON entry part of the annotation file with metadata such as contact information, references, and accession numbers of related databases (BioProject, BioSample, SRA).

3. DDBJ offers [data validation tools](https://www.ddbj.nig.ac.jp/sub/mss/massSub-e.html#tool) to submitters. Use the tools to validate your submission files before sending them to DDBJ.

#### Metadata
By default, the COMMON entry section of an annotation file is left blank, so that users can fill with the information later. Alternatively, if you provide a text file describing metadata to `metadata_file` in the config, DFAST can generate a fully qualified DDBJ submission file including metadata.
Alternatively, you can do it with the `--metadata_file` option.
* `metadata_file`
 A metadata file is a tab-separated text file with a key and a value in each line. It is recommended to use a sample metadata file (`$DFAST_APP_ROOT/example/sample.metadata.txt`) as a template.  
 See also `$DFAST_APP_ROOT/example/description.metadata.txt` for more information.

#### Example    

The below is an example of a command to generate compliant submission files for a complete genome consisting of one chromosome and 2 plasmids.  
Execute this in `$DFAST_APP_ROOT`.
```
dfast -g example/sample.lactobacillus.fna --complete t --organism "Lactobacillus hokkaidonensis" \
--strain LOOC260 --seq_names "Chromosome,pXXXX,pYYYY" --seq_topologies c,c,l --seq_types c,p,p \
--additional_modifiers "culture_collection=JCM:18461; isolation_source=silage; note=You can add a comment; note=You can add another comment; collection_date=2017-06-26" --metadata_file example/sample.metadata.txt \
--locus_tag_prefix LH260 --step 10 --use_separate_tags t --out LHLOOC
```


## GenBank submission
1. Register the genome project with the BioProject and the BioSample databases.
2. Create a submission template file [here](https://submit.ncbi.nlm.nih.gov/genbank/template/submission/).
3. Run DFAST to create two input files for the tbl2asn program, a feature table file (.tbl) and a sequence file (.fsa).
Set `center_name` (or `--center_name`) to specify your genome center tags, or you can change them later.
4.  Run tbl2asn to generate a .sqn file for submission to GenBank.
```
tbl2asn -a s -t template.sbt -i xxxxx.fsa -w xxxxx.cmt -V vb -Z discrepancy.txt
# xxxxx.cmt is an optional structured-comment file
```
5. Upload the files from [NCBI Submission Portal](https://submit.ncbi.nlm.nih.gov/subs/genome/).


For more detail, please see the [Submission Guideine](https://www.ncbi.nlm.nih.gov/genbank/genomesubmit/), [instruction manual](https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/) for tbl2asn and the [Genbank Submissions Handbook](https://www.ncbi.nlm.nih.gov/books/NBK51157/).

#### Example

The below is an example of a command to generate GenBank submission files for a draft genome of E. coli.
```
dfast -g your_genome.fa --organism "Escherichia coli" --strain Sakai \
--additional_modifiers "serovar=O157:H7; sub_strain=RIMD 0509952; country=Japan:Osaka, Sakai" \
--locus_tag_prefix ECS --center_name NIG 
```