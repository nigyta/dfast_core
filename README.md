# DFAST - DDBJ Fast Annotation and Submission Tool
DFAST is a flexible and customizable pipeline for prokaryotic genome annotation as well as data submission to the INSDC. It is originally developed as the background engine for the [DFAST web service](https://dfast.nig.ac.jp) and is also available as a stand-alone command-line tool.
The stand-alone version of DFAST is also refered to as DFAST-core to differentiate it from its on-line version.  
For inquiry and request, please contact us at `dfast @ nig.ac.jp`.

## Important Notice 2025 Feb
The reference data for DFAST is normally available from our web service (https://dfast.ddbj.nig.ac.jp). However, due to a system replacement on our institute’s supercomputer, the web service will be unavailable from mid-February to early March 2025. During this period, the `dfast_file_downloader.py` script will not work. Instead, please manually download the data (`dfast_core_db.tar.gz`) from https://dfast.annotation.jp and follow the instructions in the README file available at the site.

---

#### Contents
* [Overview](#overview)
* [Installation](#installation)
* [Installation via conda](#installation-via-conda)
* [How to run](#howto)
* [Default workflow](#workflow)
* [Options](#options)
* [Software distribution](#distribution)
* [Trouble shoot](#trouble_shoot)
* [How to run within a Docker container](#howtorundocker)
* [Citation](#citation)
* [FAQ](docs/FAQ.md)

#### Advanced contents
* [Workflow details](docs/workflow.md)
* [INSDC submission](docs/insdc_submission.md)
* [Cookbook](docs/cookbook.md)


## <a id="overview"></a>Overview

* **Easy install**  
DFAST is implemented in Python and runs on Mac and Linux. No additional modules are required other than BioPython. It comes with external binaries for the default workflow.  
Bioconda package is also available.

* **Flexible and customizable**  
You can customize the pipeline as you like by specifying parameters, gene prediction tools, and reference databases in the configuraition file.
Each of annotation processes is defined as a Python module with common interfaces,
which facilitates future development and incorporation of new tools.

* **Fast and rich annotation**  
DFAST can annotate a typical-sized bacterial genome within several minutes. In addition to the conventional homology search, it features unique functions such as orthologous gene assignment between reference genomes, pseudo/frameshifted gene prediction, and conserved domain search.
* **INSDC submission**  
As its name suggested, DFAST is intended to support rapid genome submission to the INSDC through DDBJ. DFAST generates submission files for DDBJ Mass Submission System (MSS).

## <a id="installation"></a>Installation
If you use Anaconda/Miniconda, see [here](#condainstallation) to install using `conda`.

### Prerequisites
* **Python (3.7-)**  
  DFAST runs both on Python 3.7 or later. Python 2 is no longer supported.  
  Python>=3.10 is recommended.
* **BioPython package**  
  You can install this with the Python package management tool `pip`:  
  ```
  (sudo) pip install biopython
  ```
  If `pip` is not available, please follow the [instruction](http://biopython.org/wiki/Download) of BioPython.  

* **Perl and Java**  
Some of the external programs called from DFAST depend on Perl or Java. Basically, they work with the pre-installed versions on your system.  
For **RedHat/CentOS/Fedora**, the Time::Piece module might be required:
  ```
  sudo yum install perl-Time-Piece
  ```

### Source code
Available from the GitHub repository [nigyta/dfast_core](https://github.com/nigyta/dfast_core).

* **Via git command** (recommended)  
  ```
  git clone https://github.com/nigyta/dfast_core.git
  cd dfast_core    # Hereafter, we call this directory $DFAST_APP_ROOT
  ```


* **Download the distribution**  
Download the latest DFAST distribution from [GitHub Tags](https://github.com/nigyta/dfast_core/tags), then unarchive it.
  ```
  wget https://github.com/nigyta/dfast_core/archive/refs/tags/x.x.2.tar.gz  
  tar xvfz x.x.x.tar.gz  # Files will be uncompressed into dfast_core-x.x.x direcotory   
  cd dfast_core-x.x.x    # Hereafter, we call this directory $DFAST_APP_ROOT
  ```

For your convenience, create links to DFAST executables in a directory specified by the `PATH` environment variable. For example,
```
ln -s $DFAST_APP_ROOT/dfast /usr/local/bin/
ln -s $DFAST_APP_ROOT/scripts/dfast_file_downloader.py /usr/local/bin/
```

### Reference databases
  After downloading the source code, prepare reference databases using the bundled utility script.  
  By default, database files will be generated into the directory under $DFAST_APP_ROOT/db/. You can also change the location of the directory by specifying either `--dbroot` option or `DFAST_DB_ROOT` environmental variable.
1. **Default protein database**
    ```
    dfast_file_downloader.py --protein dfast
    ```
    File downloading and database indexing for GHOSTX and BLASTP will be performed. 
2. **HMMer and RPS-BLAST databases (this may take time)**
    ```
    dfast_file_downloader.py --cdd Cog --hmm TIGR
    ```
    DFAST default workflow requires COG database for RPS-BLAST and TIGRFAM database for hmmerscan.
* **See help for more information.**
    ```
    dfast_file_downloader.py -h
    ```

The `dfast_file_downloader.py` script downloads the reference data from our web service (https://dfast.ddbj.nig.ac.jp). If file downloads fail due to server maintenance or other issues, please manually obtain the reference data from [this site](https://dfast.annotation.jp).

## Installation via conda
DFAST is also available from [Bioconda](https://bioconda.github.io/recipes/dfast/README.html). Install with:
```
conda install -c bioconda -c conda-forge dfast
```
We recommend specifying the latest version. See available versions from [here](https://github.com/nigyta/dfast_core/tags).
```
conda install -c bioconda -c conda-forge dfast=1.X.XX
```
If this does not work, please try to install DFAST into the fresh conda environment.

DFAST executables are added to the `PATH` environmental variable, and the software package is installed in the `opt` directory under the Anaconda/Miniconda root directory. (e.g. /home/USER/miniconda3/opt/dfast-X.X.X/)  

After installing DFAST, download the reference databases:
```
dfast_file_downloader.py --protein dfast --cdd Cog --hmm TIGR
```

## <a id="howto"></a>How to run
1. **Help**  
    ```
    dfast -h
    ```
    or by specifying the Python interpreter,
    ```
    python $DFAST_APP_ROOT/dfast -h
    ```

2. **Test run**  
    ```
    dfast --config $DFAST_APP_ROOT/example/test_config.py
    ```
    This minimum workflow includes CDS prediction and database search against the default protein database using the GHOSTX aligner. The result will be generated in `RESULT_TEST` dierctory.  
    If not working properly, please check if the default database is installed. Normally, it finishes within a minute.

3.  **Basic usage**  
    ```
    dfast --genome path/to/your_genome.fna(.gz)
    ```
    This invokes the DFAST pipeline with the default workflow defined in $DFAST_APP_ROOT/dfc/default_config.py. DFAST accepts a FASTA-formatted genome sequence file as a query.  

4. **Advanced usage**  
    By providing command line options, you can override the default settings described in the configuration file.
    ```
    dfast --genome your_genome.fna(.gz) --organism "Escherichia coli" --strain "str. xxx" \
    --locus_tag_prefix ECXXX --minimum_length 200 --references EC_ref_genome.gbk \
    --aligner blastp --out OUT_ECXXX
    ```
    'locus tag prefix' is required if you want your genome to be submitted to the INSDC (use `--locus_tag_prefix` option). DFAST generates DDBJ submission files. For more information, please refer to [INSDC submission](docs/insdc_submission.md).
     If you set `--references` option, OrthoSearch (orthologous gene assignment) is enabled, which conducts all-against-all protein alignments between given reference genomes to infer orthologous genes.  
     `--aligner blastp` will let DFAST use BLASTP for protein alignments instead of default GHOSTX.

     These optional values can be specified in a configuration file, saving you from providing them as command line options. See the following step. 

5. **More advanced usage: Creating your own workflow**  
An easy way to do this is to copy and edit the default configuration file, which is located in $DFAST_APP_ROOT/dfc/default_config.py.
The configuration file is a self-explanatory Python script, in which the workflow is defined using basic Python objects like lists and dictionaries.

    You can call your original configuration file with the `--config` option.
    ```
    dfast --genome your_genome.fna(.gz) --config your_config.py
    ```

## <a id="workflow"></a>Default workflow
DFAST default annotation workflow accepts a genomic FASTA file (draft or complete) as an input and includes following processes. Read [Workflow](docs/workflow.md) to learn more.

### Structural annotation
The following tools are run in parallel to predict biological features (e.g. CDSs and RNAs). After that, partial and overlapping features will be cleaned up.
* CDS prediction (MetaGeneAnnotator)
* rRNA prediction (Barrnap)
* tRNA/tmRNA prediction (Aragorn)
* CRISPR prediction (CRT)
* Assembly gaps within sequences

Optionally, you can choose Prodigal/GeneMarkS2, RNAmmer, tRNAscan-SE to predict CDS, rRNA, tRNA, respectively. See [FAQ](docs/FAQ.md). (You need to install them manually.)

### Functional annotation
1. OrthoSearch (Optional. Set `--references` option to enable this.)
2. DBsearch using the Ghostx aligner against the DFAST default database
3. PseudoGeneDetection (internal stop codons and frameshifts)
4. HMMscan against the profile HMM database of TIGRFAM
5. CDDsearch against COG database from NCBI Conserved Domain Database

By default, GHOSTX is used to align protein sequences. Diamond/BLASTP can be used optionally. See [FAQ](docs/FAQ.md). (Diamond needs to be installed manually.) 

### Output
* Sequence and annotation data in GFF3 and GenBank format
* Sequence data in FASTA format
* Statistics for genome sequences and annotated features
* DDBJ submission files

## <a id="options"></a>Options

```
Basic usage:
  usage: dfast -g your_genome.fna [options]
  
Basic options:
  -g PATH, --genome PATH
                        Genomic FASTA file for input. Can be gzipped.
  -o PATH, --out PATH   Output directory (default:OUT)
  -c PATH, --config PATH
                        Configuration file (default config will be used if not specified)
  --organism STR        Organism name
  --strain STR          Strain name

Genome settings:
  --complete BOOL       Treat the query as a complete genome. Not required unless you need INSDC submission files. [t|f(=default)]
  --use_original_name BOOL
                        Use original sequence names in a query FASTA file [t|f(=default)]
  --sort_sequence BOOL  Sort sequences by length [t(=default)|f]
  --minimum_length INT  Minimum sequence length (default:200)
  --fix_origin          Rotate/flip the chromosome so that the dnaA gene comes first. (ONLY FOR A FINISHED GENOME)
  --offset INT          Offset from the start codon of the dnaA gene. (for --fix_origin option, default=0)

Locus_tag settings:
  --locus_tag_prefix STR
                        Locus tag prefix (defaut:LOCUS)
  --step INT            Increment step of locus tag (default:10)
  --use_separate_tags BOOL
                        Use separate tags according to feature types [t(=default)|f]

Workflow options:
  --threshold STR       Thresholds for default database search (format: "pident,q_cov,s_cov,e_value", default: "0,75,75,1e-6")
  --database PATH       Additional reference database to be searched against prior to the default database. (format: db_path[,db_name[,pident,q_cov,s_cov,e_value]])
  --references PATH     Reference file(s) for OrthoSearch. Use semicolons for multiple files, e.g. 'genome1.faa;genome2.gbk'
  --aligner STR         Aligner to use [ghostx(=default)|blastp|diamond]
  --use_prodigal        Use Prodigal to predict CDS instead of MGA
  --use_genemarks2 STR  Use GeneMarkS2 to predict CDS instead of MGA. [auto|bact|arch]
  --use_trnascan STR    Use tRNAscan-SE to predict tRNA instead of Aragorn. [bact|arch]
  --use_rnammer STR     Use RNAmmer to predict rRNA instead of Barrnap. [bact|arch]
  --gcode INT           Genetic code [11(=default),4(=Mycoplasma)]
  --no_func_anno        Disable all functional annotation steps
  --no_hmm              Disable HMMscan
  --no_cdd              Disable CDDsearch
  --no_cds              Disable CDS prediction
  --no_rrna             Disable rRNA prediction
  --no_trna             Disable tRNA prediction
  --no_crispr           Disable CRISPR prediction
  --metagenome          Set options of MGA/Prodigal for metagenome contigs
  --amr                 [Preliminary implementation] Enable AMR/VFG annotation and identification of plasmid-derived contigs
  --gff GFF             [Preliminary implementation] Read GFF to import structural annotation. Ignores --use_original_name, --sort_sequence, --fix_origin.

Genome source modifiers and metadata [advanced]:
  These values are only used to create INSDC submission files and do not affect the annotation result. See documents for more detail.

  --seq_names STR       Sequence names for each sequence (for complete genome)
  --seq_types STR       Sequence types for each sequence (chromosome/plasmid, for complete genome)
  --seq_topologies STR  Sequence topologies for each sequence (linear/circular, for complete genome)
  --additional_modifiers STR
                        Additional modifiers for source features
  --metadata_file PATH  Path to a metadata file (optional for DDBJ submission file)

Run options:
  --cpu INT             Number of CPUs to use
  --use_locustag_as_gene_id
                        Use locustag as gene ID for FASTA and GFF. (Useful when providing DFAST results to other tools such as Roary)
  --dbroot PATH         DB root directory (default:APP_ROOT/db
  --force               Force overwriting output
  --debug               Run in debug mode (Extra logging and retaining temporary files)
  --show_config         Show pipeline configuration and exit
  --version             Show program version
  -h, --help            Show this help message
```

## <a id="distribution"></a>Software distribution
DFAST is freely available as open-source under the GPLv3 license (See [LICENSE](docs/LICENSE)).

This distribution contains following external programs.
* [MetaGeneAnnotator](http://metagene.cb.k.u-tokyo.ac.jp/) (© Hideki Noguchi)  
 Redistributed by courtesy of Hideki Noguchi at National Institute of Genetics.
* [Aragorn](http://mbio-serv2.mbioekol.lu.se/ARAGORN/) (GPLv3)
* [Barrnap](https://github.com/tseemann/barrnap) (GPLv3)
* [CRT](http://www.room220.com/crt/) (Public domain)
* [GHOSTX](http://www.bi.cs.titech.ac.jp/ghostx/) (BSD-2-Clause)
* [GHOSTZ](http://www.bi.cs.titech.ac.jp/ghostz/) (CC BY 4.0)
* blastp, makeblastdb, blastdbcmd, rpsblast, rpsbproc from [NCBI-BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) package. (Public domain)
* hmmpress, hmmscan from [HMMer](http://hmmer.org/) package (GPLv3)
* [LAST](http://last.cbrc.jp/) (GPLv3)

## <a id="trouble_shoot"></a>Trouble shoot
* DBsearch is slow  
The default aligner GHOSTX is fast but requires a large amount of memory. In our environment, it uses 1.8Gbyte memory per process.  
If your machine does not have enough memory, decrease the number of CPUs (`--cpu 2` or `--cpu 1`) or use BLASTP instead (`--aligner blastp`). 
* GLIBCXX not found error on Linux system  
If your system is old, DFAST will abort with the message "/usr/lib64/libstdc++.so.6: version 'GLIBCXX_3.4.15' not found".  
In this case, you need to update "libstdc++.so.6". (You might need to install a newer version of GCC.)  
Please check the file as following: `strings /usr/lib64/libstdc++.so.6 | grep GLIBCXX`
* libidn-11 on ArchLinux  
According to the report from users, DFAST fails on ArchLinux due to `libidn-11` required for BLASTP. You may need to install `libidn-133-compat` from the AUR repository.

## <a id="howtorundocker"></a>How to run DFAST within a Docker container.
The Docker container image is available from [Dockerhub:nigyta/dtast_core](https://hub.docker.com/r/nigyta/dfast_core/tags) and [quay.io:biocontainers/dfast](https://quay.io/repository/biocontainers/dfast?tab=tags).  
Use `--dbroot` to specity the location of the reference data.
Download the reference data:
```
docker run --rm -v PATH/TO/DB:/dfast_db nigyta/dfast_core:latest dfast_file_downloader.py --protein dfast --cdd Cog --hmm TIGR --dbroot /dfast_db
```

Invoke DFAST:
```
docker run --rm -v PATH/TO/DB:/dfast_db -v PATH/TO/YOUR/DATA:/data nigyta/dfast_core:latest dfast --genome /data/your_genome.fa --out /data/your_result --dbroot /dfast_db
```

## Experimental work
### Annotation for antibiotic registance genes and virulence fators  
`blastn` and `PlasmidFinder` are required. Please install it by yourself.  
__Usage__
1. Prepare [CARD](https://card.mcmaster.ca), [VFDB](http://www.mgc.ac.cn/VFs/), and [PlasmidFinder](https://bitbucket.org/genomicepidemiology/plasmidfinder_db/src/master/) reference data.
```
scripts/dfast_file_downloader.py --plasmidfinder
scripts/reference_util_for_nucl.py --card --card_version 3.2.9 --vfdb --vfdb_update_date 2024-05-10
```
Since VFDB provides only the latest version, the value specified with `--vfdb_update_date` is used only as a timestamp for the reference data.  
See [CARD/download](https://card.mcmaster.ca/download) for the latest version of CARD and [VFDB/download](http://www.mgc.ac.cn/VFs/download.htm) for the updated date of VFDB.

2. Run
Invoke DFAST with `--amr` to enable `NuclSearch` for CARD/VFDB and `ContigAnnotation` using `PlasmidFinder`
```
dfast -g example/pOXA-48.fa --amr
```

### DFAST Record JSON  
DFAST Record JSON is a JSON-formatted file that contains annotation results of DFAST, and will be used for DDBJ data submission in the future. To generate this, Python>=3.10 is required and [DR_tools](https://github.com/ddbj/dr_tools) must be installed. 
```
pip install "git+https://github.com/ddbj/dr_tools.git"
```


## <a id="citation"></a>Citation
* on-line version of DFAST  
    DFAST and DAGA: web-based integrated genome annotation tools and resources  
    *Biosci Microbiota Food Health*. 2016; 35(4): 173–184.  
    Yasuhiro TANIZAWA, Takatomo FUJISAWA, Eli KAMINUMA, Yasukazu NAKAMURA, and Masanori ARITA  
* stand-alone version (DFAST-core)  
    DFAST: a flexible prokaryotic genome annotation pipeline for faster genome publication.  
    *Bioinformatics*; 2018; 34(6): 1037–1039.  
    Yasuhiro TANIZAWA, Takatomo FUJISAWA, Yasukazu NAKAMURA  
    https://academic.oup.com/bioinformatics/article/34/6/1037/4587587


