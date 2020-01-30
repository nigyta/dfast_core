#! /usr/bin/env python
# coding:utf8

import multiprocessing


class Config:

    # This is a configuration file for the default workflow.
    # Some of the values specified in this file can be overridden by command line options.
    # You can customize the pipeine by defining your own configuration file, which can be called by -c or --config option.

    GENOME_FASTA = None  # Normally this value is left None, can be specified by -g or --genome option
    WORK_DIR = "OUT"  # Overridden by --out or -o option 

    CPU = min(max(int(multiprocessing.cpu_count() / 2), 1), 4)  # Number of CPUs to use. Default: Half of the number of processor cores.
    # Note that this is the global setting of CPU. You can also set different numbers of CPUs to each functional annotation processes. 

    FORCE_OVERWRITE = False  # Overridden to True by specifying --force option. 
    DEBUG = False  # OVerridden to True by specifying --debug option.

    GENOME_CONFIG = {
        # Set 'complete' to True to annotate a complete genome.
        # If 'use_original_name' is set to True, the leading non-space letters in the Fasta header line will be used as a sequence name.
        # Otherwise, sequences will be renamed as 'Sequence001, Sequence002, ...'.
        # If 'sort_by_length' is set to True, sequences will be sorted so that longer sequences come first.
        # In a draft genome, sequences shorter than 'minimum_length' will be eliminated.

        "complete": False,
        "use_original_name": False, # If set to True, the first word in the Fasta header line will be used as a sequence name.
        "sort_by_length": True,
        "minimum_length": 200
    }

    GENOME_SOURCE_INFORMATION = {
        # These attributes are reflected in the source feature,
        # and do not affect the annotation results.

        "organism": "",
        "strain": "",
        
        "seq_names": "",  # FOR COMPLETE GENOME. Use colons or semicolons to separate values. e.g. Chromosome,pXXX,pYYY

        "seq_types": "",  # c or chromosome for chromosomal sequences (default), p or plasmid for plasmidal sequences.
        # FOR COMPLETE GENOME. Use colons or semicolons to separate values. e.g. c,p,p

        "seq_topologies": "",  # l or linear for linear sequences (default), c or circular for circular sequences
        # FOR COMPLETE GENOME. Use colons or semicolons to separate values. e.g. c,c,l

        "additional_modifiers": "",
    }

    LOCUS_TAG_SETTINGS = {
        "locus_tag_prefix": "LOCUS",
        "step": 10,
        "use_separate_tags": True,  # If set to `True`, locus_tags are assigned separately according to feature type. 
        "symbols": {"CDS": "", "rRNA": "r", "tRNA": "t", "tmRNA": "tm"}
        # "symbols": {"CDS": "", "rRNA": "r", "tRNA": "t", "tmRNA": "tm", "nc_rna": "nc", "misc_rna": "misc"}
    }

    FEATURE_ADJUSTMENT = {
        "remove_partial_features": True,  # True: enabled, False: disabled

        "remove_overlapping_features": True,  # True: enabled, False: disabled
        "feature_type_priority": ["assembly_gap", "CRISPR", ("tmRNA", "tRNA", "rRNA"), "CDS"],

        "merge_cds": False,  # True: enabled, False: disabled
        "tool_type_priority": {"MGA": 0, "Prodigal": 1},
    }

    OUTPUT_RESULT = {
        # output verbosity level for .gff and .gbk
        "verbosity": 3  # 1: minimum, 2: standard, 3: rich
    }

    DDBJ_SUBMISSION = {
        "enabled": True,
        "output_verbosity": 1,
        "metadata_file": None
    }

    GENBANK_SUBMISSION = {
        "enabled": True,
        "center_name": "",  # Genome Center tag for GenBank submission
        "output_verbosity": 1,
    }


    STRUCTURAL_ANNOTATION = [
        {
            "tool_name": "GFF_import",
            "enabled": False,
            "options": {
                "targets": ["CDS"],
                "gff_file_name": None,
            },
        },
        {
            # GAP is a Gap Annotation Process that identifies gap regions (N's or n's runs) in the sequence.
            "tool_name": "GAP",
            "enabled": True,
            "target": "assembly_gap",
            "options": {
                "len_cutoff": 5,  # Gaps shorter than len_cutoff are ignored.
                "linkage_evidence": "paired-ends",  # You can change this as you like.
                "gap_type": "within scaffold"  # You can change this as you like.
            },
        },
        {
            # MetaGeneAnnotator for CDS prediction
            "tool_name": "MGA",
            "target": "CDS",
            "enabled": True,
            "options": {"cmd_options": "-s"},  # -s for single species, -m for multiple species
        },
        {
            # Aragorn for tRNA and tmRNA prediction
            "tool_name": "Aragorn",
            "enabled": True,
            "target": "tRNA",
            "options": {
                "gcode": "-gcbact",  # -gcbact for Bacterial/Plant Chloroplast genetic code, "-gcstd" for standard genetic code.
                "cmd_options": "-l",  # DFAST assumes linear sequences as input.
            },
        },
        {
            # tRNAscan-SE for tRNA prediction. By default, this is disabled.
            # Please insall tRNAscan-SE and put it in your PATH to enable this.
            "tool_name": "tRNAscan",
            "enabled": False,
            "target": "tRNA",
            "options": {
                "model": "--bact",  # --bact, --arch, --organ, --general
                "cmd_options": ""
            },
        },
        {
            # Barrnap for rRNA prediction
            "tool_name": "Barrnap",
             "enabled": True,
            "target": "rRNA",
            "options": {
                 # Currently, Barrnap will run with default settings.
                 # You can set parameters such as --reject and --lencutoff to cmd_options.
                 # "cmd_options": "--reject 0.4 --lencutoff 0.6"
             },
        },
        {
            # RNAmmer for rRNA prediction. By default, this is disabled.
            # Please insall RNAmmer and put it in your PATH to enable this.
            "tool_name": "RNAmmer",
            "enabled": False,
            "target": "rRNA",
            "options": {
                "model": "bac",  # arc/bac/euk
                "cmd_options": ""
            },
        },
        {
            # CRT for CRISPR detection
            # CRT runs on JavaVM. Depending on your system, it requires additional parameters to run.
            # Example
            #     "java_options": "-Xmx256m"
            "tool_name": "CRT",
            "enabled": True,
            "target": "CRISPR",
            "options": {
                "jar_file": "@@APP_ROOT@@/bin/common/CRT1.2-CLI.jar",
                "java_options": "",  # options passed to java. e.g. -Xmx256m
                "cmd_options": "",  # cmd_options passed to CRT.
            },
        },
        {
            # Prodigal for CDS prediction
            # By default Prodigal is disabled. To enable this, also set MGA disabled or enable merge_cds in FEATURE_ADJUSTMENT.
            "tool_name": "Prodigal",
            "enabled": False,
            "target": "CDS",
            "options": {
                "transl_table": 11,
                "cmd_options": "",
            },
         },
    ]

    FUNCTIONAL_ANNOTATION = [
        # Fucntional annotation steps will be conducted in the order specified in this list.
        # You can switch enabled/disabled, change the order, or add new steps.

        {
            # OrthoSearch (All-vs-all pairwise alignment between each reference genome to assign orthologous genes)
            # Normally, this should be run before other annotation steps.
            # In the default workflow, it is disabled. You can enable this by using the "--references" option. 
            "component_name": "OrthoSearch",
            "enabled": False,
            "options": {
                # "cpu": 2,  # Uncomment this to set the component-specific number of CPUs, or the global setting is used.
                "skipAnnotatedFeatures": False,
                "evalue_cutoff": 1e-6,
                "qcov_cutoff": 75,
                "scov_cutoff": 75,
                "aligner": "ghostx",  # ghostx, ghostz, or blastp
                "aligner_options": {},  # Normally, leave this empty. (Current version does not use this option.)
                "references": [
                ]
            },
        },
        {
            # By default, this is disabled. 
            # If you want to add your original databases to be searched prior to default DB,
            # set 'enabled' to True and specify 'database', which you can do with '--database' option.
            # The database file must be in a DFAST reference format,
            # and DB indexing must be done using the 'reference_util.py' script.
            "component_name": "DBsearch",
            "enabled": False,
            "options": {
                # "cpu": 2,  # Uncomment this to set the component-specific number of CPUs, or the global setting is used.
                "skipAnnotatedFeatures": True,
                "evalue_cutoff": 1e-6,
                "qcov_cutoff": 75,
                "scov_cutoff": 75,
                "aligner": "ghostx",  # ghostx, ghostz or blastp
                "aligner_options": {},  # Normally, leave this empty. (Current version does not use this option.)
                "database": "",
            },
        },
        {
            # Sequence homology search against the default DB.
            "component_name": "DBsearch",
            "enabled": True,
            "options": {
                # "cpu": 2,  # Uncomment this to set the component-specific number of CPUs, or the global setting is used.
                "skipAnnotatedFeatures": True,
                "evalue_cutoff": 1e-6,
                "qcov_cutoff": 75,
                "scov_cutoff": 75,
                "aligner": "ghostx",  # ghostz, ghostx or blastp
                "aligner_options": {},  # Normally, leave this empty. (Current version does not use this option.)
                "database": "@@APP_ROOT@@/db/protein/DFAST-default.ref",
            },
        },
        {
            # This is for homology search against a large-sized reference database,
            # such as pre-formatted Blast databases like RefSeq Protein.
            "component_name": "BlastSearch",
            "enabled": False,
            "options": {
                # "cpu": 2,  # Uncomment this to set the component-specific number of CPUs, or the global setting is used.
                "skipAnnotatedFeatures": False,
                "evalue_cutoff": 1e-6,
                "qcov_cutoff": 75,
                "scov_cutoff": 75,
                "aligner": "blastp",  # must be blastp
                "aligner_options": {},  # Normally, leave this empty. (Current version does not use this option.)
                "dbtype": "auto",  # Must be either of auto/ncbi/uniprot/plain
                "database": "",
            },
        },
        {
            # Pseudo/frame-shifted genes detection based on the upstream homology search results.
            # This also detects translation exceptions to selenocystein/pyrrolysine.
            "component_name": "PseudoGeneDetection",
            "enabled": True,
            "options": {
                "cpu": 1,  # Current version uses only 1 cpu, not supporting multi threading/processing.
                "skipAnnotatedFeatures": False,
                "extension": 300,
                "scov_cutoff": 85,
            },
        },
        {
            # Search against a profile-HMM database using hmmscan.
            # In the standard workflow, TIGRFAM HMM library will be searched.
            "component_name": "HMMscan",
            "enabled": True,
            "options": {
                # "cpu": 2,  # Uncomment this to set the component-specific number of CPUs, or the global setting is used.
                "skipAnnotatedFeatures": True,
                "evalue_cutoff": 1e-6,
                "database": "@@APP_ROOT@@/db/hmm/TIGRFAMs_15.0_HMM.LIB",
                "db_name": "TIGR",
            },
        },
        {
            # You can add another HMMscan, if you want.
            "component_name": "HMMscan",
            "enabled": False,
            "options": {
                # "cpu": 2,  # Uncomment this to set the component-specific number of CPUs, or the global setting is used.
                "skipAnnotatedFeatures": False,
                "evalue_cutoff": 1e-6,
                "db_name": "",  # eg 'Pfam',
                "database": ""  # eg '@@APP_ROOT@@/db/hmm/Pfam-A.hmm'
            },
        },
        {
            # Search against Conserved Domain Database (CDD) using RPS-BLAST
            # In the standard workflow, COG functional categories will be assigned by this step.
            "component_name": "CDDsearch",
            "enabled": True,
            "options": {
                # "cpu": 2,  # Uncomment this to set the component-specific number of CPUs, or the global setting is used.
                "skipAnnotatedFeatures": False,
                "evalue_cutoff": 1e-6,
                "database": "@@APP_ROOT@@/db/cdd/Cog",
                "rpsbproc_data": "@@APP_ROOT@@/bin/common/rpsbproc_data",  # Do not change this.
            },
        },
        {
            # You can add another CDDsearch, if you want.
            "component_name": "CDDsearch",
            "enabled": False,
            "options": {
                # "cpu": 2,  # Uncomment this to set the component-specific number of CPUs, or the global setting is used.
                "skipAnnotatedFeatures": False,
                "evalue_cutoff": 1e-6,
                "database": "",  # eg @@APP_ROOT@@/db/cdd/Prk
                "rpsbproc_data": "@@APP_ROOT@@/bin/common/rpsbproc_data",  # Do not change this.
            },
        },

    ]
