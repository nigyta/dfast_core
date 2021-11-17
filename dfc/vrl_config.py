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
        "use_separate_tags": False,  # If set to `True`, locus_tags are assigned separately according to feature type. 
        "symbols": {"CDS": "" , "stem_loop": ""}
        # "symbols": {"CDS": "", "rRNA": "r", "tRNA": "t", "tmRNA": "tm", "nc_rna": "nc", "misc_rna": "misc"}
    }

    FEATURE_ADJUSTMENT = {
        "remove_partial_features": False,  # True: enabled, False: disabled

        "remove_overlapping_features": False,  # True: enabled, False: disabled
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
        "output_verbosity": 3,
        "metadata_file": None
    }

    GENBANK_SUBMISSION = {
        "enabled": False,
        "center_name": "",  # Genome Center tag for GenBank submission
        "output_verbosity": 1,
    }


    STRUCTURAL_ANNOTATION = [
        {
            "tool_name": "BlastFeatureN",
            "enabled": True,
            "options": {
                # "targets": ["CDS", "tRNA", "rRNA"],
                # "targets": ["CDS"],  # feature types to be imported
                "targets": ["CDS", "5'UTR", "3'UTR", "stem_loop", "mat_peptide"],  # feature types to be imported
                "gbk_file_name": "/Users/tanizawa/projects/newpipe/covid19/ref/genome.gb",
            },
        },
        {
            # GAP is a Gap Annotation Process that identifies gap regions (N's or n's runs) in the sequence.
            "tool_name": "GAP",
            "target": "assembly_gap",
            "enabled": True,
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
            "enabled": False,
            "options": {"cmd_options": "-s"},  # -s for single species, -m for multiple species
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
                "pident_cutoff": 0,
                "aligner": "ghostx",  # ghostx, ghostz or blastp
                "aligner_options": {},  # Normally, leave this empty. (Current version does not use this option.)
                "database": "",
                "db_name": "",
            },
        },
        {
            # Sequence homology search against the default DB.
            "component_name": "DBsearch",
            "enabled": False,
            "default": True,
            "options": {
                # "cpu": 2,  # Uncomment this to set the component-specific number of CPUs, or the global setting is used.
                "skipAnnotatedFeatures": True,
                "evalue_cutoff": 1e-6,
                "qcov_cutoff": 75,
                "scov_cutoff": 75,
                "pident_cutoff": 0,
                "aligner": "ghostx",  # ghostz, ghostx or blastp
                "aligner_options": {},  # Normally, leave this empty. (Current version does not use this option.)
                "database": "@@APP_ROOT@@/db/protein/DFAST-default.ref",
                "db_name": "",
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
                "pident_cutoff": 0,
                "aligner": "blastp",  # must be blastp
                "aligner_options": {},  # Normally, leave this empty. (Current version does not use this option.)
                "dbtype": "auto",  # Must be either of auto/ncbi/uniprot/plain
                "database": "",
                "db_name": "",
            },
        },
        {
            # Pseudo/frame-shifted genes detection based on the upstream homology search results.
            # This also detects translation exceptions to selenocystein/pyrrolysine.
            "component_name": "PseudoGeneDetection",
            "enabled": False,
            "options": {
                "cpu": 1,  # Current version uses only 1 cpu, not supporting multi threading/processing.
                "skipAnnotatedFeatures": False,
                "extension": 300,
                "scov_cutoff": 85,
                "transl_table": 11,
                # "genetic_code_file": "@@APP_ROOT@@/bin/common/transl_table_11.txt"
            },
        },

    ]
