#! /usr/bin/env python
# coding: UTF8

import sys
import os
import pprint
from logging import getLogger

logger = getLogger(__name__)


def load_config(app_root, config_file):
    logger.info("Running on Python {}.".format(sys.version))
    logger.info("Loading a config file from {}".format(config_file))

    config = open(config_file).read()
    config = config.replace("@@APP_ROOT@@", app_root)
    exec(config, globals())  # Config object will be imported
    return Config


# def load_config(app_root, config_file):
#     logger.info("Running on Python {}.".format(sys.version))
#     logger.info("Loading a config file from {}".format(config_file))

#     version = sys.version_info.major
#     if version == 3:
#         from importlib import machinery
#         import tempfile
#         config = open(config_file).read()
#         config = config.replace("@@APP_ROOT@@", app_root)

#         with tempfile.NamedTemporaryFile(delete=True) as tf:
#             tf.write(config.encode('utf-8'))
#             temp_file_name = tf.name
#             module_ = machinery.SourceFileLoader('Config', temp_file_name).load_module()

#         return module_.Config

#     elif version == 2:
#         """
#         Python 2 does not have importlib.machinery
#         So far, a naive implementation for loading module.
#         """
#         config = open(config_file).read()
#         config = config.replace("@@APP_ROOT@@", app_root)
#         exec(config, globals())  # Config object will be imported
#         return Config
#     else:

#         raise AssertionError



def show_config(config):
    pp = pprint.PrettyPrinter(indent=1)

    attr_list = ['GENOME_FASTA', 'WORK_DIR', 'CPU', 'FORCE_OVERWRITE', 
             'GENOME_CONFIG', 'GENOME_SOURCE_INFORMATION', 'LOCUS_TAG_SETTINGS', 
             'OUTPUT_RESULT', 'DDBJ_SUBMISSION', 'GENBANK_SUBMISSION', 
             'STRUCTURAL_ANNOTATION', 'FEATURE_ADJUSTMENT', 'FUNCTIONAL_ANNOTATION']
    for attr_name in attr_list:
        attr = getattr(config, attr_name)
        if attr_name in ['STRUCTURAL_ANNOTATION', 'FUNCTIONAL_ANNOTATION']:
            sys.stdout.write("\n")
            sys.stdout.write(attr_name + ":" + "\n")
            settings = [setting for setting in attr if setting.get("enabled")]
            pp.pprint(settings)
        elif type(attr) == dict or type(attr) == list:
            sys.stdout.write("\n")
            sys.stdout.write(attr_name + ":" + "\n")
            pp.pprint(attr)
        else:
            sys.stdout.write(attr_name + ": " + str(attr) + "\n")


def set_references(config, references):
    for setting in config.FUNCTIONAL_ANNOTATION:
        if setting.get("component_name", "") == "OrthoSearch":
            references = references.strip(";, ").replace(",", ";").split(";")
            references = [x.strip() for x in references]
            setting["options"]["references"] = references
            setting["enabled"] = True
            return


def set_database(config, database):
    for setting in config.FUNCTIONAL_ANNOTATION:
        if setting.get("component_name", "") == "DBsearch":
            args = database.split(",")
            db_path = args[0]
            if len(args) > 1:
                db_name = args[1]
                setting["options"]["db_name"] = db_name
            else:
                setting["options"]["db_name"] = os.path.basename(db_path)
            if len(args) > 2:
                _set_threshold(args[2:], setting["options"])
            setting["options"]["database"] = db_path
            setting["enabled"] = True
            return


def set_threshold(config, args):
    args = args.split(",")
    for setting in config.FUNCTIONAL_ANNOTATION:
        if setting.get("component_name", "") == "DBsearch" and setting.get("default"):
            _set_threshold(args, setting["options"])


def _set_threshold(args, options):
    # helper function for set_threashold and set_database
    if len(args) > 0:
        arg = args[0]
        if arg:
            pident = float(arg)
            options["pident_cutoff"] = pident
    if len(args) > 1:
        arg = args[1]
        if arg:
            q_cov = float(arg)
            options["qcov_cutoff"] = q_cov
    if len(args) > 2:
        arg = args[2]
        if arg:
            s_cov = float(arg)
            options["scov_cutoff"] = s_cov
    if len(args) > 3:
        arg = args[3]
        if arg:
            e_value = float(arg)
            options["evalue_cutoff"] = e_value


def set_aligner(config, aligner):
    for setting in config.FUNCTIONAL_ANNOTATION:
        if "aligner" in setting.setdefault("options", {}):
            setting["options"]["aligner"] = aligner


def disable_cdd_search(config):
    for setting in config.FUNCTIONAL_ANNOTATION:
        if setting.get("component_name", "") == "CDDsearch":
            setting["enabled"] = False


def disable_hmm_scan(config):
    for setting in config.FUNCTIONAL_ANNOTATION:
        if setting.get("component_name", "") == "HMMscan":
            setting["enabled"] = False

def disable_cds_prediction(config):
    for setting in config.STRUCTURAL_ANNOTATION:
        if setting.get("target", "") == "CDS":
            setting["enabled"] = False

def disable_trna_prediction(config):
    for setting in config.STRUCTURAL_ANNOTATION:
        if setting.get("target", "") == "tRNA":
            setting["enabled"] = False

def disable_rrna_prediction(config):
    for setting in config.STRUCTURAL_ANNOTATION:
        if setting.get("target", "") == "rRNA":
            setting["enabled"] = False

def disable_crispr_prediction(config):
    for setting in config.STRUCTURAL_ANNOTATION:
        if setting.get("target", "") == "CRISPR":
            setting["enabled"] = False

def enable_trnascan(config, model):
    for setting in config.STRUCTURAL_ANNOTATION:
        if setting.get("tool_name", "") == "tRNAscan":
            setting["enabled"] = True
            model = "-A" if model == "arch" else "-B"
            setting["options"]["model"] = model
        if setting.get("tool_name", "") == "Aragorn":
            setting["enabled"] = False

def enable_rnammer(config, model):
    # model should be arc/bac
    if model == "bact":
        model = "bac"
    elif model == "arch":
        model = "arc"
    for setting in config.STRUCTURAL_ANNOTATION:
        if setting.get("tool_name", "") == "RNAmmer":
            setting["enabled"] = True
            setting["options"]["model"] = model
        if setting.get("tool_name", "") == "Barrnap":
            setting["enabled"] = False


def _select_CDS_prediction_tool(config, tool_name):
    for setting in config.STRUCTURAL_ANNOTATION:
        if setting.get("target", "") == "CDS":
            if setting.get("tool_name", "") == tool_name:
                setting["enabled"] = True
            else:
                setting["enabled"] = False

def enable_prodigal(config):
    _select_CDS_prediction_tool(config, "Prodigal")

def enable_mga(config):
    _select_CDS_prediction_tool(config, "MGA")

def enable_genemarks2(config, genome_type):
    # genome type must be bacteria or archaea
    if genome_type == "bact":
        genome_type = "bacteria"
    elif genome_type == "arch":
        genome_type = "archaea"
    _select_CDS_prediction_tool(config, "GeneMarkS2")
    gms2_config = [conf for conf in config.STRUCTURAL_ANNOTATION if conf["tool_name"] == "GeneMarkS2"][0]
    gms2_config["options"]["genome_type"] = genome_type

def set_gff(config, gff_file_name):
    targets = []
    for setting in config.STRUCTURAL_ANNOTATION:
        if setting.get("tool_name", "") == "GFF_import":
            setting["enabled"] = True
            setting["options"]["gff_file_name"] = gff_file_name
            targets = setting["options"]["targets"]
    for setting in config.STRUCTURAL_ANNOTATION:
        if setting.get("tool_name", "") != "GFF_import":
            target = setting.get("target", "-")
            if target in targets:
                setting["enabled"] = False

def set_genetic_code(config, value):
    if value not in [4, 11]:
        logger.warning("[WARNING] Unsupported genetic code {} is specified.".format(value))
    if value == 11:
        pass
    else:
        # genetic code 4 can be specified when using Prodigal/GeneMarkS2/Aragorn
        logger.warning("Genetic code is set to {}.".format(value))
        for setting in config.STRUCTURAL_ANNOTATION:
            if "transl_table" in setting["options"]:
                setting["options"]["transl_table"] = value
        for setting in config.FUNCTIONAL_ANNOTATION:
            if setting.get("component_name", "") == "PseudoGeneDetection":
                setting["options"]["transl_table"] = value




def set_values_from_metadata(config):
    metadata_file = config.DDBJ_SUBMISSION.get("metadata_file")
    if not metadata_file or not os.path.exists(metadata_file):
        return
    logger.info("Loading metadata from {}. (Will be overidden by command-line options)".format(metadata_file))
    D = {}
    with open(metadata_file) as f:
        for line in f:
            key, value = line.strip("\n").split("\t")
            if value:
                D[key] = value

    keys = ["organism", "strain", "seq_names", "seq_types", "seq_topologies", "additional_modifiers"]
    for key in keys:
        value = D.get(key)
        if value:
            config.GENOME_SOURCE_INFORMATION[key] = value

    keys = ["locus_tag_prefix"]
    for key in keys:
        value = D.get(key)
        if value:
            config.LOCUS_TAG_SETTINGS[key] = value

    complete = D.get("complete", "false")
    complete = True if complete.lower() in ["t", "true"] else False
    project_type = D.get("project_type") or D.get("projectType", "")
    if project_type == "gnm":
        complete = True    
    config.GENOME_CONFIG["complete"] = complete