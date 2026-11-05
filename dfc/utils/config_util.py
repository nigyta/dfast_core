#! /usr/bin/env python
# coding: UTF8

import sys
import os
import shutil
import importlib.util
import pprint
from logging import getLogger

logger = getLogger(__name__)


def load_config(app_root, config_file, db_root=None):
    logger.info("Running on Python {}.".format(sys.version))
    logger.info("Loading a config file from {}".format(config_file))

    config = open(config_file).read()
    config = config.replace("@@APP_ROOT@@", app_root)
    db_root_env = os.getenv('DFAST_DB_ROOT')
    if db_root:
        logger.info("DB_ROOT is specified by a command-line option --dbroot %s", db_root)
    elif db_root_env:
        logger.info("DFAST_DB_ROOT is specified [%s]", db_root_env)
        db_root = db_root_env
    else:
        db_root = os.path.join(app_root, "db")
    config = config.replace("@@DB_ROOT@@", db_root)
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


def set_hmm_db(config, hmm_db):
    """Choose the NCBI HMM collection to search: "tigr" (TIGR models only, default) or "ncbifam" (all models)."""
    target = {"tigr": "TIGR", "ncbifam": "NCBIfam"}[hmm_db]
    settings = [setting for setting in config.FUNCTIONAL_ANNOTATION if setting.get("component_name", "") == "HMMscan"
                and setting.get("options", {}).get("db_name") in ("TIGR", "NCBIfam")]
    if not any(setting["options"]["db_name"] == target for setting in settings):
        logger.error("HMMscan with db_name '{}' is not defined in the config file. Aborting...".format(target))
        exit(1)
    for setting in settings:
        setting["enabled"] = setting["options"]["db_name"] == target
    logger.info("HMMscan will search the {} HMM database.".format(target))

def disable_hmm_scan(config):
    for setting in config.FUNCTIONAL_ANNOTATION:
        if setting.get("component_name", "") == "HMMscan":
            setting["enabled"] = False

def disable_cds_prediction(config):
    for setting in config.STRUCTURAL_ANNOTATION:
        if _is_cds_tool(setting):
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

def set_rrna_model(config, model):
    """Choose the rRNA profile HMMs for Barrnap: "barrnap" (Barrnap 0.9 models, default) or "rfam"."""
    for setting in config.STRUCTURAL_ANNOTATION:
        if setting.get("tool_name", "") == "Barrnap":
            setting.setdefault("options", {})["model"] = model
    logger.info("Barrnap will use the {} rRNA models.".format(model))


def set_rrna_kingdoms(config, kingdoms):
    """Choose the prokaryotic rRNA models for Barrnap: ["bac"] (default), ["arc"] or ["bac", "arc"]."""
    for setting in config.STRUCTURAL_ANNOTATION:
        if setting.get("tool_name", "") == "Barrnap":
            setting.setdefault("options", {})["kingdoms"] = list(kingdoms)
    logger.info("Barrnap will search rRNA models of {}.".format(", ".join(kingdoms)))


CDS_TOOLS = ("MGA", "Prodigal", "GeneMarkS2")


def _is_cds_tool(setting):
    # some configs (e.g. example/test_config.py) have no "target" keys, so also match by tool name
    return setting.get("target", "") == "CDS" or setting.get("tool_name", "") in CDS_TOOLS


def _select_CDS_prediction_tool(config, tool_name):
    for setting in config.STRUCTURAL_ANNOTATION:
        if _is_cds_tool(setting):
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

def disable_func_annotation(config):
    for setting in config.FUNCTIONAL_ANNOTATION:
        setting["enabled"] = False

def enable_metagenome_mode(config):
    logger.info("Metagenome mode is enabled.")
    for setting in config.STRUCTURAL_ANNOTATION:
        if _is_cds_tool(setting):
            tool_name = setting.get("tool_name", "")
            if tool_name == "MGA":
                if setting.get("enabled"):
                    logger.info("MGA will be invoked with '-m' option (for multiple species).")
                setting["options"]["cmd_options"] = "-m"
            if tool_name == "Prodigal":
                if setting.get("enabled"):
                    logger.info("Prodigal will be invoked with '-p meta' option.")
                setting["options"]["cmd_options"] = "-p meta"

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
    if value not in [4, 11, 25]:
        logger.warning("Unsupported genetic code {} is specified.".format(value))
    if value != 11:
        logger.info("Genetic code is set to {}.".format(value))
        for setting in config.STRUCTURAL_ANNOTATION:
            # CDS tools are included even without a "transl_table" key (MGA in example configs), so check_genetic_code sees it
            if "transl_table" in setting["options"] or _is_cds_tool(setting):
                setting["options"]["transl_table"] = value
        for setting in config.FUNCTIONAL_ANNOTATION:
            if setting.get("component_name", "") in ["PseudoGeneDetection", "NuclMapping"]:
                setting["options"]["transl_table"] = value


def check_genetic_code(config):
    """Exit if the enabled CDS predictor cannot handle the genetic code: MGA and metagenome mode support table 11 only."""
    for setting in config.STRUCTURAL_ANNOTATION:
        if not _is_cds_tool(setting) or not setting.get("enabled"):
            continue
        options = setting.get("options", {})
        tool_name = setting.get("tool_name", "")
        transl_table = options.get("transl_table", 11)
        cmd_options = options.get("cmd_options", "")
        if tool_name == "MGA" and transl_table != 11:
            logger.error("MGA supports only the genetic code 11, but {} is specified. Please use the default CDS predictor (Prodigal).".format(transl_table))
            exit(1)
        metagenome = (tool_name == "Prodigal" and "-p meta" in cmd_options) or (tool_name == "MGA" and "-m" in cmd_options.split())
        if metagenome and transl_table != 11:
            logger.error("Metagenome mode supports only the genetic code 11, but {} is specified.".format(transl_table))
            exit(1)


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

    keys = ["organism", "strain", "isolate", "seq_names", "seq_types", "seq_topologies", "additional_modifiers"]
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
    if project_type in ("gnm", "mag"):
        complete = True
    config.GENOME_CONFIG["complete"] = complete
    config.GENOME_CONFIG["project_type"] = project_type


MODULES = ("amr", "vf", "plasmid", "mge", "integron", "defense")

# module -> (component_name or tool_name of the config entry, db_name of the entry or None)
_MODULE_ENTRIES = {"amr": ("NuclMapping", "CARD"), "vf": ("NuclMapping", "VFDB"), "plasmid": ("Plasmidfinder", None),
                   "mge": ("MobileElementFinder", None), "integron": ("IntegronFinder", None),
                   "defense": ("DefenseFinder", None)}


def _find_entry(config, name, db_name):
    for setting in config.FUNCTIONAL_ANNOTATION + config.CONTIG_ANNOTATION:
        if name in (setting.get("component_name"), setting.get("tool_name")) \
                and (db_name is None or setting.get("options", {}).get("db_name") == db_name):
            yield setting


def _module_entry(config, module):
    """Return the config entry of a module, or None."""
    return next(_find_entry(config, *_MODULE_ENTRIES[module]), None)


def _nuclsearch_entry(config, module):
    """Return the enabled NuclSearch entry of the database of amr/vf, or None."""
    return next((s for s in _find_entry(config, "NuclSearch", _MODULE_ENTRIES[module][1]) if s.get("enabled")), None)


def _plasmidfinder_available():
    # Same alternatives as dfc/tools/dfast_plasmidfinder.py; its version (3.x) is checked when the tool starts.
    return bool(shutil.which("plasmidfinder.py") or importlib.util.find_spec("plasmidfinder"))


def module_unavailable(config, module):
    """Return why a module cannot run (what is missing and how to install it), or None if it can."""
    name, db_name = _MODULE_ENTRIES[module]
    if module in ("amr", "vf") and _nuclsearch_entry(config, module):
        # ponytail: a CARD/VFDB NuclSearch entry enabled in the config is used instead of the mapping; its reference
        # files are not checked here, NuclSearch reports them when it starts.
        return None
    entry = _module_entry(config, module)
    if entry is None:
        return "The config has no {0} entry.".format(name + (" " + db_name if db_name else ""))
    options = entry.get("options", {})
    if module in ("amr", "vf"):
        database = options.get("database", "")
        missing = [f for f in (database + ".nucl.ref", database + ".nucl.fasta") if not os.path.exists(f)]
        if missing:
            return ("{0} reference files are not found ({1}). Download them with 'dfast_file_downloader.py --modules {2}' "
                    "or 'scripts/reference_util_for_nucl.py --{3}'.".format(db_name, ", ".join(missing), module, db_name.lower()))
        return None
    problems = []
    data = None  # data directory the module needs, if any
    if module == "plasmid":
        if not _plasmidfinder_available():
            problems.append("PlasmidFinder 3.x is not found. Install it with 'pip install plasmidfinder'.")
        data = options.get("db_path", "")
    elif module == "mge":
        if not shutil.which("mefinder"):
            problems.append("'mefinder' is not found. Install it with 'pip install MobileElementFinder'.")
        data = options.get("db_path", "")
    elif module == "integron":
        if not shutil.which("integron_finder"):
            problems.append("'integron_finder' is not found. Install it with 'conda install -c bioconda integron_finder'.")
    else:
        if not shutil.which("defense-finder"):
            problems.append("'defense-finder' is not found. Install it with 'conda install -c bioconda defense-finder'.")
        data = options.get("models_dir", "")
    if data is not None and not os.path.exists(data):
        problems.append("The data of {0} are not found in '{1}'. Run 'dfast_file_downloader.py --modules {2}' to download them."
                        .format(name, data, module))
    return " ".join(problems) or None


def enable_modules(config, names):
    """Enable the accessory annotation modules named in `names` (--modules); "all" stands for every module.
    Returns the enabled modules in MODULES order. A module that cannot run is skipped with a warning when it is only
    selected by "all", and stops DFAST with an error when it is named explicitly."""
    explicit = set(names) - {"all"}
    selected = [m for m in MODULES if "all" in names or m in explicit]
    enabled, failed = [], False
    for module in selected:
        reason = module_unavailable(config, module)
        if reason is None:
            enabled.append(module)
        elif module in explicit:
            logger.error("--modules {0} cannot be used. {1}".format(module, reason))
            failed = True
        else:
            logger.warning("Module '{0}' is skipped. {1}".format(module, reason))
    if failed:
        exit(1)
    for module in enabled:
        if module in ("amr", "vf") and _nuclsearch_entry(config, module):
            logger.info("{} is searched with NuclSearch as set in the config, not mapped with NuclMapping.".format(_MODULE_ENTRIES[module][1]))
        else:
            _module_entry(config, module)["enabled"] = True
    logger.info("Modules enabled: {}".format(", ".join(enabled) or "none"))
    return enabled


def enable_mag(config, complete=False):
    project_type = "mag" if complete else "mag-wgs"
    config.GENOME_CONFIG["project_type"] = project_type
    if complete:
        config.GENOME_CONFIG["complete"] = True
    logger.info("MAG mode enabled (project_type={})".format(project_type))
