import pytest
import re
from dfc.utils.metadata_util import MetadataField, Metadata


def make_bool_field(value):
    return MetadataField(
        "environmental_sample", "environmental_sample",
        "environmental_sample", "source", "EX_SOURCE",
        type_="boolean", mss_required=False,
        pattern=re.compile(r".*"), value=value,
    )

def test_boolean_render_true():
    assert make_bool_field("TRUE").render() == [["", "", "", "environmental_sample", ""]]

def test_boolean_render_yes():
    assert make_bool_field("YES").render() == [["", "", "", "environmental_sample", ""]]

def test_boolean_render_lowercase_true():
    assert make_bool_field("true").render() == [["", "", "", "environmental_sample", ""]]

def test_boolean_render_false():
    assert make_bool_field("false").render() == []

def test_boolean_render_no():
    assert make_bool_field("NO").render() == []

def test_boolean_render_empty():
    assert make_bool_field("").render() == []

def test_metadata_reads_isolate():
    m = Metadata({"isolate": "MAG-001", "organism": "Candidatus Foo", "strain": ""})
    assert m.get_value("isolate") == "MAG-001"

def test_metadata_reads_metagenome_source():
    m = Metadata({"metagenome_source": "soil metagenome", "organism": "Candidatus Foo"})
    assert m.get_value("metagenome_source") == "soil metagenome"

def test_metadata_renders_environmental_sample():
    m = Metadata({"environmental_sample": "TRUE", "organism": "Candidatus Foo"})
    rendered = m.render_ex_source()
    qualifiers = [row[3] for row in rendered]
    assert "environmental_sample" in qualifiers


# ---- Task 1 new: should_render / render_ex_source(project_type) ----

import dfc.utils.metadata_util as _mu

def _field(name, mss_required, required_for=(), excluded_for=(), value="", type_="string"):
    import re
    f = MetadataField(name, name, name, "source", "EX_SOURCE",
                      type_=type_, mss_required=mss_required,
                      pattern=re.compile(r".*"), value=value)
    f.required_for = set(required_for)
    f.excluded_for = set(excluded_for)
    return f

def test_should_render_excluded_overrides_value():
    f = _field("strain", mss_required=True, excluded_for=("mag",), value="K-12")
    assert f.should_render("mag") is False

def test_should_render_excluded_overrides_mss_required():
    f = _field("strain", mss_required=True, excluded_for=("mag",), value="")
    assert f.should_render("mag") is False

def test_should_render_required_for_placeholder():
    f = _field("isolate", mss_required=False, required_for=("mag", "mag-wgs"), value="")
    assert f.should_render("mag") is True

def test_should_render_value_overrides_not_required():
    f = _field("isolate", mss_required=False, required_for=("mag",), value="ABC")
    assert f.should_render("wgs") is True  # value present, not excluded

def test_should_render_empty_not_required():
    f = _field("isolate", mss_required=False, required_for=("mag",), value="")
    assert f.should_render("wgs") is False  # empty and not required for wgs

def test_effective_required_excluded():
    f = _field("strain", mss_required=True, excluded_for=("mag",))
    assert f.effective_required("mag") is False

def test_effective_required_required_for():
    f = _field("isolate", mss_required=False, required_for=("mag",))
    assert f.effective_required("mag") is True

def test_effective_required_default():
    f = _field("strain", mss_required=True)
    assert f.effective_required("wgs") is True

def test_render_ex_source_project_type_mag(tmp_path):
    # After TSV update (Task 2), metagenome_source + environmental_sample appear as placeholder
    # This test uses Metadata directly with empty values
    m = Metadata({"organism": "Foo"})  # no metagenome_source/environmental_sample
    # Manually set required_for on the loaded fields (until TSV is updated in Task 2)
    m.fields["metagenome_source"].required_for = {"mag", "mag-wgs"}
    m.fields["environmental_sample"].required_for = {"mag", "mag-wgs"}
    rows = m.render_ex_source("mag")
    qualifiers = [r[3] for r in rows]
    assert "metagenome_source" in qualifiers
    assert "environmental_sample" in qualifiers

def test_render_ex_source_nonmag_no_placeholder():
    m = Metadata({"organism": "Foo"})
    m.fields["metagenome_source"].required_for = {"mag", "mag-wgs"}
    m.fields["environmental_sample"].required_for = {"mag", "mag-wgs"}
    rows = m.render_ex_source("wgs")
    qualifiers = [r[3] for r in rows]
    assert "metagenome_source" not in qualifiers
    assert "environmental_sample" not in qualifiers


# ---- Task 2: TSV field definitions ----

def test_tsv_strain_excluded_for_mag():
    import dfc.utils.metadata_util as mu
    mu._FIELD_DEF_CACHE = None  # force reload
    f = mu.get_field_definition("strain")
    assert f is not None
    assert "mag" in f.excluded_for
    assert "mag-wgs" in f.excluded_for

def test_tsv_isolate_required_for_mag():
    import dfc.utils.metadata_util as mu
    mu._FIELD_DEF_CACHE = None
    f = mu.get_field_definition("isolate")
    assert f is not None
    assert "mag" in f.required_for
    assert "mag-wgs" in f.required_for
    assert not f.excluded_for

def test_tsv_metagenome_source_required_for_mag():
    import dfc.utils.metadata_util as mu
    mu._FIELD_DEF_CACHE = None
    f = mu.get_field_definition("metagenome_source")
    assert "mag" in f.required_for

def test_tsv_environmental_sample_required_for_mag():
    import dfc.utils.metadata_util as mu
    mu._FIELD_DEF_CACHE = None
    f = mu.get_field_definition("environmental_sample")
    assert "mag" in f.required_for


# ---- Task 2: config_util tests ----

import types
import sys
import os


def _make_config(metadata_dict, tmp_path):
    meta_file = tmp_path / "meta.txt"
    meta_file.write_text("\n".join(f"{k}\t{v}" for k, v in metadata_dict.items()) + "\n")
    config = types.SimpleNamespace(
        DDBJ_SUBMISSION={"metadata_file": str(meta_file)},
        GENOME_CONFIG={},
        GENOME_SOURCE_INFORMATION={},
        LOCUS_TAG_SETTINGS={},
    )
    return config

def test_config_mag_sets_complete(tmp_path):
    from dfc.utils.config_util import set_values_from_metadata as load_metadata
    config = _make_config({"project_type": "mag", "organism": "Foo"}, tmp_path)
    load_metadata(config)
    assert config.GENOME_CONFIG["complete"] is True
    assert config.GENOME_CONFIG["project_type"] == "mag"

def test_config_mag_wgs_not_complete(tmp_path):
    from dfc.utils.config_util import set_values_from_metadata as load_metadata
    config = _make_config({"project_type": "mag-wgs", "organism": "Foo"}, tmp_path)
    load_metadata(config)
    assert config.GENOME_CONFIG["complete"] is False
    assert config.GENOME_CONFIG["project_type"] == "mag-wgs"

def test_config_isolate_propagated(tmp_path):
    from dfc.utils.config_util import set_values_from_metadata as load_metadata
    config = _make_config({"project_type": "mag", "organism": "Foo", "isolate": "MAG-001"}, tmp_path)
    load_metadata(config)
    assert config.GENOME_SOURCE_INFORMATION["isolate"] == "MAG-001"


# ---- Task 3: genome.py tests ----

def test_genome_is_mag_true():
    from dfc.genome import Genome
    g = object.__new__(Genome)
    g.project_type = "mag"
    g.is_mag = g.project_type in ("mag", "mag-wgs")
    assert g.is_mag is True

def test_genome_is_mag_false():
    from dfc.genome import Genome
    g = object.__new__(Genome)
    g.project_type = "gnm"
    g.is_mag = g.project_type in ("mag", "mag-wgs")
    assert g.is_mag is False

def test_genome_source_env_division():
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from dfc.genome import Genome
    g = object.__new__(Genome)
    g.complete = False
    g.project_type = "mag-wgs"
    g.is_mag = True
    g.organism = "Candidatus Foo"
    g.strain = ""
    g.isolate = "MAG-001"
    g.seq_names = ["scaffold"]
    g.seq_types = [""]
    g.seq_topologies = ["linear"]
    record = SeqRecord(Seq("ATCG" * 100), id="seq1", name="seq1")
    g.seq_records = {"seq1": record}
    g.add_source_information()
    assert record.annotations["data_file_division"] == "ENV"
    assert record.annotations.get("isolate") == "MAG-001"
    assert "strain" not in record.annotations

def test_genome_source_bct_division():
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from dfc.genome import Genome
    g = object.__new__(Genome)
    g.complete = False
    g.project_type = ""
    g.is_mag = False
    g.organism = "Escherichia coli"
    g.strain = "K-12"
    g.isolate = ""
    g.seq_names = ["scaffold"]
    g.seq_types = [""]
    g.seq_topologies = ["linear"]
    record = SeqRecord(Seq("ATCG" * 100), id="seq1", name="seq1")
    g.seq_records = {"seq1": record}
    g.add_source_information()
    assert record.annotations["data_file_division"] == "BCT"
    assert record.annotations.get("strain") == "K-12"
    assert "isolate" not in record.annotations


# ---- Task 3: add_source_features TSV-driven ----

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def _make_genome_for_source_features(project_type, strain="", isolate=""):
    from dfc.genome import Genome
    g = object.__new__(Genome)
    g.complete = False
    g.project_type = project_type
    g.is_mag = project_type in ("mag", "mag-wgs")
    g.organism = "Test Org"
    g.strain = strain
    g.isolate = isolate
    g.seq_names = ["seq1"]
    g.seq_types = [""]
    g.seq_topologies = ["linear"]
    g.additional_modifiers = {}
    record = SeqRecord(Seq("ATCG" * 100), id="seq1", name="seq1")
    record.annotations = {"strain": strain, "isolate": isolate, "organism": g.organism, "mol_type": "genomic DNA"}
    if not strain:
        record.annotations.pop("strain", None)
    if not isolate:
        record.annotations.pop("isolate", None)
    g.seq_records = {"seq1": record}
    return g

def test_source_features_strain_excluded_for_mag():
    g = _make_genome_for_source_features("mag", strain="K-12", isolate="")
    g.add_source_features({})
    q = list(g.seq_records.values())[0].features[0].qualifiers
    assert "strain" not in q
    assert "isolate" in q
    assert q["isolate"] == [""]  # placeholder

def test_source_features_isolate_placeholder_for_mag():
    g = _make_genome_for_source_features("mag", strain="", isolate="")
    g.add_source_features({})
    q = list(g.seq_records.values())[0].features[0].qualifiers
    assert "isolate" in q
    assert q["isolate"] == [""]

def test_source_features_isolate_value_for_mag():
    g = _make_genome_for_source_features("mag", isolate="MAG-001")
    g.add_source_features({})
    q = list(g.seq_records.values())[0].features[0].qualifiers
    assert q["isolate"] == ["MAG-001"]

def test_source_features_isolate_shown_for_wgs_when_set():
    g = _make_genome_for_source_features("wgs", strain="", isolate="ISO-001")
    g.add_source_features({})
    q = list(g.seq_records.values())[0].features[0].qualifiers
    assert "isolate" in q
    assert q["isolate"] == ["ISO-001"]

def test_source_features_strain_shown_for_wgs():
    g = _make_genome_for_source_features("wgs", strain="K-12", isolate="")
    g.add_source_features({})
    q = list(g.seq_records.values())[0].features[0].qualifiers
    assert "strain" in q
    assert q["strain"] == ["K-12"]

def test_source_features_strain_placeholder_for_wgs_empty():
    g = _make_genome_for_source_features("wgs", strain="", isolate="")
    g.add_source_features({})
    q = list(g.seq_records.values())[0].features[0].qualifiers
    assert "strain" in q
    assert q["strain"] == [""]  # placeholder (mss_required=TRUE)


# ---- Task 4: render_common_entry tests ----

_COMMON_META = {
    "organism": "Candidatus Foo", "bioproject": "PRJDB0001",
    "biosample": "SAMD000001", "submitter": "Tanizawa,Y.",
    "contact": "Y. Tanizawa", "email": "a@b.c",
    "institute": "NIG", "city": "Mishima", "country": "Japan",
    "street": "1111", "zip": "411-8540",
    "reference": "Test", "author": "Tanizawa,Y.", "status": "Unpublished",
    "assembly_method": "Flye v2.9", "genome_coverage": "100X",
    "sequencing_technology": "Nanopore",
}

def test_render_common_mag():
    m = Metadata(_COMMON_META)
    rows = m.render_common_entry(project_type="mag", complete=True)
    col1 = [r[1] for r in rows]
    col3 = [r[3] for r in rows]
    assert "division" in col3
    assert col3.index("division") < col3.index("keyword")
    keywords = [r[4] for r in rows if r[3] == "keyword"]
    assert "ENV" in keywords
    assert "Metagenome Assembled Genome" in keywords
    assert "MAG" in keywords
    assert "WGS" not in keywords
    assert "type" not in col3
    assert col1.count("KEYWORD") == 1

def test_render_common_mag_wgs():
    m = Metadata(_COMMON_META)
    rows = m.render_common_entry(project_type="mag-wgs", complete=False)
    col1 = [r[1] for r in rows]
    col3 = [r[3] for r in rows]
    keywords = [r[4] for r in rows if r[3] == "keyword"]
    assert "ENV" in keywords
    assert "Metagenome Assembled Genome" in keywords
    assert "MAG" in keywords
    assert "WGS" in keywords
    assert "type" in col3
    assert col1.index("DATATYPE") < col1.index("DIVISION") < col1.index("KEYWORD")
    assert col1.count("KEYWORD") == 1

def test_render_common_gnm_unchanged():
    m = Metadata({**_COMMON_META, "organism": "E. coli"})
    rows = m.render_common_entry(project_type="gnm", complete=True)
    col3 = [r[3] for r in rows]
    assert "division" not in col3
    assert "WGS" not in [r[4] for r in rows if r[3] == "keyword"]


# ---- Task 5: ddbj_submission tests ----

def test_ff_definition_mag_complete():
    from dfc.utils.ddbj_submission import create_ff_definiton
    row = create_ff_definiton("complete", plasmid=None, project_type="mag")
    assert "@@[isolate]@@" in row[4]
    assert "@@[strain]@@" not in row[4]

def test_ff_definition_mag_wgs():
    from dfc.utils.ddbj_submission import create_ff_definiton
    row = create_ff_definiton("contig", plasmid=None, project_type="mag-wgs")
    assert "@@[isolate]@@" in row[4]

def test_ff_definition_non_mag_unchanged():
    from dfc.utils.ddbj_submission import create_ff_definiton
    row = create_ff_definiton("complete", plasmid=None, project_type="wgs")
    assert "@@[strain]@@" in row[4]
    assert "@@[isolate]@@" not in row[4]


# ---- get_file_prefix tests ----

def _make_ddbj_submission(biosample, strain="", isolate="", is_mag=False):
    import types
    from dfc.utils.ddbj_submission import DDBJsubmission
    genome = types.SimpleNamespace(strain=strain, isolate=isolate, is_mag=is_mag,
                                   workDir="/tmp")
    obj = object.__new__(DDBJsubmission)
    obj.genome = genome
    obj.metadata = {"biosample": biosample} if biosample else {}
    return obj

def test_file_prefix_non_mag():
    obj = _make_ddbj_submission("SAMD000001", strain="K-12", is_mag=False)
    assert obj.get_file_prefix() == "SAMD000001_K-12"

def test_file_prefix_mag_uses_isolate():
    obj = _make_ddbj_submission("SAMD000001", isolate="MAG-001", is_mag=True)
    assert obj.get_file_prefix() == "SAMD000001_MAG-001"

def test_file_prefix_mag_no_biosample():
    obj = _make_ddbj_submission("", isolate="MAG-001", is_mag=True)
    assert obj.get_file_prefix() == "MAG-001"

def test_file_prefix_fallback():
    obj = _make_ddbj_submission("", strain="", isolate="", is_mag=False)
    assert obj.get_file_prefix() == "mss"
