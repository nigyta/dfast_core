#! /usr/bin/env python
# coding: UTF8


from collections import OrderedDict
import re
import os.path


class MetadataField(object):

    def __init__(self, name, description, qualifier, feature, entry, type_="string", mss_required=False, pattern=re.compile(r".*"), value=""):
        """
            name: a name for a variable used in DFAST (case sensitive, must be unique)
            description: description of this field
                ex) name: locusTagPrefix  -> description: Locus Tag Prefix
            qualifier is a term used in MSS format. ex) ab_name
            feature: a feature label where this metadata are described. ex) SUBMITTER, REFERENCE, ST_COMMENT
            entry: an entry label where this metadata are described.  COMMON|SEQUENCE|OTHER
            array: Is set True if multiple values are acceptable. Must be splitted with semicolons. ex) "spam; ham; eggs"
            mss_required: Is used to check if a generated file is in the valid DDBJ-MSS format.
            pattern: a regex pattern that must MATCH the value
            value: a value for this field
        """
        if type_ == "array":
            value = [x.strip() for x in value.split(";") if x.strip()]
        self.name, self.description, self.qualifier, self.feature, self.entry, self.type, \
            self.mss_required, self.pattern, self.value = \
            name, description, qualifier, feature, entry, type_, mss_required, pattern, value

    def __repr__(self):
        if self.type == "array":
            values = "; ".join(self.value)
            return "{self.name}\t{values}".format(self=self, values=values)
        else:
            return "{self.name}\t{self.value}".format(self=self)

    def append_value(self, value):
        assert self.type == "array"
        self.value.append(value)

    def set_value(self, value):
        if self.type == "array":
            self.value = [x.strip() for x in value.split(";") if x.strip()]
        else:
            self.value = value

    def set_default(self, value, default_value):
        value = value or default_value
        self.set_value(value)

    def getAdditionalField(self, newName, newValue=""):
        if ":" not in newName:
            print("no colon : in the parameter name")
            raise AssertionError
        num = newName.split(":")[-1]
        newFeature = "{0}:{1}".format(self.feature, newName.split(":")[-1])
        newObj = MetadataField(newName, self.description, self.qualifier, newFeature, self.entry,
                               self.type, self.mss_required, self.pattern, newValue)
        return newObj

    def validate(self, ignoreNull=True):
        # if ignoreNull and (not self.mss_required) and (self.value == "" or self.value == []):
        if self.value == "" or self.value == []:
            if ignoreNull or (not self.mss_required):
                return {}
            else:
                # return {self.name: {"msg": "Empty value", "qualifier": self.qualifier, "feature": self.feature}}
                return {"MISSING_VALUES": {"key": self.name}}

        failedValues = {}
        if self.type == "array":
            for eachVal in self.value:
                if self.pattern.match(eachVal) is None:
                    failedValues.setdefault("values", []).append(eachVal)
                    failedValues["pattern"] = self.pattern.pattern.replace("^(", "").replace(")$", "")

        else:
            if self.pattern.match(self.value) is None:
                failedValues.setdefault("values", []).append(self.value)
                failedValues["pattern"] = self.pattern.pattern.replace("^(", "").replace(")$", "")

        if len(failedValues) > 0:
            failedValues.update({"key": self.name})
            return {"INCORRECT_VALUES": failedValues}
        else:
            return {}

    def render(self):
        RET = []

        if self.value == "" or self.value == []:
            RET.append(["", "", "", self.qualifier, ""])
            if self.qualifier == "ab_name":
                # if ab_name has no value, render additional 2 empty lines into the template
                RET.append(["", "", "", self.qualifier, ""])
                RET.append(["", "", "", self.qualifier, ""])
        elif self.type == "array":
            for eachVal in self.value:
                if eachVal:
                    RET.append(["", "", "", self.qualifier, eachVal])
        elif self.type == "boolean" and self.value:
            if self.value != "NO":
                RET.append(["", "", "", self.qualifier, ""])
        else:
            RET.append(["", "", "", self.qualifier, self.value])
        return RET


class Metadata(object):
    METADATA_DEFINITION_FILE = os.path.dirname(os.path.abspath(__file__)) + "/metadata_definition.tsv"

    def __init__(self, metadata_dict):
        self.fields = OrderedDict()
        self.refNumber = 0
        self.commentNumber = 0

        for line in open(self.__class__.METADATA_DEFINITION_FILE):
            if line.startswith("#"):
                continue
            else:
                name, description, qualifier, feature, entry, type_, mss_required, pattern, value = line.strip("\n").split("\t")[0:9]
                mss_required = True if mss_required == "TRUE" else False
                pattern = re.compile(r"^({0})$".format(pattern))
                field = MetadataField(name, description, qualifier, feature, entry, type_, mss_required, pattern, value)
                self.fields[name] = field

        for key, value in metadata_dict.items():
            field = self.fields.get(key)
            if field:
                if field.feature not in ["COMMENT"]:
                    field.set_value(value)

        self.addComments(metadata_dict)
        # self.addReferences(metadata_dict)

    @staticmethod
    def load(file_name):
        D = {}
        for line in open(file_name):
            key, value = line.strip("\n").split("\t")
            if value:
                D[key] = value
        return Metadata(D)

    def get_value(self, field_name, default_value="", idx=0):
        field = self.fields.get(field_name)
        if field is None:
            return default_value
        if field.type == "array":
            if len(field.value) > idx:
                return field.value[idx] or default_value
            else:
                return default_value
        else:
            return field.value or default_value

    def set_value(self, field_name, value):
        field = self.fields.get(field_name)
        if field:
            field.set_value(value)
        else:
            raise IndexError("Field name not found")


    def getFields(self):
        return list(self.fields.keys())  # self.fields is a dictionary of metadata

    def toTSV(self, outputFileName):
        Buffer = ""
        for field in self.fields.values():
            Buffer += repr(field) + "\n"
        with open(outputFileName, "w") as f:
            f.write(Buffer)

    def validateValues(self, ignoreNull=True):

        errors = {"status": "success",
                  "MISSING_VALUES": [],
                  "INCORRECT_VALUES": [],
                  "INCONSISTENT_VALUES": [],
                  }
        for field in self.fields.values():
            ret = field.validate(ignoreNull)
            # ret should be like: {"INCORRECT_VALUE": {"key": hoge, "value": [ham, spam, egg]}
            for key, val in ret.items():
                errors[key].append(val)

        if len(errors["MISSING_VALUES"]) > 0 or len(errors["INCORRECT_VALUES"]) > 0 or len(errors["INCONSISTENT_VALUES"]) > 0:
            errors["status"] = "fail"
        return errors

    def render_common_entry(self, dfast_version=None, complete=False):

        def render_category():
            if complete:
                pass
            else:
                outputData.append(["", "DATATYPE", "", "type", "WGS"])
                outputData.append(["", "KEYWORD", "", "keyword", "WGS"])
                outputData.append(["", "KEYWORD", "", "keyword", self.get_value("keyword", "STANDARD_DRAFT")])

        def renderFeature(outputData, featureType):
            # assert featureType in ["DBLINK", "SUBMITTER", "REFERENCE"]
            RET = []
            for field in self.fields.values():
                if field.feature == featureType:
                    if field.value or field.mss_required:
                        RET += field.render()
            if len(RET) > 0:
                RET[0][1] = featureType.split(":")[0]
                outputData += RET

        def addDfastComment(outputData):
            if not dfast_version is None:
                outputData.append(["", "COMMENT", "", "line", "Annotated by DFAST v.{} (https://dfast.nig.ac.jp)".format(dfast_version)])

        outputData = []

        render_category()
        renderFeature(outputData, "DBLINK")
        renderFeature(outputData, "SUBMITTER")
        renderFeature(outputData, "REFERENCE")   # render REFERENCE regardless of refNumber
        for i in range(1, self.refNumber):   # render additional REFERENCES
            renderFeature(outputData, "REFERENCE:" + str(i + 1))
        for i in range(self.commentNumber):
            if i == 0:
                renderFeature(outputData, "COMMENT")
            else:
                renderFeature(outputData, "COMMENT:" + str(i + 1))
        addDfastComment(outputData)
        renderFeature(outputData, "ST_COMMENT")
        renderFeature(outputData, "DATE")

        outputData[0][0] = "COMMON"
        return outputData


    def addComments(self, metadataDict):
        checkNext = True
        while checkNext:
            if self.commentNumber == 0:
                fieldName = "comment"
                value = metadataDict.get("comment") or metadataDict.get("comment:1")

                if value:
                    self.fields[fieldName].set_value(value)
                    self.commentNumber += 1
                else:
                    self.fields[fieldName].value = []
                    self.commentNumber += 1

            else:
                fieldName = "comment:{}".format(str(self.commentNumber + 1))
                value = metadataDict.get(fieldName)
                originalField = self.fields["comment"]
                if fieldName in metadataDict:
                    commentField = originalField.getAdditionalField(fieldName, value)
                    self.fields[fieldName] = commentField
                    self.commentNumber += 1
                else:
                    checkNext = False



if __name__ == '__main__':
    # metadataFlle = "/Users/ytanizaw/project/labrep_dev/jobs/53eeb0b5-edc5-427d-99c2-b213d6c187a0/result/metadata.txt"


    D = {"keyword": "HIGH_QUALITY_DRAFT", "comment": "OK", "comment:1": "COM1", "comment:2": "COM2; COM2-2", "comment:3": "COMMENT3 is here.; COM3-2; COM3-3",
         "reference": "REF1", "reference:2": "REF2"}
    metadata = Metadata(D)
    out = metadata.render_common_entry()
    for row in out:
        print("\t".join(row))

