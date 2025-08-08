import os
import ete3

from pathlib import Path

from importlib.resources import files


BUSCO_LINEAGES = Path(files('GAQET').joinpath("docs/busco_lineages.txt"))
BULLET_OK = "\t✓\t"
BULLET_FIX = "\tERROR!\t"
HEADER = "-"*5


def check_required_inputs(yaml):
    errors = []
    required_inputs = ["Assembly", "Annotation"]
    if "ID" not in yaml:
        errors.append(BULLET_FIX + "Field {} is required".format("ID"))
    elif not yaml["ID"]:
        errors.append(BULLET_FIX + "Field {} is empty".format("ID"))
    for input in required_inputs:
        if input not in yaml:
            errors.append(BULLET_FIX + "Field {} is required".format(input))
        else:
            path = Path(yaml[input])
            if not path.is_file():
                errors.append(BULLET_FIX + "Path for field {} doesn't exists: {}".format(input, str(path)))
    if len(errors) == 0:
        errors = [BULLET_OK + "All required inputs are present"]
    return errors

def check_required_reviewer_inputs(yaml):
    errors = []
    required_inputs = ["Assembly", "Annotation"]
    if "Species" not in yaml:
        errors.append(BULLET_FIX + "Field {} is required".format("ID"))
    elif not yaml["Species"]:
        errors.append(BULLET_FIX + "Field {} is empty".format("ID"))
    for input in required_inputs:
        if input not in yaml:
            errors.append(BULLET_FIX + "Field {} is required".format(input))
        else:
            path = Path(yaml[input])
            if not path.is_file():
                errors.append(BULLET_FIX + "Path for field {} doesn't exists: {}".format(input, str(path)))
    if "Basedir" not in yaml:
        errors.append(BULLET_FIX + "Field {} is required".format("Basedir"))
    else:
        path = Path(yaml[input])
        if not path.is_dir():
            errors.append(BULLET_FIX + "Path for field {} doesn't exists: {}".format(input, str(path)))
    if len(errors) == 0:
        errors = [BULLET_OK + "All required inputs are present"]
    return errors


def check_busco_lineages(yaml):
    errors = []
    available_lineages = [lineage.split()[0].split(".")[0] for lineage in open(BUSCO_LINEAGES)]
    for lineage in yaml["BUSCO_lineages"]:
        if lineage not in available_lineages:
            errors.append(BULLET_FIX + "BUSCO lineage {} doesn't exists")
    if len(errors) == 0:
        errors.append(BULLET_OK + "BUSCO lineages are valid")
    return errors


def check_available_analysis(yaml):
    errors = []
    available_analysis = ["AGAT", "BUSCO", "PSAURON",
                              "DETENGA", "OMARK", "PROTHOMOLOGY"]
    if not yaml["Analysis"]:
        return [BULLET_FIX + "No analysis found in YAML config file"]
    else:
        for analysis in yaml["Analysis"]:
            if analysis.strip() not in available_analysis:
                errors.append(BULLET_FIX + "{} is not a valid analysis".format(analysis))
    if len(errors) == 0:
        errors = [BULLET_OK + "All analysis are valid"]
    return errors


def check_taxid(yaml):
    not_defined = False
    if "OMARK_taxid" not in yaml:
        not_defined = True
    elif not yaml["OMARK_taxid"]:
        not_defined = True
    if not_defined:
        return [BULLET_FIX + "OMARK_taxid field is not defined"]
    else:
        taxid = yaml["OMARK_taxid"]
        try:
            ncbi = ete3.NCBITaxa()
            linid = ncbi.get_lineage(taxid)
        except ValueError:
            return [BULLET_FIX + "NCBI taxid {} is not valid".format(taxid)]
    return [BULLET_OK + "Taxid for OMARK is valid"]


def check_prothomology_dbs(yaml):
    errors = []
    not_defined = False
    if "PROTHOMOLOGY_tags" not in yaml:
        not_defined = True
    elif not yaml["PROTHOMOLOGY_tags"]:
        not_defined = True
    if not_defined:
        return [BULLET_FIX + "PROTHOMOLOGY_tags field is not defined"]
    else:
        for db in yaml["PROTHOMOLOGY_tags"]:
            for tag, path in db.items():
                if not Path(path).is_file():
                    errors.append(BULLET_FIX + "Protein database for tag {} doesn't exists".format(tag))
    if len(errors) == 0:
        errors = [BULLET_OK +  "All protein databases are valid"]
    return errors
    

def check_OMARK_db(yaml):
    errors = []
    not_defined = False
    if "OMARK_db" not in yaml:
        not_defined = True
    elif not yaml["OMARK_db"]:
        not_defined = True
    if not_defined:
        return [BULLET_FIX + "OMARK_db field is not defined"]
    else:
        path = Path(yaml["OMARK_db"])
        if not path.is_file():
            errors.append(BULLET_FIX + "OMARK_db database {} doesn't exists".format(str(path)))
    if len(errors) == 0:
        errors = [BULLET_OK +  "OMARK_db {} found".format(str(path))]
    return errors

def check_detenga_db(yaml):
    errors = []
    not_defined = False
    if "OMARK_db" not in yaml:
        not_defined = True
    elif not yaml["DETENGA_db"]:
        not_defined = True
    if not_defined:
        return [BULLET_FIX + "DETENGA_db field is not defined"]
    else:
        detenga_dbs = ["rexdb-plant", "rexdb-metazoa", "rexdb"]
        if yaml["DETENGA_db"] not in detenga_dbs:
            errors.append(BULLET_FIX + "DETENGA_db database {} doesn't exists. Available options are {}".format(yaml["DETENGA_db"], ",".join(detenga_dbs)))
    if len(errors) == 0:
        errors = [BULLET_OK +  "DETENGA_db {} found".format(yaml["DETENGA_db"])]
    return errors


def report_yaml_file(yaml):
    report = []
    report += [HEADER + "Checking if all required inputs are present" + HEADER]
    report += check_required_inputs(yaml)
    report += [HEADER + "Checking if all analysis are valid" + HEADER]
    report += check_available_analysis(yaml)
    if yaml["Analysis"]:
        if "BUSCO" in yaml["Analysis"]:
            report += [HEADER + "Checking if BUSCO lineages are valid" + HEADER]
            report += check_busco_lineages(yaml)
        if "OMARK" in yaml["Analysis"]:
            report += [HEADER + "Checking if OMARK taxid is valid" + HEADER]
            report += check_taxid(yaml)
            report += [HEADER + "Checking if OMARK db is available"]
            report += check_OMARK_db(yaml)
        if "DETENGA" in yaml["Analysis"]:
            report += [HEADER + "Checking if DeTEnGA db is available"]
            report += check_detenga_db(yaml)
        if "PROTHOMOLOGY" in yaml["Analysis"]:
            report += [HEADER + "Checking if protein databases exists" + HEADER]
            report += check_prothomology_dbs(yaml)
    return "\n".join(report)


def report_yaml_reviewer_file(yaml):
    report = []
    report += [HEADER + "Checking if all required inputs are present" + HEADER]
    report += check_required_reviewer_inputs(yaml)
    report += [HEADER + "Checking if all analysis are valid" + HEADER]
    return "\n".join(report)