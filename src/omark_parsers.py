from src.error_check import operation_failed


def omark_stats(omark):
    error = operation_failed(omark)
    results = {}
    if error:
        results["OMArk Consistency Results"] = error
        results["OMArk Completeness Results"] = error
        results["OMArk Species Composition"] = error
        return results
    clades = []
    with open(omark["OMARK"]["outfile"]) as fhand:
        for line in fhand:
            if "The clade used was" in line:
                clade = line.strip().split()[-1]+"-HOGs"
            if "Number of conserved HOGs" in line:
                hogs = line.strip().split()[-1]
            if "Single:" in line:
                single = line.strip().split()[-1]
            if "Duplicated:" in line:
                dup = line.rstrip().split()[-1]
            if "Duplicated, Unexpected:" in line:
                dup_un = line.strip().split()[-1]
            if "Duplicated, Expected:" in line:
                dup_exp = line.strip().split()[-1]
            if "Missing:" in line:
                missing = line.strip().split()[-1]
            if "Total Consistent" in line:
                consistent = line.strip().split()[-1]
            if "Consistent, partial hits" in line:
                cons_partial = line.strip().split()[-1]
            if "Consistent, fragmented" in line:
                cons_frag = line.strip().split()[-1]
            if "Total Inconsistent" in line:
                inconsistent = line.strip().split()[-1]
            if "Inconsistent, partial hits" in line:
                inconsistent_partial =  line.strip().split()[-1]
            if "Inconsistent, fragmented" in line:
                inconsistent_fragmented = line.strip().split()[-1]
            if "Total Contaminants" in line:
                contaminants = line.strip().split()[-1]
            if "Contaminants, partial hits" in line:
                contaminants_partial = line.strip().split()[-1]
            if "Contaminants, fragmented" in line: 
                contaminants_fragmented = line.strip().split()[-1]
            if "Total Unknown" in line:
                unkown = line.strip().split()[-1]
            if "Clade" in line:
                clade_comp = line.split(":")[-1].strip()
            if "associated query proteins" in line:
                clade_per = line.strip().split()[-1]
                clades.append("{}: {}".format(clade_comp, clade_per))
    consistency_results = "Cons:{}[P:{};F:{}],Inco:{}[P:{},F:{}],Cont:{},Unkn:{}".format(consistent,
                                                                                 cons_partial,
                                                                                 cons_frag,
                                                                                 inconsistent,
                                                                                 inconsistent_partial,
                                                                                 inconsistent_fragmented,
                                                                                 contaminants,
                                                                                 unkown).replace("(", "").replace(")", "")
    results["OMArk Consistency Results"] = consistency_results

    completness_results = "{}: {}; S:{},D:{}[U:{},E:{}],M:{}".format(clade, hogs, single,
                                                                     dup, dup_un, dup_exp,
                                                                     missing).replace("(", "").replace(")", "")
    results["OMArk Completeness Results"] = completness_results
    results["OMArk Species Composition"] = "; ".join(clades).replace("(", "").replace(")", "")
    return results
            
