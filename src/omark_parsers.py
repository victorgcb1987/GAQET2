def omark_stats(omark):
    results = {}
    clades = []
    print(omark["OMARK"]["outfile"])
    with open(omark["OMARK"]["outfile"]) as fhand:
        for line in fhand:
            if "The clade used was" in line:
                clade = line.strip().split()[-1]
                print(clade)
            if "Number of conserved HOGs" in line:
                hogs = line.strip().split()[-1]
                print(hogs)   
            if "Single:" in line:
                single = line.strip().split()[-1]
                print(single)
            if "Duplicated:" in line:
                dup = line.rstrip().split()[-1]
                print(dup)
            if "Duplicated, Unexpected:" in line:
                dup_un = line.strip().split()[-1]
                print(dup_un)
            if "Duplicated, Expected:" in line:
                dup_exp = line.strip().split()[-1]
                print(dup_exp)
            if "Missing:" in line:
                missing = line.strip().split()[-1]
                print(missing)
            if "Total Consistent" in line:
                consistent = line.strip().split()[-1]
                print(consistent)
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
            if "Contaminants, partial hits":
                contaminants_partial = line.strip().split()[-1]
            if "Contaminants, fragmented": 
                contaminants_fragmented = line.strip().split()[-1]
            if "Total Unknown" in line:
                unkown = line.strip().split()[-1]
            if "Clade" in line:
                clade = line.split(":")[-1].strip()
            if "associated query proteins" in line:
                clade_per = line.strip().split()[-1]
                clades.append("{}: {}".format(clade, clade_per))
    results["OMArk Consistency Results"] = "Cons:{}[P:{};F:{}],Inco:{}[P:{},F:{}],Cont:{},Unkn:{}".format(consistent,
                                                                                                          cons_partial,
                                                                                                          cons_frag,
                                                                                                          inconsistent,
                                                                                                          inconsistent_partial,
                                                                                                          inconsistent_fragmented,
                                                                                                          contaminants,
                                                                                                          unkown)
    results["OMArk Completeness Results"] = "{}: {}; S:{},D:{}[U:{},E:{}],M:{}".format(clade, hogs, single,
                                                                                             dup, dup_un, dup_exp,
                                                                                             missing)
    results["OMArk Species Composition"] = "; ".join(clades)
    return results
            
