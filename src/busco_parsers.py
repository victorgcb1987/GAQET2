from error_check import error_check

def busco_stats(busco):
    error = error_check(busco)
    if error:
        return error
    results = {}
    for lineage, report in busco.items():
        with open(report["outfile"]) as fhand:
            for line in fhand:
                if "%" in line:
                    results[lineage] = line.strip()
    return results