from error_check import operation_failed

def busco_stats(busco):
    error = operation_failed(busco)
    if error:
        return error
    results = {}
    for lineage, report in busco.items():
        with open(report["outfile"]) as fhand:
            for line in fhand:
                if "%" in line:
                    results[lineage] = line.strip()
    return results