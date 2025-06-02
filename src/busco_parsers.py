from src.error_check import operation_failed

def busco_stats(busco):
    results = {}
    for lineage, report in busco.items():
        error = operation_failed(report)
        if error:
             results[lineage] = error
        else:
            with open(report["outfile"]) as fhand:
                for line in fhand:
                    if "%" in line:
                        results[lineage] = line.strip()
    return results