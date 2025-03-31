def busco_stats(busco):
    results = {}
    for lineage, results in busco.items():
        with open(results["outfile"]) as fhand:
            for line in fhand:
                if "%" in line:
                    results[lineage] = line.strip()
    return results