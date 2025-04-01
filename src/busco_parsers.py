def busco_stats(busco):
    results = {}
    for lineage, report in busco.items():
        print(lineage)
        with open(report["outfile"]) as fhand:
            for line in fhand:
                if "%" in line:
                    results[lineage] = line.strip()
    return results