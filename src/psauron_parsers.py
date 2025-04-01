def psauron_stats(psauron):
    with open(psauron["outfile"]) as fhand:
        for line in fhand:
            if "psauron score" in line:
                return {"PSAURON SCORE": line.strip().split()[-1]}        