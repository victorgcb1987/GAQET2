from error_check import error_check

def psauron_stats(psauron):
    error = error_check(psauron)
    if error_check:
        return {"PSAURON SCORE": error}
    with open(psauron["outfile"]) as fhand:
        for line in fhand:
            if "psauron score" in line:
                return {"PSAURON SCORE": line.strip().split()[-1]}        