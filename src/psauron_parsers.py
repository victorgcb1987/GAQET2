from src.error_check import operation_failed

def psauron_stats(psauron):
    print(psauron)
    error = operation_failed(psauron)
    if operation_failed:
        return {"PSAURON SCORE": error}
    print(psauron["outfile"])
    with open(psauron["outfile"]) as fhand:
        for line in fhand:
            if "psauron score" in line:
                return {"PSAURON SCORE": line.strip().split()[-1]}        