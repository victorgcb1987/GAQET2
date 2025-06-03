from src.error_check import operation_failed

def psauron_stats(psauron):
    error = operation_failed(psauron)
    if error:
        return {"PSAURON SCORE": error}
    with open(psauron["outfile"]) as fhand:
        for line in fhand:
            if "psauron score" in line:
                return {"PSAURON SCORE": line.strip().split()[-1]}   
    return  {"PSAURON SCORE": "FAILED"}  