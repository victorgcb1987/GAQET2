def error_check(results):
    if "Failed" in results["status"]:
        return "FAILED"
    else:
        return False