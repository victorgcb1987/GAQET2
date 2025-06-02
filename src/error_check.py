def operation_failed(results):
    if "Failed" in results["status"]:
        return "FAILED"
    else:
        return False