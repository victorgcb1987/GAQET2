import subprocess


def operation_failed(results):
    if "Failed" in results["status"]:
        return "FAILED"
    else:
        return False
    
def correct_fasta_length(arguments):
    cmd = "sed -n '2p' {} | wc -c".format(arguments["Assembly"])
    seq_length = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
    print(seq_length.stdout.decode().strip())
