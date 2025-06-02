from error_check import error_check

def protein_homology_stats(homology, num_transcripts):
    results = {}
    for tag, values in homology.items():
        prots = set()
        error = error_check(values)
        if error:
            results[f"ProteinsWith{tag}Hits (%)"] = error
        else:
            with open(values["outfile"]) as results_fhand:
                for line in results_fhand:
                    parts = line.rstrip().split()
                    if len(parts) > 10:  # para evitar errores por líneas mal formateadas
                        try:
                            evalue = float(parts[10])
                            prot_id = parts[0]
                            if evalue < 1e-20:
                                prots.add(prot_id)
                        except ValueError:
                            continue  # en caso de que evalue no sea convertible
            percentage = round((len(prots) / num_transcripts) * 100, 2)
            results[f"ProteinsWith{tag}Hits (%)"] = percentage
    return results