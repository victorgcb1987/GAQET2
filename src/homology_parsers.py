def protein_homology_stats(homology, num_transcripts):
    results = {}
    for tag, values in homology.items():
        prots = []
        with open(values["outfile"]) as results_fhand:
            for line in results_fhand:
                line = line.rstrip().split()
                evalue = float(line[10])
                prot_id = line[0]
                if evalue < 1e-20 and prot_id not in prots:
                    prots.append(line[0])
        results["ProteinsWith{}Hits (%)".format(tag)] = round((float(len(prots)/num_transcripts) * 100), 2)
    return results