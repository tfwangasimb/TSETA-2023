import sys

if __name__ != '__main__':
    sys.exit(1)

# comp = ["Rad51_vs_Sae2","spo11_rad51_vs_Sae2","spo11_rad51_vs_spo11_sae2","spo11_sae2_vs_Sae2"]
# cutoff = [2, 3, 4]

comp = sys.argv[1].split(",")
cutoff = [int(i) for i in sys.argv[2].split(",")]

if len(comp) == 0 or len(cutoff) == 0:
    sys.exit(1)

for cp in comp:
    for cf in cutoff:
        start = 0
        frag_start = 0
        frag_ct = 0
        fct = 0
        fsum = 0
        fhigh = 0
        fname1 = "Diff_peaks_lg"+ str(cf) + "_" + cp + ".txt"
        fname2 = "Diff_peak_bins_lg"+ str(cf) + "_" + cp + ".txt"
        with open(fname1, 'w') as fout, open(fname2, 'w') as fout2:
            fout.write("Peak\tChr\tStart\tEnd\tMean\tHighest\n")
            fout2.write("Peak\tChr\tStart\tEnd\tValue\n")
            fname3 = "diff_bins_scaled_"+ cp +".tsv"
            with open(fname3, 'r') as f:
                for line in f:
                    ln = line.rstrip('\n')
                    ar = ln.split("\t")
                    if start == 0 and float(ar[3]) > cf: ## start of peak
                        frag_start = ar[1]
                        start = 1
                        fsum += float(ar[3])
                        fhigh = float(ar[3])
                        fct += 1
                        frag_ct += 1
                        frag_name = "peak" + str(frag_ct)
                        fout2.write(frag_name + "\t" + ar[0] + "\t" + str(ar[1]) + "\t" + str(ar[2]) + "\t" + str(ar[3]) + "\n")
                    elif start == 1 and float(ar[3]) <= cf: ## end of peak
                        fout.write(frag_name + "\t" + ar[0] + "\t" + str(frag_start) + "\t" + str(ar[1]) + "\t" + str(round(fsum/fct, 2)) + "\t" + str(fhigh) + "\n")
                        frag_start = 0
                        fct = 0
                        fsum = 0
                        fhigh = 0
                        start = 0
                    elif start == 1 and float(ar[3]) > cf:    ## continue of peak
                        fsum += float(ar[3])
                        fct += 1
                        if float(ar[3]) > float(fhigh):
                            fhigh = ar[3]
                        fout2.write(frag_name + "\t" + ar[0] + "\t" + str(ar[1]) + "\t" + str(ar[2]) + "\t" + str(ar[3]) + "\n")
            with open(fname3, 'r') as f:
                for line in f:
                    ln = line.rstrip('\n')
                    ar = ln.split("\t")
                    if start == 0 and float(ar[3]) < -cf: ## start of peak
                        frag_start = ar[1]
                        start = 1
                        fsum += float(ar[3])
                        fhigh = float(ar[3])
                        fct += 1
                        frag_ct += 1
                        frag_name = "peak" + str(frag_ct)
                        fout2.write(frag_name + "\t" + ar[0] + "\t" + str(ar[1]) + "\t" + str(ar[2]) + "\t" + str(ar[3]) + "\n")
                    elif start == 1 and float(ar[3]) >= -cf: ## end of peak
                        fout.write(frag_name + "\t" + ar[0] + "\t" + str(frag_start) + "\t" + str(ar[1]) + "\t" + str(round(fsum/fct, 2)) + "\t" + str(fhigh) + "\n")
                        frag_start = 0
                        fct = 0
                        fsum = 0
                        fhigh = 0
                        start = 0
                    elif start == 1 and float(ar[3]) < -cf:    ## continue of peak
                        fsum += float(ar[3])
                        fct += 1
                        if float(ar[3]) < float(fhigh):
                            fhigh = ar[3]
                        fout2.write(frag_name + "\t" + ar[0] + "\t" + str(ar[1]) + "\t" + str(ar[2]) + "\t" + str(ar[3]) + "\n")
        fout.close()
        fout2.close()
        f.close()
        print(fname3)
