###########################################################################################################################################################
# The script prepares multiple manifest files based on experiment design from researchers
# It takes a pairwise comparison file as the input, and output the following two manifest files for VAF and Mutsig plotting respectively:
#       * items_to_plots.txt
#       * items_to_plots_MutSig.txt
###########################################################################################################################################################
import sys

comparison = sys.argv[1]
sequencer = sys.argv[2]
output_vaf = sys.argv[3]
output_sig = sys.argv[4]

parental = ""
treated = []


with open(comparison,"r") as fin:
    for line in fin:
        line = line.strip()
        arr = line.split()
        if parental == "":
            parental = arr[0]
        treated.append(arr[1])


# Prepare the items_to_plot
change = ["gain","loss","total"]
type = ["nonsynonymous","synonymous","frameshift","missense"]

with open(output_vaf, "w") as vaf:
    vaf.write("cats\tcomparison\tvariant_type\tfiles\n")
    for sample in treated:
        for c in change:
            for t in type:
                file = ""
                if c == "gain":
                    file = "\t".join([sample,c,t,sequencer+"/"+parental+"_vs_"+sample+".VEP.ann.altAD4_AF0.01_MBQ20_MMQ50_MPOS5."+ t + "_variant_only.query"])
                elif c == "loss":
                    file = "\t".join([sample,c,t,sequencer+"/"+sample+"_vs_"+parental+".VEP.ann.altAD4_AF0.01_MBQ20_MMQ50_MPOS5."+ t + "_variant_only.query"])
                elif c == "total":
                    file = "\t".join([sample,c,t,sequencer+"/"+sample+"/"+sample+".aligned.duplicates_marked.recalibrated.filtered.w.BALCnormal.rm_MGP_snps_indels.singleSample.PASSOnly.VEP.ann.altAD4_AF0.01_MBQ20_MMQ50_MPOS5."+t+"_variant_only.query"])
                vaf.write(file+"\n") 
    

# Prepare the items_to_plot_MutSig
change = ["gain","loss","total"]

with open(output_sig, "w") as sig:
    sig.write("cats\tcomparison\tfiles\n")
    for sample in treated:
        for c in change:
            file = ""
            if c == "gain":
                file = "\t".join([sample,c,sequencer+"/"+parental+"_vs_"+sample+".VEP.ann.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.vcf"])
            elif c == "loss":
                file = "\t".join([sample,c,sequencer+"/"+sample+"_vs_"+parental+".VEP.ann.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.vcf"])
            elif c == "total":
                file = "\t".join([sample,c,sequencer+"/"+sample+"/"+sample+".aligned.duplicates_marked.recalibrated.filtered.w.BALCnormal.rm_MGP_snps_indels.singleSample.PASSOnly.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.VEP.ann.vcf"])
            sig.write(file+"\n") 
