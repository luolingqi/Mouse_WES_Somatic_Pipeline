#BSUB -o Myjob.%J.log
#BSUB -e Myjob.%J.err
#BSUB -n 2 -R "rusage[mem=8]"
#BSUB -W 48:00

module load singularity/3.0.1-to-ve-removed
module load java/1.8.0_31

project=PROJECT
subject=SUBJECT
sample=TEST_SAMPLE

#data_path="/data/ldiaz/luol2/Projects/Jones_Lee_Exercise_Project"
data_path=DATA_PATH
#ref_fasta="/data/ldiaz/luol2/References/GRCm38/GRCm38.fa"
ref_fasta="/data/ldiaz/luol2/References/Mus_musculus_balbcj.BALB_cJ_v1.dna.toplevel.fa"

# frameshift mutation as input
InputVcf="${data_path}/${project}/${subject}/${sample}/${sample}.aligned.duplicates_marked.recalibrated.filtered.w.BALCnormal.rm_MGP_snps_indels.singleSample.PASSOnly.VEP.ann.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.frameshift_variant_only.vcf"
OutputVcf=${InputVcf/.vcf/.LiftOver_to_BALB_cJ.vcf}
RejectVcf=${InputVcf/.vcf/.rejected_variants.vcf}
chain="/data/ldiaz/luol2/References/mm10-to-BALB_cJ.quality.reformatted.chain"

# modify chromesome names for X, Y, MT to 20, 21, 22
cp ${InputVcf} ${InputVcf/frameshift_variant_only.vcf/frameshift_variant_only.chomXYMT_modified.vcf}
perl -i -pe 's/contig=<ID=X/contig=<ID=20/g;s/^X/20/g' ${InputVcf/frameshift_variant_only.vcf/frameshift_variant_only.chomXYMT_modified.vcf}
perl -i -pe 's/contig=<ID=Y/contig=<ID=21/g;s/^Y/21/g' ${InputVcf/frameshift_variant_only.vcf/frameshift_variant_only.chomXYMT_modified.vcf}
perl -i -pe 's/contig=<ID=MT/contig=<ID=22/g;s/^MT/22/g' ${InputVcf/frameshift_variant_only.vcf/frameshift_variant_only.chomXYMT_modified.vcf}


java -jar /data/ldiaz/luol2/picard.jar LiftoverVcf \
I=${InputVcf/frameshift_variant_only.vcf/frameshift_variant_only.chomXYMT_modified.vcf} \
O=${OutputVcf} \
CHAIN=${chain} \
REJECT=${RejectVcf} \
R=${ref_fasta}
