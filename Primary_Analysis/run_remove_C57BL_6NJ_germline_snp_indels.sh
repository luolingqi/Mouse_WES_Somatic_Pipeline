#BSUB -o Myjob.%J.log
#BSUB -e Myjob.%J.err
#BSUB -n 2 -R "rusage[mem=8]"
#BSUB -W 48:00

module load java/1.8.0_31
ref_fasta="/data/ldiaz/luol2/References/GRCm38/GRCm38.fa"
ref_MGP_snps_indels="/data/ldiaz/luol2/References/MGP_snp_indel_db/strain_specific/C57BL_6NJ.mgp.v5.combined.snps_indels.dbSNP142.vcf"

project=PROJECT
subject=SUBJECT
sample=TEST_SAMPLE

data_path=DATA_PATH

unfiltered_vcf=${data_path}/${project}/${subject}/${sample}/${sample}.aligned.duplicates_marked.recalibrated.filtered.vcf # actually mutect2 filtered vcf
command_mem="8000"

#/data/ldiaz/luol2/gatk-4.0.1.2/gatk  --java-options "-Xmx${command_mem}m" SelectVariants -R ${ref_fasta} -V ${unfiltered_vcf} --discordance ${ref_MGP_snps_indels} -O ${unfiltered_vcf/.vcf/.w.BALCnormal.rm_MGP_snps_indels.vcf}
/data/ldiaz/luol2/gatk-4.1.2.0/gatk --java-options "-Xmx${command_mem}m" SelectVariants -R ${ref_fasta} -V ${unfiltered_vcf} --discordance ${ref_MGP_snps_indels} -O ${unfiltered_vcf/.vcf/.w.C57BL_6NJ.normal.rm_MGP_snps_indels.vcf}
