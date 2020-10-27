#BSUB -o Myjob.%J.log
#BSUB -e Myjob.%J.err
#BSUB -n 2 -R "rusage[mem=8]"
#BSUB -W 48:00

module load java/1.8.0_31
ref_fasta="/data/ldiaz/luol2/References/GRCm38/GRCm38.fa"
ref_MGP_snps_indels="/data/ldiaz/luol2/References/MGP_snp_indel_db/strain_specific/BALB_cJ.mgp.v5.combined.snps_indels.dbSNP142.vcf"
#ref_dbsnp="/data/ldiaz/luol2/References/GRCm38/Mus_musculus/snpdb/mouse_dbsnp.sort.vcf.gz"
#ref_MGP_indels="/data/ldiaz/luol2/References/GRCm38/Mus_musculus/MGP_indels/mgp.v5.indels.pass.sort.vcf.gz"
#ref_Ensembl_variant="/data/ldiaz/luol2/References/GRCm38/Mus_musculus/Ensembl/known_germline_variants/mus_musculus.vcf.gz"

project=PROJECT
subject=SUBJECT
sample=TEST_SAMPLE

#data_path="/data/ldiaz/luol2/Projects/Jones_Lee_Exercise_Project"
data_path=DATA_PATH

#unfiltered_vcf=${data_path}/${sample}/${sample}.aligned.duplicates_marked.recalibrated.unfiltered.vcf
unfiltered_vcf=${data_path}/${project}/${subject}/${sample}/${sample}.aligned.duplicates_marked.recalibrated.filtered.vcf # actually mutect2 filtered vcf
command_mem="8000"

#ln -s /home/luol2/lingqi_workspace/gatk-4.0.1.2/gatk-package-4.0.1.2-local.jar gatk-package-4.0.1.2-local.jar
#java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=1 "-Xmx${command_mem}m" gatk-package-4.0.1.2-local.jar  SelectVariants -R ${ref_fasta} -V ${unfiltered_vcf} --discordance ${ref_MGP_snps_indels} -O ${unfiltered_vcf/.vcf/.w.BALCnormal.rm_MGP_snps_indels.vcf}

#/data/ldiaz/luol2/gatk-4.0.1.2/gatk  --java-options "-Xmx${command_mem}m" SelectVariants -R ${ref_fasta} -V ${unfiltered_vcf} --discordance ${ref_MGP_snps_indels} -O ${unfiltered_vcf/.vcf/.w.BALCnormal.rm_MGP_snps_indels.vcf}
/data/ldiaz/luol2/gatk-4.1.2.0/gatk --java-options "-Xmx${command_mem}m" SelectVariants -R ${ref_fasta} -V ${unfiltered_vcf} --discordance ${ref_MGP_snps_indels} -O ${unfiltered_vcf/.vcf/.w.BALCnormal.rm_MGP_snps_indels.vcf}
