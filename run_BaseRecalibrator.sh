#BSUB -o Myjob.%J.log
#BSUB -e Myjob.%J.err
#BSUB -n 4 -R "rusage[mem=16]"
#BSUB -W 48:00

# Generate Base Quality Score Recalibration (BQSR) model

module load java/1.8.0_31

compression_level="1"
java_heap_memory_initial="4g"
tool_path="/data/ldiaz/luol2/GATK_workflows/up-to-date-broad-genomes-in-the-cloud-tools"

ref_fasta="/data/ldiaz/luol2/References/GRCm38/GRCm38.fa"
ref="/data/ldiaz/luol2/References/GRCm38"
project=PROJECT
subject=SUBJECT
sample=TEST_SAMPLE

#data_path="/home/luol2/lingqi_workspace/Projects/Jones_Lee_Exercise_Project"
data_path=DATA_PATH

input_bam="${data_path}/${project}/${subject}/${sample}/${sample}.aligned.duplicate_marked.sorted.bam"
input_bai="${data_path}/${project}/${subject}/${sample}/${sample}.aligned.duplicate_marked.sorted.bai"
recalibration_report_filename="${data_path}/${project}/${subject}/${sample}/${sample}.recal_data.csv"
#dbSNP_vcf="/home/luol2/lingqi_workspace/References/GRCm38/Mus_musculus/snpdb/mouse_dbsnp.sort.vcf.gz"
#dbMGP_indels_vcf="/home/luol2/lingqi_workspace/References/GRCm38/Mus_musculus/MGP_indels/mgp.v5.indels.pass.sort.vcf.gz"
#dbENSEMBL_vcf="/home/luol2/lingqi_workspace/References/GRCm38/Mus_musculus/Ensembl/known_germline_variants/mus_musculus.vcf.gz"
dbMGP_snp_indels_vcf="/data/ldiaz/luol2/References/MGP_snp_indel_db/strain_specific/BALB_cJ.mgp.v5.combined.snps_indels.dbSNP142.vcf"

${tool_path}/gatk4/gatk-launch --javaOptions "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
      -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
      -Xloggc:gc_log.log -Dsamjdk.compression_level=${compression_level} -Xmx${java_heap_memory_initial}" \
      BaseRecalibrator \
      -R ${ref_fasta} \
      -I ${input_bam} \
      --useOriginalQualities \
      -O ${recalibration_report_filename} \
      -knownSites ${dbMGP_snp_indels_vcf} \
      -L ${ref}/mouse_all_exon_mm10.interval_list


