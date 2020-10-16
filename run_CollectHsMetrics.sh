#BSUB -o Myjob.%J.log
#BSUB -e Myjob.%J.err
#BSUB -n 4 -R "rusage[mem=16]"
#BSUB -W 48:00

module load java/1.8.0_31
module load R/R-3.6.0

compression_level="1"
java_heap_memory_initial="5g"
tool_path="/data/ldiaz/luol2/GATK_workflows/up-to-date-broad-genomes-in-the-cloud-tools"

ref_fasta="/data/ldiaz/luol2/References/GRCm38/GRCm38.fa"
ref="/data/ldiaz/luol2/References/GRCm38"


project=PROJECT
subject=SUBJECT
sample=TEST_SAMPLE

#data_path="/home/luol2/lingqi_workspace/Projects/Jones_Lee_Exercise_Project"
data_path=DATA_PATH

input_bam="${data_path}/${project}/${subject}/${sample}/sorted.${sample}.bam"
input_bai="${data_path}/${project}/${subject}/${sample}/sorted.${sample}.bai"
output_bam_prefix="${data_path}/${project}/${subject}/${sample}/sorted.${sample}.aligned"

java -Xmx${java_heap_memory_initial} -jar ${tool_path}/picard.jar \
    CollectHsMetrics \
    INPUT=${input_bam} \
    OUTPUT=${output_bam_prefix}.hs_metrics.txt \
    R=${ref_fasta} \
    BAIT_INTERVALS=${ref}/mouse_all_exon_mm10.interval_list \
    TARGET_INTERVALS=${ref}/mouse_all_exon_mm10.interval_list


