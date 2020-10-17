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
#data_path="/home/luol2/lingqi_workspace/Projects/Ben_Projects/WES_mouse"

project=PROJECT
subject=SUBJECT
sample=TEST_SAMPLE

#data_path="/home/luol2/lingqi_workspace/Projects/Jones_Lee_Exercise_Project"
data_path=DATA_PATH

input_bam="${data_path}/${project}/${subject}/${sample}/${sample}.aligned.duplicate_marked.sorted.bam"
input_bai="${data_path}/${project}/${subject}/${sample}/${sample}.aligned.duplicate_marked.sorted.bai"
recalibration_report="${data_path}/${project}/${subject}/${sample}/${sample}.recal_data.csv"

output_bam_basename="${data_path}/${project}/${subject}/${sample}/${sample}.aligned.duplicates_marked.recalibrated"

${tool_path}/gatk4/gatk-launch --javaOptions "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
      -XX:+PrintGCDetails -Xloggc:gc_log.log \
      -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Dsamjdk.compression_level=${compression_level} -Xmx${java_heap_memory_initial}" \
      ApplyBQSR \
      --createOutputBamMD5 \
      --addOutputSAMProgramRecord \
      -R ${ref_fasta} \
      -I ${input_bam} \
      --useOriginalQualities \
      -O ${output_bam_basename}.bam \
      -bqsr ${recalibration_report} \
      -SQQ 10 -SQQ 20 -SQQ 30 \
      -L ${ref}/mouse_all_exon_mm10.interval_list

