#BSUB -o Myjob.%J.log
#BSUB -e Myjob.%J.err
#BSUB -n 16 -R "rusage[mem=8]"
#BSUB -W 48:00

module load java/1.8.0_31

compression_level="1"
java_heap_memory_initial="2g"
mem_limit="2G"
samtools_threads="32"

tool_path="/data/ldiaz/luol2/GATK_workflows/up-to-date-broad-genomes-in-the-cloud-tools"

project=PROJECT
subject=SUBJECT
sample=TEST_SAMPLE

#data_path="/home/luol2/lingqi_workspace/Projects/Jones_Lee_Exercise_Project"
data_path=DATA_PATH

INPUT="${data_path}/${project}/${subject}/${sample}/${sample}.bam"
OUTPUT="${data_path}/${project}/${subject}/${sample}/${sample}.aligned.unsorted.duplicates_marked"
METRICS="${data_path}/${project}/${subject}/${sample}/${sample}.duplicate_metrics"

SORT_BAM="${data_path}/${project}/${subject}/${sample}/${sample}.aligned.duplicate_marked.sorted.bam"
SORT_BAI="${data_path}/${project}/${subject}/${sample}/${sample}.aligned.duplicate_marked.sorted.bai"

java -Dsamjdk.compression_level=${compression_level} -Xmx${java_heap_memory_initial} -jar ${tool_path}/picard.jar \
    MarkDuplicates \
    INPUT=${INPUT} \
    OUTPUT= ${OUTPUT} \
    METRICS_FILE= ${METRICS} \
    VALIDATION_STRINGENCY=SILENT \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
    ASSUME_SORT_ORDER="queryname" \
    CLEAR_DT="false" \
    ADD_PG_TAG_TO_READS=false

# Sort the deduplicated sam files
${tool_path}/samtools sort -m ${mem_limit} --threads ${samtools_threads} -l ${compression_level} ${OUTPUT} -o ${SORT_BAM} 
${tool_path}/samtools index ${SORT_BAM} ${SORT_BAI}


