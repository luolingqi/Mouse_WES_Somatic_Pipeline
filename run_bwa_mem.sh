#BSUB -o Myjob.%J.log
#BSUB -e Myjob.%J.err
#BSUB -n 8 -R "rusage[mem=8]"
#BSUB -W 48:00

module load java/1.8.0_31

bash_ref_fasta="/data/ldiaz/luol2/References/GRCm38"
bwa_threads=16
tool_path="/data/ldiaz/luol2/GATK_workflows/up-to-date-broad-genomes-in-the-cloud-tools"
compression_level="1"
java_heap_memory_initial="6g"

project=PROJECT
subject=SUBJECT
sample=TEST_SAMPLE

#data_path="/home/luol2/lingqi_workspace/Projects/Jones_Lee_Exercise_Project"
data_path=DATA_PATH #e.g.    /home/luol2/lingqi_workspace/Projects/Jones_Lee_Exercise_Project

input_bam=$(ls ${data_path}/${project}/${subject}/${sample}/*unmapped.bam)
#input_bam=${data_path}/${sample}/CT26_P_3m_IGO_10212_1_S137.unmapped.bam  # original unmapped BAM file
bwa_commandline="bwa mem -K 100000000 -p -v 3 -t $bwa_threads -Y ${bash_ref_fasta}/GRCm38.fa"
#bwa_commandline="bwa mem -K 100000000 -p -v 3 -t $bwa_threads -Y ${bash_ref_fasta}/Mus_musculus_balbcj.BALB_cJ_v1.dna.toplevel.fa"

java -Dsamjdk.compression_level=${compression_level} -Xmx${java_heap_memory_initial} -jar ${tool_path}/picard.jar \
        SamToFastq \
        VALIDATION_STRINGENCY=SILENT \
        INPUT=${input_bam} \
        FASTQ=/dev/stdout \
        INTERLEAVE=true \
        NON_PF=true | \
      ${tool_path}/${bwa_commandline} /dev/stdin - | \
 java -Dsamjdk.compression_level=${compression_level} -Xmx${java_heap_memory_initial} -jar ${tool_path}/picard.jar \
        MergeBamAlignment \
        VALIDATION_STRINGENCY=SILENT \
        EXPECTED_ORIENTATIONS=FR \
        ATTRIBUTES_TO_RETAIN=X0 \
        ATTRIBUTES_TO_REMOVE=NM \
        ATTRIBUTES_TO_REMOVE=MD \
        ALIGNED_BAM=/dev/stdin \
        UNMAPPED_BAM=${input_bam} \
        OUTPUT=${data_path}/${project}/${subject}/${sample}/${sample}.bam \
        REFERENCE_SEQUENCE=${bash_ref_fasta}/GRCm38.fa \
        PAIRED_RUN=true \
        SORT_ORDER="unsorted" \
        IS_BISULFITE_SEQUENCE=false \
        ALIGNED_READS_ONLY=false \
        CLIP_ADAPTERS=false \
        MAX_RECORDS_IN_RAM=2000000 \
        ADD_MATE_CIGAR=true \
        MAX_INSERTIONS_OR_DELETIONS=-1 \
        PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
        PROGRAM_RECORD_ID="bwamem" \
        PROGRAM_GROUP_VERSION="0.7.15-r1140" \
        PROGRAM_GROUP_COMMAND_LINE="bwa mem -K 100000000 -p -v 3 -t $bwa_threads -Y ${bash_ref_fasta}/GRCm38.fa" \
        PROGRAM_GROUP_NAME="bwamem" \
        UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
        ALIGNER_PROPER_PAIR_FLAGS=true \
        UNMAP_CONTAMINANT_READS=true \
        ADD_PG_TAG_TO_READS=false






