#BSUB -o Myjob.%J.log
#BSUB -e Myjob.%J.err
#BSUB -n 4 -R "rusage[mem=16]"
#BSUB -W 48:00

module load java/1.8.0_31
module load samtools/1.7

project=PROJECT
subject=SUBJECT
sample=TEST_SAMPLE

data_path="/home/luol2/lingqi_workspace/Projects/Ben_Projects/WES_mouse_01062020"

prob=PROB
#cov=COV

input_bam="${data_path}/${project}/${subject}/${sample}/${sample}.aligned.duplicates_marked.recalibrated.bam"
output_bam="${data_path}/${project}/${subject}/${sample}/${sample}.aligned.duplicates_marked.recalibrated_${prob}X.bam"

/home/luol2/lingqi_workspace/gatk-4.0.1.2/gatk DownsampleSam \
    -I ${input_bam} \
    -O ${output_bam} \
    -P ${prob}

samtools index ${output_bam}
