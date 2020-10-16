#BSUB -o Myjob.%J.log
#BSUB -e Myjob.%J.err
#BSUB -n 8 -R "rusage[mem=8]"
#BSUB -W 48:00

module load samtools/1.10

project1=PROJ1
project2=PROJ2
sample=SAMPLE
data_path="/home/luol2/lingqi_workspace/Projects/Ben_Projects/WES_mouse_Project_10212_E"

samtools merge ${data_path}/${sample}.merged.bam \
${data_path}/${project1}/${sample}/${sample}.aligned.duplicates_marked.recalibrated.bam \
${data_path}/${project2}/${sample}/${sample}.aligned.duplicates_marked.recalibrated.bam

samtools index ${data_path}/${sample}.merged.bam
