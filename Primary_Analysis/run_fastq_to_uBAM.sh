#BSUB -o Myjob.%J.log
#BSUB -e Myjob.%J.err
#BSUB -n 1 -R "rusage[mem=32]"
#BSUB -W 48:00

module load java/1.8.0_31

project=PROJECT
subject=SUBJECT
sample=TEST_SAMPLE
readGroup=READGROUP # e.g.   H2JFYBBXY.5
data_path=DATA_PATH #e.g.    /home/luol2/lingqi_workspace/Projects/Jones_Lee_Exercise_Project
mkdir `pwd`/tmp_${sample}
#R1_fa=$(ls ${data_path}/${project}/${subject}/${sample}/${sample/Sample_/}_S*_R1_001_trim.fastq.gz)
#R2_fa=$(ls ${data_path}/${project}/${subject}/${sample}/${sample/Sample_/}_S*_R2_001_trim.fastq.gz)

R1_fa=$(ls ${data_path}/${project}/${subject}/${sample}/${sample/Sample_/}_S*_R1_001.fastq.gz)
R2_fa=$(ls ${data_path}/${project}/${subject}/${sample}/${sample/Sample_/}_S*_R2_001.fastq.gz)

java -Djava.io.tmpdir=`pwd`/tmp_${sample} -jar /home/luol2/lingqi_workspace/picard.jar FastqToSam \
      F1=${R1_fa} \
      F2=${R2_fa} \
      O=${data_path}/${project}/${subject}/${sample}/${sample}.unmapped.bam \
      SAMPLE_NAME=${sample} \
      READ_GROUP_NAME=${readGroup} \
      LIBRARY_NAME=Lib.${sample} \
      PLATFORM_UNIT=${readGroup} \
      PLATFORM=illumina


