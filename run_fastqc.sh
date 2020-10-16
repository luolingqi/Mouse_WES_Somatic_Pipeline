#BSUB -o Myjob.%J.log
#BSUB -e Myjob.%J.err
#BSUB -n 1 -R "rusage[mem=16]"
#BSUB -W 48:00

module load singularity/3.0.1-to-ve-removed

trim=TRIM
project=PROJECT
subject=SUBJECT
sample=TEST_SAMPLE

#data_path="/home/luol2/lingqi_workspace/Projects/Jones_Lee_Exercise_Project"
data_path=DATA_PATH

# Sample_Folder would be replacable by any folder path to a sample

fastq1=$(basename ${data_path}/${project}/${subject}/${sample}/*R1_001${trim}.fastq.gz)
fastq2=$(basename ${data_path}/${project}/${subject}/${sample}/*R2_001${trim}.fastq.gz)

singularity run --bind ${data_path}/${project}/${subject}/${sample}:/data /data/ldiaz/luol2/trim_galore_v0.6.0.sif fastqc /data/$fastq1 /data/$fastq2 -o /data



