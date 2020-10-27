#BSUB -o Myjob.%J.log
#BSUB -e Myjob.%J.err
#BSUB -n 1 -R "rusage[mem=8]"
#BSUB -W 48:00

module load singularity/3.1.1
project=PROJECT
subject=SUBJECT
sample=TEST_SAMPLE

#data_path="/home/luol2/lingqi_workspace/Projects/Jones_Lee_Exercise_Project"
data_path=DATA_PATH

# Sample_Folder would be replacable by any folder path to a sample
mkdir -p ${data_path}/${project}/${subject}/${sample}/trim_galore_output;
/bin/cp /home/luol2/bin/trim_galore_docker_code.sh ${data_path}/${project}/${subject}/${sample}/trim_galore_output;
singularity run \
    --bind ${data_path}/${project}/${subject}/${sample}:/data \
    --bind ${data_path}/${project}/${subject}/${sample}/trim_galore_output:/data2 \
    /data/ldiaz/luol2/trim_galore_v0.6.0.sif /bin/bash /data2/trim_galore_docker_code.sh

