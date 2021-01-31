#BSUB -o Myjob.%J.log
#BSUB -e Myjob.%J.err
#BSUB -n 1 -R "rusage[mem=8]"
#BSUB -W 48:00

module load singularity/3.6.0
#module load java/1.8.0_31

path=PATH
manifest=MANIFEST

#docker run -v ${PWD}:/data marsluo/deconstructsigshg38 Rscript /data/deconstructSigs.R /data 2>&1 > mutsig.log
singularity run \
    --bind ${path}:/data \
    /home/luol2/lingqi_workspace/deconstructsigshg38.sif \
           Rscript /data/data_analysis/deconstructSigs.R /data ${manifest} 2>&1 > mutsig.log
