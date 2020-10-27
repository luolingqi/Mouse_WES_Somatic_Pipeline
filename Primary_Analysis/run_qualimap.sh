#BSUB -o Myjob.%J.log
#BSUB -e Myjob.%J.err
#BSUB -n 16 -R "rusage[mem=8]"
#BSUB -W 48:00

module load singularity/3.1.1
module load samtools/1.7

#bn=$(basename "Sample_Folder")

project=PROJECT
subject=SUBJECT
sample=TEST_SAMPLE

#data_path="/home/luol2/lingqi_workspace/Projects/Jones_Lee_Exercise_Project"
data_path=DATA_PATH

# Sample_Folder would be replacable by any folder path to a sample
mkdir -p ${data_path}/${project}/${subject}/${sample}/qualimap
mem_limit=2G
samtools_threads=32
compression_level=1

# sort bam
samtools sort -m ${mem_limit} --threads ${samtools_threads} -l ${compression_level} -o ${data_path}/${project}/${subject}/${sample}/sorted.${sample}.bam ${data_path}/${project}/${subject}/${sample}/${sample}.bam
samtools index ${data_path}/${project}/${subject}/${sample}/sorted.${sample}.bam ${data_path}/${project}/${subject}/${sample}/sorted.${sample}.bai

#qualimap bamqc
singularity run --bind ${data_path}/${project}/${subject}/${sample}/:/data \
                --bind ${data_path}/${project}/${subject}/${sample}/qualimap:/data2 \
                --bind /data/ldiaz/luol2/References/GRCm38:/ref \
                /data/ldiaz/luol2/qualimap_v2.2.1.sif qualimap bamqc \
                        -bam /data/sorted.${sample}.bam \
                        --feature-file /ref/nochr_S0276129_Regions_mm10.gff \
                        --java-mem-size=64G \
                        -outfile ${sample}.exome.pdf \
                        -outdir /data2

singularity run --bind ${data_path}/${project}/${subject}/${sample}/:/data \
                --bind ${data_path}/${project}/${subject}/${sample}/qualimap:/data2 \
                --bind /data/ldiaz/luol2/References/GRCm38:/ref \
                /data/ldiaz/luol2/qualimap_v2.2.1.sif qualimap rnaseq \
                        -bam /data/sorted.${sample}.bam \
                        -gtf /ref/Mus_musculus/Ensembl/gtf/Mus_musculus.GRCm38.97.gtf \
                        --java-mem-size=64G \
                        -outfile ${sample}.rnaseq.pdf \
                        -outdir /data2




