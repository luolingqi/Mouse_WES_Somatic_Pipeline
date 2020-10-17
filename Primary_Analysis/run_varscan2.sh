#BSUB -o Myjob.%J.log
#BSUB -e Myjob.%J.err
#BSUB -n 8 -R "rusage[mem=8]"
#BSUB -W 48:00

module load samtools/1.7
module load java/1.8.0_31
bash_ref_fasta="/home/luol2/lingqi_workspace/References/GRCm38"
tool_path="/home/luol2/lingqi_workspace/GATK_workflows/up-to-date-broad-genomes-in-the-cloud-tools"

data_path="/home/luol2/lingqi_workspace/Projects/Ben_Projects/WES_mouse"

normal="Sample_BALBc_normal_C3_IGO_07035_B_4"
normal_bam=${data_path}/${normal}/${normal}.aligned.duplicates_marked.recalibrated.bam  # recalibrated BAM file

tumor=TUMOR_SAMPLE
tumor_bam=${data_path}/${tumor}/${tumor}.aligned.duplicates_marked.recalibrated.bam


if [[ ! -e ${normal_bam/.bam/.mpileup} ]];then
samtools mpileup -f ${bash_ref_fasta}/GRCm38.fa ${normal_bam} > ${normal_bam/.bam/.mpileup}
fi

if [[ ! -e ${tumor_bam/.bam/.mpileup} ]];then
samtools mpileup -f ${bash_ref_fasta}/GRCm38.fa  ${tumor_bam} >${tumor_bam/.bam/.mpileup}
fi

java -jar /home/luol2/lingqi_workspace/VarScan.v2.3.9.jar somatic ${normal_bam/.bam/.mpileup} ${tumor_bam/.bam/.mpileup} ${data_path}/${tumor}/${tumor} --min-var-freq 0.01



