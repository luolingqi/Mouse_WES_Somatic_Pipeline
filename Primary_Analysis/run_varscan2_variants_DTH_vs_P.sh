#BSUB -o Myjob.%J.log
#BSUB -e Myjob.%J.err
#BSUB -n 1 -R "rusage[mem=8]"
#BSUB -W 48:00

module load java/1.8.0_31
ref_fasta="/home/luol2/lingqi_workspace/References/GRCm38/GRCm38.fa"

data_path="/home/luol2/lingqi_workspace/Projects/Ben_Projects/WES_mouse"
sample_ctl="Sample_CT26_P_3m_IGO_10212_1"
sample_treatment="Sample_CT26_DTH_8W_IGO_10212_2"

vcf_ctl=${data_path}/${sample_ctl}/${sample_ctl}.snp_indel.Somatic.hc.fpfilter.vcf
vcf_treatment=${data_path}/${sample_treatment}/${sample_treatment}.snp_indel.Somatic.hc.fpfilter.vcf

command_mem="100"

/home/luol2/lingqi_workspace/gatk-4.0.1.2/gatk  --java-options "-Xmx${command_mem}m" SelectVariants -R ${ref_fasta} -V ${vcf_treatment} --discordance ${vcf_ctl} -O ${data_path}/data_analysis/${sample_ctl}_vs_${sample_treatment}.snp_indel.Somatic.hc.fpfilter.vcf

