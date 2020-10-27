#BSUB -o Myjob.%J.log
#BSUB -e Myjob.%J.err
#BSUB -n 8 -R "rusage[mem=8]"
#BSUB -W 48:00

### Refer to  ~/lingqi_workspace/Projects/Mitesh/WGS/intel-gatk4-somatic-with-preprocessing/mutect2_nodocker.wdl

module load java/1.8.0_31
ref_fasta="/data/ldiaz/luol2/References/GRCm38/GRCm38.fa"

#data_path="/data/ldiaz/luol2/Projects/Jones_Lee_Exercise_Project"
data_path=DATA_PATH
project=PROJECT
subject=SUBJECT
sample=TEST_SAMPLE
normal_sample=NORMAL_SAMPLE
#normal_sample=${data_path}/${project}/Sample_BALBc_Norm_IGO_10212_E_1/Sample_BALBc_Norm_IGO_10212_E_1

tumor_bam=${data_path}/${project}/${subject}/${sample}/${sample}.aligned.duplicates_marked.recalibrated.bam
normal_bam=${data_path}/${project}/${subject}/${normal_sample}/${normal_sample}.aligned.duplicates_marked.recalibrated.bam
tumor_name=${data_path}/${project}/${subject}/${sample}/tumor_name.txt
normal_name=${data_path}/${project}/${subject}/${normal_sample}/normal_name.txt

unfiltered_output_vcf=${tumor_bam/.bam/.unfiltered.vcf}
output_vcf=${tumor_bam/.bam/.filtered.vcf}
command_mem="8000"

# index the bam
if [ ! -e ${tumor_bam/.bam/.bai} ];then
module load samtools/1.7
samtools index ${tumor_bam}
fi

# Run Mutect2

        # We need to create these files regardless, even if they stay empty
        touch ${data_path}/${project}/${subject}/${sample}/bamout.bam
        echo "" > ${normal_name}

        /data/ldiaz/luol2/gatk-4.1.2.0/gatk --java-options "-Xmx${command_mem}m" GetSampleName -R ${ref_fasta} -I ${tumor_bam} -O ${tumor_name} -encode
        tumor_command_line="-I ${tumor_bam} -tumor `cat ${tumor_name}`"

	if [[ -f "${normal_bam}" ]]; then
            /data/ldiaz/luol2/gatk-4.1.2.0/gatk --java-options "-Xmx${command_mem}m" GetSampleName -R ${ref_fasta} -I ${normal_bam} -O ${normal_name} -encode
            normal_command_line="-I ${normal_bam} -normal `cat ${normal_name}`"
        fi

	# copy normal name to sample folder also
	copy ${normal_name} ${data_path}/${project}/${subject}/${sample}

        /data/ldiaz/luol2/gatk-4.1.2.0/gatk --java-options "-Xmx${command_mem}m" Mutect2 \
            -R ${ref_fasta} \
            $tumor_command_line \
            $normal_command_line \
            -L ${data_path}/data_analysis/mouse_all_exon_mm10.interval_list \
            -O "${unfiltered_output_vcf}" \
            -pairHMM AVX_LOGLESS_CACHING \
            --native-pair-hmm-threads 1 \
            --smith-waterman AVX_ENABLED 
            
           


# Run variant Filtration
        /data/ldiaz/luol2/gatk-4.1.2.0/gatk --java-options "-Xmx${command_mem}m" FilterMutectCalls -V ${unfiltered_output_vcf} \
                                    -R ${ref_fasta} \
                                    -O ${output_vcf}




