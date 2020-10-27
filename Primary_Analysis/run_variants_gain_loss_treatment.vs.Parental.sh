#BSUB -o Myjob.%J.log
#BSUB -e Myjob.%J.err
#BSUB -n 1 -R "rusage[mem=8]"
#BSUB -W 48:00

module load java/1.8.0_31
ref_fasta="/home/luol2/lingqi_workspace/References/GRCm38/GRCm38.fa"


data_path="/home/luol2/lingqi_workspace/Projects/Ben_Projects/WES_mouse_01062020"
project=PROJECT
subject=SUBJECT
sample_ctl="PARENTAL_SAMPLE"
sample_treatment="TREATED_SAMPLE"

vcf_ctl=${data_path}/${project}/${subject}/${sample_ctl}/${sample_ctl}.aligned.duplicates_marked.recalibrated.filtered.w.BALCnormal.rm_MGP_snps_indels.singleSample.PASSOnly.vcf
vcf_treatment=${data_path}/${project}/${subject}/${sample_treatment}/${sample_treatment}.aligned.duplicates_marked.recalibrated.filtered.w.BALCnormal.rm_MGP_snps_indels.singleSample.PASSOnly.vcf


#######################################################################
# Variant gain against Parental
####################################################################### 

# grab difference as treatment variant gain from Parental
command_mem="100"
/home/luol2/lingqi_workspace/gatk-4.0.1.2/gatk  --java-options "-Xmx${command_mem}m" SelectVariants -R ${ref_fasta} -V ${vcf_treatment} --discordance ${vcf_ctl} -O ${data_path}/${project}/${subject}/${sample_ctl}_vs_${sample_treatment}.w.BALCnormal.rm_MGP_snps_indels.singleSample.PASSOnly.vcf

module load singularity/3.1.1
VEPInputVcf=${sample_ctl}_vs_${sample_treatment}.w.BALCnormal.rm_MGP_snps_indels.singleSample.PASSOnly.vcf
VEPOutputVcf=${sample_ctl}_vs_${sample_treatment}.w.BALCnormal.rm_MGP_snps_indels.singleSample.PASSOnly.VEP.ann.vcf
VEPOutputStat=${sample_ctl}_vs_${sample_treatment}.w.BALCnormal.rm_MGP_snps_indels.singleSample.PASSOnly.VEP.ann.summary.html
#Perform VEP using local cache (Mouse)
singularity run --bind ${data_path}/${project}/${subject}:/data --bind /home/luol2/lingqi_workspace/vep_data:/vep_data \
/home/luol2/lingqi_workspace/ensembl-vep.sif \
vep -i /data/${VEPInputVcf} \
    -o /data/${VEPOutputVcf} \
    --vcf \
    --force_overwrite \
    --cache \
    --dir_cache /vep_data \
    --species mouse \
    --stats_file /data/${VEPOutputStat}


#Perform VEP summarization after manual filtration
bcftools view -i 'FORMAT/AD[0:1] > 4 & FORMAT/AF[0:0] > 0.01 & INFO/MBQ[1] > 20 & INFO/MMQ[1] > 50 & INFO/MPOS > 5' ${data_path}/${project}/${subject}/${VEPInputVcf} > ${data_path}/${project}/${subject}/${VEPInputVcf/.vcf/.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.vcf} 
singularity exec --bind ${data_path}/${project}/${subject}:/data --bind /home/luol2/lingqi_workspace/vep_data:/vep_data \
/data/ldiaz/luol2/ensembl-vep.sif \
vep -i /data/${VEPInputVcf/.vcf/.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.vcf} \
    -o /data/${VEPInputVcf/.vcf/.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.VEP.ann.vcf} \
    --vcf \
    --force_overwrite \
    --cache \
    --dir_cache /vep_data \
    --species mouse \
    --stats_file /data/${VEPInputVcf/.vcf/.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.VEP.ann.summary.html}

# VEP output in tab format
singularity exec --bind ${data_path}/${project}/${subject}/${sample}:/data --bind /home/luol2/lingqi_workspace/vep_data:/vep_data \
/data/ldiaz/luol2/ensembl-vep.sif \
vep -i /data/${VEPInputVcf/.vcf/.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.vcf} \
    -o /data/${VEPInputVcf/.vcf/.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.VEP.ann.txt} \
    --tab \
    --fields "Uploaded_variation,Location,Allele,Gene,CANONICAL,SYMBOL,SYMBOL_SOURCE,Feature,Feature_type,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,IMPACT,DISTANCE,STRAND" \
    --force_overwrite \
    --cache \
    --dir_cache /vep_data \
    --species mus_musculus \
    --symbol


# Perform VEP filtrations:
# Generate a grand VCF file with the manual filtration
bcftools view -i 'FORMAT/AD[0:1] > 4 & FORMAT/AF[0:0] > 0.01 & INFO/MBQ[1] > 20 & INFO/MMQ[1] > 50 & INFO/MPOS > 5' ${data_path}/${project}/${subject}/${VEPOutputVcf} > ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.vcf}
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT{0}\t%DP\t[%GT\t%AD{0}\t%AD{1}\t%AF\t]\t%MBQ[1]\t%MMQ[1]\t%MPOS\t\n" ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.vcf} > ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.query}

singularity exec --bind ${data_path}/${project}/${subject}:/data \
/home/luol2/lingqi_workspace/ensembl-vep.sif \
filter_vep -i /data/${VEPOutputVcf} --filter "Consequence in missense_variant,frameshift_variant,inframe_deletion,coding_sequence_variant,stop_gained,inframe_insertion,splice_region_variant" | \
bcftools view -i 'FORMAT/AD[0:1] > 4 & FORMAT/AF[0:0] > 0.01 & INFO/MBQ[1] > 20 & INFO/MMQ[1] > 50 & INFO/MPOS > 5' > ${data_path}//${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.non_synonymous_variant_only.vcf}
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT{0}\t%DP\t[%GT\t%AD{0}\t%AD{1}\t%AF\t]\t%MBQ[1]\t%MMQ[1]\t%MPOS\t\n" ${data_path}//${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.non_synonymous_variant_only.vcf} > ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.non_synonymous_variant_only.query}

singularity exec --bind ${data_path}/${project}/${subject}:/data \
/home/luol2/lingqi_workspace/ensembl-vep.sif \
filter_vep -i /data/${VEPOutputVcf} --filter "Consequence is missense_variant" | \
bcftools view -i 'FORMAT/AD[0:1] > 4 & FORMAT/AF[0:0] > 0.01 & INFO/MBQ[1] > 20 & INFO/MMQ[1] > 50 & INFO/MPOS > 5' > ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.missense_variant_only.vcf}
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT{0}\t%DP\t[%GT\t%AD{0}\t%AD{1}\t%AF\t]\t%MBQ[1]\t%MMQ[1]\t%MPOS\t\n" ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.missense_variant_only.vcf} > ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.missense_variant_only.query}

singularity exec --bind ${data_path}/${project}/${subject}:/data \
/home/luol2/lingqi_workspace/ensembl-vep.sif \
filter_vep -i /data/${VEPOutputVcf} --filter "Consequence is frameshift_variant" | \
bcftools view -i 'FORMAT/AD[0:1] > 4 & FORMAT/AF[0:0] > 0.01 & INFO/MBQ[1] > 20 & INFO/MMQ[1] > 50 & INFO/MPOS > 5' > ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.frameshift_variant_only.vcf}
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT{0}\t%DP\t[%GT\t%AD{0}\t%AD{1}\t%AF\t]\t%MBQ[1]\t%MMQ[1]\t%MPOS\t\n" ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.frameshift_variant_only.vcf} > ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.frameshift_variant_only.query}

singularity exec --bind ${data_path}/${project}/${subject}:/data \
/home/luol2/lingqi_workspace/ensembl-vep.sif \
filter_vep -i /data/${VEPOutputVcf} --filter "Consequence is synonymous_variant" | \
bcftools view -i 'FORMAT/AD[0:1] > 4 & FORMAT/AF[0:0] > 0.01 & INFO/MBQ[1] > 20 & INFO/MMQ[1] > 50 & INFO/MPOS > 5' > ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.synonymous_variant_only.vcf}
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT{0}\t%DP\t[%GT\t%AD{0}\t%AD{1}\t%AF\t]\t%MBQ[1]\t%MMQ[1]\t%MPOS\t\n" ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.synonymous_variant_only.vcf} > ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.synonymous_variant_only.query}


#######################################################################
# Variant loss against Parental
####################################################################### 

# grab difference as treatment variant gain from Parental
command_mem="100"
/home/luol2/lingqi_workspace/gatk-4.0.1.2/gatk  --java-options "-Xmx${command_mem}m" SelectVariants -R ${ref_fasta} -V ${vcf_ctl} --discordance ${vcf_treatment} -O ${data_path}/${project}/${subject}/${sample_treatment}_vs_${sample_ctl}.w.BALCnormal.rm_MGP_snps_indels.singleSample.PASSOnly.vcf

module load singularity/3.1.1
VEPInputVcf=${sample_treatment}_vs_${sample_ctl}.w.BALCnormal.rm_MGP_snps_indels.singleSample.PASSOnly.vcf
VEPOutputVcf=${sample_treatment}_vs_${sample_ctl}.w.BALCnormal.rm_MGP_snps_indels.singleSample.PASSOnly.VEP.ann.vcf
VEPOutputStat=${sample_treatment}_vs_${sample_ctl}.w.BALCnormal.rm_MGP_snps_indels.singleSample.PASSOnly.VEP.ann.summary.html
#Perform VEP using local cache (Mouse)
singularity run --bind ${data_path}/${project}/${subject}:/data --bind /home/luol2/lingqi_workspace/vep_data:/vep_data \
/home/luol2/lingqi_workspace/ensembl-vep.sif \
vep -i /data/${VEPInputVcf} \
    -o /data/${VEPOutputVcf} \
    --vcf \
    --force_overwrite \
    --cache \
    --dir_cache /vep_data \
    --species mouse \
    --stats_file /data/${VEPOutputStat}

#Perform VEP summarization after manual filtration
bcftools view -i 'FORMAT/AD[0:1] > 4 & FORMAT/AF[0:0] > 0.01 & INFO/MBQ[1] > 20 & INFO/MMQ[1] > 50 & INFO/MPOS > 5' ${data_path}/${project}/${subject}/${VEPInputVcf} > ${data_path}/${project}/${subject}/${VEPInputVcf/.vcf/.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.vcf} 
singularity exec --bind ${data_path}/${project}/${subject}:/data --bind /home/luol2/lingqi_workspace/vep_data:/vep_data \
/data/ldiaz/luol2/ensembl-vep.sif \
vep -i /data/${VEPInputVcf/.vcf/.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.vcf} \
    -o /data/${VEPInputVcf/.vcf/.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.VEP.ann.vcf} \
    --vcf \
    --force_overwrite \
    --cache \
    --dir_cache /vep_data \
    --species mouse \
    --stats_file /data/${VEPInputVcf/.vcf/.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.VEP.ann.summary.html}

# VEP output in tab format
singularity exec --bind ${data_path}/${project}/${subject}/${sample}:/data --bind /home/luol2/lingqi_workspace/vep_data:/vep_data \
/data/ldiaz/luol2/ensembl-vep.sif \
vep -i /data/${VEPInputVcf/.vcf/.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.vcf} \
    -o /data/${VEPInputVcf/.vcf/.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.VEP.ann.txt} \
    --tab \
    --fields "Uploaded_variation,Location,Allele,Gene,CANONICAL,SYMBOL,SYMBOL_SOURCE,Feature,Feature_type,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,IMPACT,DISTANCE,STRAND" \
    --force_overwrite \
    --cache \
    --dir_cache /vep_data \
    --species mus_musculus \
    --symbol


# Perform VEP filtrations:
singularity exec --bind ${data_path}/${project}/${subject}:/data \
/home/luol2/lingqi_workspace/ensembl-vep.sif \
filter_vep -i /data/${VEPOutputVcf} --filter "Consequence in missense_variant,frameshift_variant,inframe_deletion,coding_sequence_variant,stop_gained,inframe_insertion,splice_region_variant" | \
bcftools view -i 'FORMAT/AD[0:1] > 4 & FORMAT/AF[0:0] > 0.01 & INFO/MBQ[1] > 20 & INFO/MMQ[1] > 50 & INFO/MPOS > 5' > ${data_path}//${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.non_synonymous_variant_only.vcf}
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT{0}\t%DP\t[%GT\t%AD{0}\t%AD{1}\t%AF\t]\t%MBQ[1]\t%MMQ[1]\t%MPOS\t\n" ${data_path}//${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.non_synonymous_variant_only.vcf} > ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.non_synonymous_variant_only.query}

singularity exec --bind ${data_path}/${project}/${subject}:/data \
/home/luol2/lingqi_workspace/ensembl-vep.sif \
filter_vep -i /data/${VEPOutputVcf} --filter "Consequence is missense_variant" | \
bcftools view -i 'FORMAT/AD[0:1] > 4 & FORMAT/AF[0:0] > 0.01 & INFO/MBQ[1] > 20 & INFO/MMQ[1] > 50 & INFO/MPOS > 5' > ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.missense_variant_only.vcf}
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT{0}\t%DP\t[%GT\t%AD{0}\t%AD{1}\t%AF\t]\t%MBQ[1]\t%MMQ[1]\t%MPOS\t\n" ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.missense_variant_only.vcf} > ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.missense_variant_only.query}

singularity exec --bind ${data_path}/${project}/${subject}:/data \
/home/luol2/lingqi_workspace/ensembl-vep.sif \
filter_vep -i /data/${VEPOutputVcf} --filter "Consequence is frameshift_variant" | \
bcftools view -i 'FORMAT/AD[0:1] > 4 & FORMAT/AF[0:0] > 0.01 & INFO/MBQ[1] > 20 & INFO/MMQ[1] > 50 & INFO/MPOS > 5' > ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.frameshift_variant_only.vcf}
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT{0}\t%DP\t[%GT\t%AD{0}\t%AD{1}\t%AF\t]\t%MBQ[1]\t%MMQ[1]\t%MPOS\t\n" ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.frameshift_variant_only.vcf} > ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.frameshift_variant_only.query}

singularity exec --bind ${data_path}/${project}/${subject}:/data \
/home/luol2/lingqi_workspace/ensembl-vep.sif \
filter_vep -i /data/${VEPOutputVcf} --filter "Consequence is synonymous_variant" | \
bcftools view -i 'FORMAT/AD[0:1] > 4 & FORMAT/AF[0:0] > 0.01 & INFO/MBQ[1] > 20 & INFO/MMQ[1] > 50 & INFO/MPOS > 5' > ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.synonymous_variant_only.vcf}
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT{0}\t%DP\t[%GT\t%AD{0}\t%AD{1}\t%AF\t]\t%MBQ[1]\t%MMQ[1]\t%MPOS\t\n" ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.synonymous_variant_only.vcf} > ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.01_MBQ20_MMQ50_MPOS5.synonymous_variant_only.query}
