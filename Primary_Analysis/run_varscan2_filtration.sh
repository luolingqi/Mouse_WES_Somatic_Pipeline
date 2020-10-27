#BSUB -o Myjob.%J.log
#BSUB -e Myjob.%J.err
#BSUB -n 8 -R "rusage[mem=8]"
#BSUB -W 48:00

module load java/1.8.0_31
module load singularity/3.1.1

data_path="/home/luol2/lingqi_workspace/Projects/Ben_Projects/WES_mouse"

tumor=TUMOR_SAMPLE

# Run processSomatic and isolate the high-confidence somatic calls
java -jar /home/luol2/lingqi_workspace/VarScan.v2.3.9.jar processSomatic \
${data_path}/${tumor}/${tumor}.snp.vcf \
--min-tumor-freq 0.01 \
--max-normal-freq 0.05 \
--p-value 0.07

java -jar /home/luol2/lingqi_workspace/VarScan.v2.3.9.jar processSomatic \
${data_path}/${tumor}/${tumor}.indel.vcf \
--min-tumor-freq 0.01 \
--max-normal-freq 0.05 \
--p-value 0.07


# Generate bam-readcount results for those HC calls only, using the TUMOR BAM file
   # Prepare site list file
   grep -v "^#" ${data_path}/${tumor}/${tumor}.snp.Somatic.hc.vcf | awk '{OFS="\t"}{print $1,$2,$2}' > ${data_path}/${tumor}/${tumor}.snp.Somatic.hc.sites
   grep -v "^#" ${data_path}/${tumor}/${tumor}.indel.Somatic.hc.vcf | awk '{OFS="\t"}{print $1,$2-20,$2+20}' > ${data_path}/${tumor}/${tumor}.indel.Somatic.hc.sites
   # Generate bam-readcount
   singularity exec --bind /home/luol2/lingqi_workspace/References/GRCm38/:/ref --bind ${data_path}/${tumor}/:/tmp  /home/luol2/lingqi_workspace/bam-readcount.sif bam-readcount -w 10 -f /ref/GRCm38.fa --site-list /tmp/${tumor}.snp.Somatic.hc.sites /tmp/${tumor}.aligned.duplicates_marked.recalibrated.bam >${data_path}/${tumor}/${tumor}.aligned.duplicates_marked.recalibrated.snp.Somatic.hc.bam-readcount.txt
   singularity exec --bind /home/luol2/lingqi_workspace/References/GRCm38/:/ref --bind ${data_path}/${tumor}/:/tmp  /home/luol2/lingqi_workspace/bam-readcount.sif bam-readcount -w 10 -f /ref/GRCm38.fa --site-list /tmp/${tumor}.indel.Somatic.hc.sites /tmp/${tumor}.aligned.duplicates_marked.recalibrated.bam >${data_path}/${tumor}/${tumor}.aligned.duplicates_marked.recalibrated.indel.Somatic.hc.bam-readcount.txt


# Run fpfilter command
java -jar /home/luol2/lingqi_workspace/VarScan.v2.3.9.jar fpfilter ${data_path}/${tumor}/${tumor}.snp.Somatic.hc.vcf \
${data_path}/${tumor}/${tumor}.aligned.duplicates_marked.recalibrated.snp.Somatic.hc.bam-readcount.txt \
--min-var-count 5 \
--min-ref-basequal 20 \
--min-var-basequal 20 \
--min-var-freq 0.01 \
--output-file ${data_path}/${tumor}/${tumor}.snp.Somatic.hc.fpfilter.vcf

java -jar /home/luol2/lingqi_workspace/VarScan.v2.3.9.jar fpfilter ${data_path}/${tumor}/${tumor}.indel.Somatic.hc.vcf \
${data_path}/${tumor}/${tumor}.aligned.duplicates_marked.recalibrated.indel.Somatic.hc.bam-readcount.txt \
--min-var-count 5 \
--min-ref-basequal 20 \
--min-var-basequal 20 \
--min-var-freq 0.01 \
--output-file ${data_path}/${tumor}/${tumor}.indel.Somatic.hc.fpfilter.vcf

# combine snp and indel variants and indexing
java -jar /home/luol2/lingqi_workspace/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T CombineVariants -R /home/luol2/lingqi_workspace/References/GRCm38/GRCm38.fa --variant:foo ${data_path}/${tumor}/${tumor}.snp.Somatic.hc.fpfilter.vcf --variant:bar ${data_path}/${tumor}/${tumor}.indel.Somatic.hc.fpfilter.vcf -o ${data_path}/${tumor}/${tumor}.snp_indel.Somatic.hc.fpfilter.vcf -genotypeMergeOptions PRIORITIZE -priority foo,bar
