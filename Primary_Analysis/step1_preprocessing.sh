######################################################################################################################
# This script automates all the preprocessing steps described in GATK Best Practice
############################################################################################################################

###############################
#parameters
###############################
data_path=$1;export data_path
project=$2;export project
subject=$3;export subject
sample=$4;export sample
#readGroup=$5;export readGroup # e.g. H2JFYBBXY.5
jobId=""
jobName=""
message=""
function checkJobSuccess {
    # track if the job stops running    
    jobPrev=""
    jobId=$(bjobs -J ${jobName} | awk '{print $1}' | grep -v JOBID)
    until [ "$jobId" == "" ];do
	sleep 10
	echo ${message}
	err=$(bjobs -J ${jobName} 2>&1)
	if [[ $err == *"not found"* ]]; then
	    jobPrev=${jobId}
	    jobId=""
	fi
	echo ${jobId}
    done

    # track if the job succeeds running. Exit immediately with a failed LSF job
    sleep 20
    success=$(grep "Successfully completed" Myjob.${jobPrev}.log)
    #success=$(bhist -l ${jobPrev} | grep  "Done successfully")
    echo ${success}
    if [ -z "${success}" ]; then
	exit 1;
    fi
}
#0. Extract/build readgroup from the sample fastq.gz file
readGroup=$(gzip -cd ${data_path}/${project}/${subject}/${sample}/*R1_001.fastq.gz | head -1 | cut -d: -f3,4 | perl -pe 's/:/./g')
export readGroup

#1. Fastqc on the raw fastq files
perl -pe 's#TEST_SAMPLE#$ENV{sample}#g;s#PROJECT#$ENV{project}#g;s#SUBJECT#$ENV{subject}#g;s#DATA_PATH#$ENV{data_path}#g;s#TRIM##g' run_fastqc.sh |bsub -J fastqc_${project}_${subject}_${sample}
jobId=$(bjobs -J fastqc_${project}_${subject}_${sample} | awk '{print $1}' | grep -v JOBID)
jobName=fastqc_${project}_${subject}_${sample}
message="fastqc going......"
checkJobSuccess
sleep 30

#2. Trimgalore to remove adapter and low quality reads
perl -pe 's#TEST_SAMPLE#$ENV{sample}#g;s#PROJECT#$ENV{project}#g;s#SUBJECT#$ENV{subject}#g;s#DATA_PATH#$ENV{data_path}#g' run_trim_galore.sh | bsub -J trimGalore_${project}_${subject}_${sample}
jobId=$(bjobs -J trimGalore_${project}_${subject}_${sample} | awk '{print $1}' | grep -v JOBID)
jobName=trimGalore_${project}_${subject}_${sample}
message="trimGalore going......"
checkJobSuccess
sleep 30
mv ${data_path}/${project}/${subject}/${sample}/trim_galore_output/*trim.fastq.gz ${data_path}/${project}/${subject}/${sample}
sleep 10

#3. Fastqc on the trimmed fastq files
perl -pe 's#TEST_SAMPLE#$ENV{sample}#g;s#PROJECT#$ENV{project}#g;s#SUBJECT#$ENV{subject}#g;s#DATA_PATH#$ENV{data_path}#g;s#TRIM#_trim#g' run_fastqc.sh |bsub -J fastqc_${project}_${subject}_${sample}_trim
jobId=$(bjobs -J fastqc_${project}_${subject}_${sample}_trim | awk '{print $1}' | grep -v JOBID)
jobName=fastqc_${project}_${subject}_${sample}_trim
message="fastqc after trimming going......"
checkJobSuccess
sleep 30

#4. Convert trimmed fastq to uBAM
perl -pe 's/DATA_PATH/$ENV{data_path}/g;s/PROJECT/$ENV{project}/g;s/SUBJECT/$ENV{subject}/g;s/TEST_SAMPLE/$ENV{sample}/g;s/READGROUP/$ENV{readGroup}/g' run_fastq_to_uBAM.sh | bsub -J uBAM_${project}_${subject}_${sample}
jobId=$(bjobs -J uBAM_${project}_${subject}_${sample} | awk '{print $1}' | grep -v JOBID)
jobName=uBAM_${project}_${subject}_${sample}
message="uBAM Conversion going......"
checkJobSuccess
sleep 30


#5. BWA MEM Alignment
perl -pe 's/DATA_PATH/$ENV{data_path}/g;s/PROJECT/$ENV{project}/g;s/SUBJECT/$ENV{subject}/g;s/TEST_SAMPLE/$ENV{sample}/g' run_bwa_mem.sh | bsub -J bwa_${project}_${subject}_${sample}
jobId=$(bjobs -J bwa_${project}_${subject}_${sample} | awk '{print $1}' | grep -v JOBID)
jobName=bwa_${project}_${subject}_${sample}
message="BWA Align going......."
checkJobSuccess
sleep 30


#6. Collect Quality Metrics
#6.1 Qualimap to gather both exome and rnaseq qc metrics
perl -pe 's/DATA_PATH/$ENV{data_path}/g;s/PROJECT/$ENV{project}/g;s/SUBJECT/$ENV{subject}/g;s/TEST_SAMPLE/$ENV{sample}/g' run_qualimap.sh | bsub -J qualimap_${project}_${subject}_${sample}
jobId=$(bjobs -J qualimap_${project}_${subject}_${sample} | awk '{print $1}' | grep -v JOBID)
jobName=qualimap_${project}_${subject}_${sample}
message="Qualimap going......."
checkJobSuccess
sleep 30
#}
#fi

#6.2 HsMetrics to collect hybrid selection metrics for WES data
perl -pe 's/DATA_PATH/$ENV{data_path}/g;s/PROJECT/$ENV{project}/g;s/SUBJECT/$ENV{subject}/g;s/TEST_SAMPLE/$ENV{sample}/g' run_CollectHsMetrics.sh | bsub -J Hs_${project}_${subject}_${sample}
jobId=$(bjobs -J Hs_${project}_${subject}_${sample} | awk '{print $1}' | grep -v JOBID)
jobName=Hs_${project}_${subject}_${sample}
message="HsMetrics going......"
checkJobSuccess
sleep 30

#6.3 MultipleMetrics
perl -pe 's/DATA_PATH/$ENV{data_path}/g;s/PROJECT/$ENV{project}/g;s/SUBJECT/$ENV{subject}/g;s/TEST_SAMPLE/$ENV{sample}/g' run_CollectQualityMetrics.sh | bsub -J Quality_${project}_${subject}_${sample}
jobId=$(bjobs -J Quality_${project}_${subject}_${sample} | awk '{print $1}' | grep -v JOBID)
jobName=Quality_${project}_${subject}_${sample}
message="QualityMetrics going......"
checkJobSuccess
sleep 30

#7. MarkDuplicate
perl -pe 's/DATA_PATH/$ENV{data_path}/g;s/PROJECT/$ENV{project}/g;s/SUBJECT/$ENV{subject}/g;s/TEST_SAMPLE/$ENV{sample}/g' run_markduplicate.sh | bsub -J markdup_${project}_${subject}_${sample}
jobId=$(bjobs -J markdup_${project}_${subject}_${sample} | awk '{print $1}' | grep -v JOBID)
jobName=markdup_${project}_${subject}_${sample}
message="Mark Duplicate going......."
checkJobSuccess

#8. Estimate Base Quality
perl -pe 's/DATA_PATH/$ENV{data_path}/g;s/PROJECT/$ENV{project}/g;s/SUBJECT/$ENV{subject}/g;s/TEST_SAMPLE/$ENV{sample}/g' run_BaseRecalibrator.sh | bsub -J BQSR_${project}_${subject}_${sample}
jobId=$(bjobs -J BQSR_${project}_${subject}_${sample} | awk '{print $1}' | grep -v JOBID)
jobName=BQSR_${project}_${subject}_${sample}
message="Estimate Base Quality going......."
checkJobSuccess

#9. Apply BQSR
perl -pe 's/DATA_PATH/$ENV{data_path}/g;s/PROJECT/$ENV{project}/g;s/SUBJECT/$ENV{subject}/g;s/TEST_SAMPLE/$ENV{sample}/g' run_ApplyBQSR.sh | bsub -J ApplyBQSR_${project}_${subject}_${sample}
jobId=$(bjobs -J ApplyBQSR_${project}_${subject}_${sample} | awk '{print $1}' | grep -v JOBID)
jobName=ApplyBQSR_${project}_${subject}_${sample}
message="Apply BQSR going......."
checkJobSuccess

# Declare mission completed
echo "Mission Accomplished!"


