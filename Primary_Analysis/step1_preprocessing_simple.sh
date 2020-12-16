######################################################################################################################
# This script automates all the preprocessing steps described in GATK Best Practice
############################################################################################################################

###############################
#parameters
###############################
data_path=$1;export data_path
project=$2
subject=$3
sample=$4
#readGroup=$5;export readGroup # e.g. H2JFYBBXY.5
jobId=""
jobName=""
message=""

# function parameters: $1: step name; $2:
function runSteps () {
    #export data_path project subject sample readGroup
    #perl -pe 's#TEST_SAMPLE#$ENV{sample}#g;s#PROJECT#$ENV{project}#g;s#SUBJECT#$ENV{subject}#g;s#DATA_PATH#$ENV{data_path}#g;s#TRIM##g' run_${1}.sh |bsub -J ${1}_${2}_${3}_${4}
    sed "s#TEST_SAMPLE#${4}#g;s#PROJECT#${2}#g;s#SUBJECT#${3}#g;s#DATA_PATH#${5}#g;s#TRIM##g" run_${1}.sh |bsub -J ${1}_${2}_${3}_${4}
    #unset data_path project subject sample readGroup
    jobId=$(bjobs -J ${1}_${2}_${3}_${4} | awk '{print $1}' | grep -v JOBID)
    jobName=${1}_${2}_${3}_${4}
    message="${1} going......"
    checkJobSuccess
    sleep 30
}


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

#1. Fastqc on the raw fastq files
runSteps fastqc $project $subject $sample $data_path

#2. Trimgalore to remove adapter and low quality reads
runSteps trim_galore $project $subject $sample $data_path
mv ${data_path}/${project}/${subject}/${sample}/trim_galore_output/*trim.fastq.gz ${data_path}/${project}/${subject}/${sample}
sleep 10

#3. Fastqc on the trimmed fastq files
#runSteps fastqc_trim

#4. Convert trimmed fastq to uBAM
runSteps fastq_to_uBAM $project $subject $sample $data_path

#5. BWA MEM Alignment
runSteps bwa_mem $project $subject $sample $data_path

#6. Collect Quality Metrics
#6.1 Qualimap to gather both exome and rnaseq qc metrics
runSteps qualimap $project $subject $sample $data_path

#6.2 HsMetrics to collect hybrid selection metrics for WES data
runSteps CollectHsMetrics $project $subject $sample $data_path

#6.3 MultipleMetrics
runSteps CollectQualityMetrics $project $subject $sample $data_path

#7. MarkDuplicate
runSteps markduplicate $project $subject $sample $data_path

#   INDELREALIGNER
runSteps indelrealigner $project $subject $sample $data_path
#8. Estimate Base Quality
runSteps BaseRecalibrator $project $subject $sample $data_path

#9. Apply BQSR
runSteps ApplyBQSR $project $subject $sample $data_path

# Declare mission completed
echo "Mission Accomplished!"


