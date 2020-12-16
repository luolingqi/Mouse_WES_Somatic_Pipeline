######################################################################################################################
# This script automates all the preprocessing steps described in GATK Best Practice
############################################################################################################################

###############################
#parameters
###############################
data_path=$1
project=$2
subject=$3
sample=$4
normal_sample=$5
#readGroup=$5;export readGroup # e.g. H2JFYBBXY.5
jobId=""
jobName=""
message=""
function runSteps () {
    sed "s#TEST_SAMPLE#${4}#g;s#PROJECT#${2}#g;s#SUBJECT#${3}#g;s#DATA_PATH#${5}#g;s#NORMAL_SAMPLE#${6}#g;s#TRIM##g" run_${1}.sh |bsub -J ${1}_${2}_${3}_${4}
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


#0. Copy mouse_all_exon_mm10.interval_list over to current folder (data_analysis)
cp /home/luol2/lingqi_workspace/Mouse_WES_Pipeline/Primary/mouse_all_exon_mm10.interval_list .

#1. Mutect2 Calling and Filtration
runSteps mutect2_and_Filter $project $subject $sample $data_path $normal_sample


#2. Remove the Germline dbsnp variants for BALB/cJ strain
runSteps remove_C57BL_6NJ_germline_snp_indels $project $subject $sample $data_path $normal_sample


#3. VEP Annotation, extra manual filtration and variant classification by type
runSteps VEP_annotation_C57BL_6NJ $project $subject $sample $data_path $normal_sample


# Declare mission completed
echo "Mission Accomplished!"




