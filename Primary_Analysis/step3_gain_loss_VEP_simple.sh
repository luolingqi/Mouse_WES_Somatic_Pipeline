######################################################################################################################
# This script automates all the preprocessing steps described in GATK Best Practice
############################################################################################################################

###############################
#parameters
###############################
data_path=$1
project=$2
subject=$3
p_sample=$4
t_sample=$5

#readGroup=$5;export readGroup # e.g. H2JFYBBXY.5
jobId=""
jobName=""
message=""
function runSteps () {
    sed "s#PARENTAL_SAMPLE#${4}#g;s#TREATED_SAMPLE#${5}#g;s#PROJECT#${2}#g;s#SUBJECT#${3}#g;s#DATA_PATH#${6}#g;" run_${1}.sh |bsub -J ${1}_${2}_${3}_${4}_vs_${5}
    #unset data_path project subject sample readGroup
    jobId=$(bjobs -J ${1}_${2}_${3}_${4}_vs_${5} | awk '{print $1}' | grep -v JOBID)
    jobName=${1}_${2}_${3}_${4}_vs_${5}
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


#1 acquire variant gain/loss against parental 
runSteps variants_gain_loss_treatment.vs.Parental.AF.0.05 $project $subject $p_sample $t_sample $data_path


# Declare mission completed
echo "Mission Accomplished!"




