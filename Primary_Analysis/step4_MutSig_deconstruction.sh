######################################################################################################################
# This script automates all the preprocessing steps described in GATK Best Practice
############################################################################################################################

###############################
#parameters
###############################
data_path=$1
project=$2
manifest="$3"

#readGroup=$5;export readGroup # e.g. H2JFYBBXY.5
jobId=""
jobName=""
message=""
function runSteps () {
    sed "s#PATH#${2}/${3}#g;s#MANIFEST#${4}#g;" run_${1}.sh |bsub -J ${1}_${2}_${3}_${4}
    #unset data_path project subject sample readGroup
    jobId=$(bjobs -J ${1}_${2}_${3}_${4} | awk '{print $1}' | grep -v JOBID)
    jobName=${1}_${2}_${3}_${4}
    message="${1} going...... ${2} ${3} ${4} "
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

#1 run deconstructSig for all the samples and variant types listed in the manifest file
runSteps MutSig_deconstruction ${data_path} ${project} ${manifest}

# Declare mission completed
echo "Mission Accomplished!"




