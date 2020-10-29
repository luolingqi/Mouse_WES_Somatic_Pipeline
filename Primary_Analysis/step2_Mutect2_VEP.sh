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
normal_sample=$5;export normal_sample
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


#0. Copy mouse_all_exon_mm10.interval_list over to current folder (data_analysis)
cp /home/luol2/lingqi_workspace/Mouse_WES_Pipeline/Primary/mouse_all_exon_mm10.interval_list .

#1. Mutect2 Calling and Filtration
perl -pe 's#TEST_SAMPLE#$ENV{sample}#g;s#PROJECT#$ENV{project}#g;s#SUBJECT#$ENV{subject}#g;s#DATA_PATH#$ENV{data_path}#g;s#NORMAL_SAMPLE#$ENV{normal_sample}#g' run_mutect2_and_Filter.sh | bsub -J mutect_${project}_${subject}_${sample}
jobId=$(bjobs -J mutect_${project}_${subject}_${sample} | awk '{print $1}' | grep -v JOBID)
jobName=mutect_${project}_${subject}_${sample}
message="Mutect2 calling and filtration going......"
checkJobSuccess
sleep 30


#2. Remove the Germline dbsnp variants for BALB/cJ strain
perl -pe 's#TEST_SAMPLE#$ENV{sample}#g;s#PROJECT#$ENV{project}#g;s#SUBJECT#$ENV{subject}#g;s#DATA_PATH#$ENV{data_path}#g' run_remove_BALB_cJ_germline_snp_indels.sh | bsub -J rmsnp_${project}_${subject}_${sample}
jobId=$(bjobs -J rmsnp_${project}_${subject}_${sample} | awk '{print $1}' | grep -v JOBID)
jobName=rmsnp_${project}_${subject}_${sample}
message="Removing BALB/cJ Germline SNPs/INDELs going......."
checkJobSuccess

#3. VEP Annotation, extra manual filtration and variant classification by type
perl -pe 's#TEST_SAMPLE#$ENV{sample}#g;s#PROJECT#$ENV{project}#g;s#SUBJECT#$ENV{subject}#g;s#DATA_PATH#$ENV{data_path}#g' run_VEP_annotation_BALB_cJ.sh | bsub -J VEP_${project}_${subject}_${sample}
jobId=$(bjobs -J VEP_${project}_${subject}_${sample} | awk '{print $1}' | grep -v JOBID)
jobName=VEP_${project}_${subject}_${sample}
message="VEP annotation going......."
checkJobSuccess

# Declare mission completed
echo "Mission Accomplished!"




