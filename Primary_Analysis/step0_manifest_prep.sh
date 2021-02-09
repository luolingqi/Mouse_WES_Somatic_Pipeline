######################################################################################################################
# The script prepares multiple manifest files based on experiment design from researchers                                                
# It takes a pairwise comparison file as the input, and output the following two manifest files for VAF and Mutsig plotting respectively: 
#       * items_to_plots.txt                                                                                                               
#       * items_to_plots_MutSig.txt                                                                                                         
############################################################################################################################

data_path=$1
project=$2
subject=$3
comparison=$4
python ${data_path}/${project}/data_analysis/prepare_manifest.py ${data_path}/${project}/${comparison} ${subject} ${data_path}/${project}/items_to_plots.txt ${data_path}/${project}/items_to_plots_MutSig.txt


# Declare mission completed
echo "Mission Accomplished!"


