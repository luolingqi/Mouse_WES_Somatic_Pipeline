# Diaz Lab Mouse WES Somatic Variant Calling Pipeline

&nbsp;&nbsp;&nbsp;&nbsp;
This is a version-controlled code repository for **Mouse Somatic Variant Calling Pipeline** in development by the Diaz group at MSKCC. The pipeline implementation is specific for execution by the members in the **Diaz group only!!!**

&nbsp;&nbsp;&nbsp;&nbsp;
The entire pipeline was built on top of GATK/Mutect against the GRcm38/mm10 reference mouse resource bundle, further reinforced by a few more manual variant filtration and legimate ENSEMBL VEP annotation

![GitHub Logo](/images/Mouse_WES_Somatic_Mutation_Calling_Pipeline.png)


**Primary Lilac Locations for the Pipeline and associated Resource Bundles**

* Top-level directory of the pipeline: **_/home/luol2/lingqi_workspace/Mouse_WES_Pipeline_** w/sub-directories: 
  - Main pipeline implementation scripts: **_`/Primary`_**
  - Downstream R analysis for variant gain/loss, MAF distribution plots, Mutation signature, etc.: **_`/Seconcary`_**
  
* Primary Mouse Reference Genome (GRCm38/mm10): **_/home/luol2/lingqi_workspace/References/GRCm38_** w/contents:
  - Reference Genome files in diverse format: fasta, dict, bwamem index
  - The Agilent SureSelect exome enrichment target bed file
  - Germline SNP/INDEL recorded by various sources: Mouse Genome Project (MGP, Ensembl, etc.): **_`/Mus_musculus`_**
  
* ENSEMBL VEP resources for both Mouse and Human: **_/home/luol2/lingqi_workspace/vep_data_**

**Primary Scripts for Automatic Pipeline Running**
  * `step1_preprocessing.sh` -- it takes fastq file and readgroup info as inputs, and runs the following as listed in the table below

  * `step2_Mutect2_VEP.sh` -- it takes the duplicate-removed and BQSR-recalibrated file outputs from `step1_preprocessing.sh`, and runs the following as listed in the table below
    
Step1: Data Quality Checking & Preprocessing  |  Step2: Variant Calling, filtering & Annotation
-------------------------------------------   |  ----------------------------------------------
Quality checking of raw fastq files <br/> **(Fastqc - run_fastqc.sh)**  |  Somatic Variant Calling and filtration <br/> **(GATK Mutect2 - run_mutect2_and_Filter.sh)**
Adapter & low quality reads trimming <br/> **(Trimgalore - run_trim_galore.sh)** |  Removing Germline SNP/INDEL variants <br/> **(For BALB/cJ strain - run_remove_BALB_cJ_germline_snp_indels.sh)**
Trimmed fastq to uBAM format conversion <br/> **(required by GATK pipeline - run_fastq_to_uBAM.sh)**  |  ENSEMBL VEP variant annotation & type filtration <br/> **(missense, frameshit, nonsynonymous, etc. - run_VEP_annotation_BALB_cJ.sh)**
BWA MEM alignment to GRCm38/mm10 <br/> **(BWA MEM - run_bwa_mem.sh)**  |  Extra manual filtrations by quality <br/> **(AD, MBQ, MMQ, MPOS5, etc. - run_VEP_annotation_BALB_cJ.sh)**
Alignment quality metrics collection <br/> **(Qualimap - run_qualimap.sh)**  |  
Hybrid selection quality metrics collection <br/> **(GATK HsMetrics - run_CollectHsMetrics.sh)**  |  
Estimate and Apply MarkDuplicate and <br/> Base Quality Score recalibration <br/> **(GATK MarkDuplicate, BQSR  - run_markduplicate.sh)**  |  

    
**Prerequisites for Running the Pipeline**<br/>

* The entire pipeline was built on MSK High Performance Computing (HPC) platform with all the individual building blocks/tools developed in worry-free encapsulated enviroment (Singularity). So, there is little dependency to the system we log in on Lilac, which means, **_anyone with an active Lilac account and basic skill of linux_** can easily run it without any headache of adjusting running environment/parameters.
* The input data structure needs to be organized as following, so that the pipeline can locate the pair-end fastq files in gz format in each sample folder.
```
DATA_PATH/
|-- PROJECT/
|   |-- SUBJECT/ # can be missed
|       |-- SAMPLE/
|           |-- *R1*fastq.gz
|           |-- *R2*fastq.gz
```


**Main Pipeline Usage (Primary Analysis Only!)**

&nbsp;&nbsp;&nbsp;&nbsp;
The pipeline is automatically implemented in two steps using the following two batch scripts respectively: 
* step1_preprocessing.sh
* step2_Mutect2_VEP.sh

&nbsp;&nbsp;&nbsp;&nbsp;
Please first copy all the shell scripts into a folder as your working directory and launch the 2 steps from there! Each step consists of multiple heavy-load tasks, which take long time to accomplish. Please be sure to wait till the step1 to be successfully accomplished before the step2 could be launched. While I don't expect end users to trouble shoot any errors that cause interruption of the pipeline, the pipeline does log the  running status and errors in a log file named like "nohup_step1_*.log" or "nohup_step2_*.log". A note message "Mission Accomplished!" at the end of the log file indicates the success of the step. Make sure you see it before you go to next step.
  
```  
  # preprocessing
  .USAGE.
  nohup sh step1_preprocessing.sh DATA_PATH PROJECT SUBJECT SAMPLE 2>&1 >nohup_step1_SAMPLE.log &
  
  .OPTIONS.
  DATA_PATH  a root directory of the entire study, required.             e.g. /home/luol2/lingqi_workspace/Projects/Ben_Projects
  PROJECT    a project name, required.                                   e.g. WES_mouse_Project_10212_E
  SUBJECT    a subject name if any.                                      e.g. any name here, if no, just use '.'
  SAMPLE     a sample name, required.                                    e.g. Sample_CT26CDDP_M1_IGO_10212_E_13
  
  
  # Mutect2 calling & VEP annotation
  nohup sh step2_Mutect2_VEP.sh DATA_PATH PROJECT SUBJECT SAMPLE NORMAL_SAMPLE 2>&1 >nohup_step2_SAMPLE.log &
  
  .OPTIONS.
  Same as the options in step 1, except the name of the normal sample (NORMAL_SAMPLE) is required here
  and the normal sample needs to be in the same SUBJECT folder as SAMPLE
  
```

## Versioning
For the versions available, see the [tags on this repository](https://github.com/luolingqi/Mouse_WES_Somatic_Primary_Analysis/releases/tag/v0.2-alpha).

## Authors
* **Lingqi Luo, PhD** - initial drafting - [luolingqi](https://github.com/luolingqi) <br/>
See also the list of [contributors](https://github.com/luolingqi/Mouse_WES_Somatic_Primary_Analysis/contributors) who participated in this project.

## License
This project is licensed by the Diaz laboratory at MSKCC. Only Diaz group is authorized to use this pipeline.

## Acknowledgements
* Team DiazLab @ MSKCC

