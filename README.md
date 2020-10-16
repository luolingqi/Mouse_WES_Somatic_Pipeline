# Diaz Lab Mouse WES Somatic Variant Calling Pipeline (Primary Analysis)

&nbsp;&nbsp;&nbsp;&nbsp;
This is a version-controlled code repository for **Mouse Somatic Variant Calling Pipeline** in development by the Diaz group at MSKCC. The pipeline implementation is specific for execution by the members in the **Diaz group only!!!**

&nbsp;&nbsp;&nbsp;&nbsp;
The entire pipeline was built on top of GATK/Mutect against the GRcm38/mm10 reference mouse resource bundle, further reinforced by a few more manual variant filtration and legimate ENSEMBL VEP annotation

**Primary Lilac Locations for the Pipeline and associated resource bundles**

* Top-level directory of the pipeline: **_/home/luol2/lingqi_workspace/Mouse_WES_Pipeline_** w/sub-directories: 
  - Main pipeline implementation scripts: **_`/Primary`_**
  - Downstream R analysis for variant gain/loss, MAF distribution plots, Mutation signature, etc.: **_`/Seconcary`_**
  
* Primary Mouse Reference Genome (GRCm38/mm10): **_/home/luol2/lingqi_workspace/References/GRCm38_** w/contents:
  - Reference Genome files in diverse format: fasta, dict, bwamem index
  - The Agilent SureSelect exome enrichment target bed file
  - Germline SNP/INDEL recorded by various sources: Mouse Genome Project (MGP, Ensembl, etc.): **_`/Mus_musculus`_**
  
* ENSEMBL VEP resources for both Mouse and Human: **_/home/luol2/lingqi_workspace/vep_data_**

**Primary Scripts for automatic pipeline running**
  * `step1_preprocessing.sh` -- it takes fastq file and readgroup info as inputs, and run the following tasks: 
    1) Quality checking of raw fastq files using **Fastqc**
    2) Adapter & low quality reads trimming using **Trimgalore**
    3) Trimmed fastq to uBAM format conversion required by GATK pipeline
    4) BWA MEM alignment to Reference Mouse Genome (mm10)
    5) Alignment quality metrics collection using **Qualimap**
    6) Hybrid selection quality metrics collection using **GATK HsMetrics**
    7) Estimate and Apply MarkDuplicate and Base Quality Score recalibration using **GATK**
    
  * `step2_Mutect2_VEP.sh` -- it takes the duplicate-removed and BQSR-recalibrated file from `step1_preprocessing.sh`, and run the following tasks:
    1) Somatic Variant Calling and filtration using **GATK Mutect2**
    2) Removing Germline SNP/INDEL variants of the BALB/cJ strain
    3) ENSEMBL VEP variant annotation, and extra manual filtrations(MAF, MPOS5, etc.)
    
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


**Main Pipeline Usage (Primary Analysis Only!)**<br/>
&nbsp;&nbsp;&nbsp;&nbsp;
The pipeline is automatically implemented as two bash scripts to cover two steps, as detailed above: 
    1) step1_preprocessing.sh
    2) step2_Mutect2_VEP.sh
    
&nbsp;&nbsp;&nbsp;&nbsp;    
  Each step consists of multiple heavy-load tasks, which take long time to accomplish. Please be sure to wait till the step1 to be successfully accomplished before the step2 could be launched. While I don't expect end users to trouble shoot any errors that cause interruption of the pipeline, the pipeline does log the  running status and errors in a log file named like "nohup_step1_*.log" or "nohup_step2_*.log". A note message "Mission Accomplished!" at the end of the log file indicates the success of the step. Make sure you see it before you go to next step.
  
```  
  # preprocessing
  .USAGE.
  nohup sh step1_preprocessing.sh DATA_PATH PROJECT SUBJECT SAMPLE READGROUP 2>&1 >nohup_step1_SAMPLE.log &
  
  .OPTIONS.
  DATA_PATH  a root directory of the entire study, required.             e.g. /home/luol2/lingqi_workspace/Projects/Ben_Projects
  PROJECT    a project name, required.                                   e.g. WES_mouse_Project_10212_E
  SUBJECT    a subject name if any.                                      e.g. any name here, if no, just use '.'
  SAMPLE     a sample name, required.                                    e.g. Sample_CT26CDDP_M1_IGO_10212_E_13
  READGROUP  a read group ID for the sample in fastq file, required.     e.g. HMVGLDRXX.2
  
  
  # Mutect2 calling & VEP annotation
  nohup sh step2_Mutect2_VEP.sh DATA_PATH PROJECT SUBJECT SAMPLE 2>&1 >nohup_step2_SAMPLE.log &
  
  .OPTIONS.
  Same as the options in step 1, except READGROUP is not required here
  
```

## Versioning
For the versions available, see the tags on this repository.

## Authors
* Lingqi Luo - initial drafting - luolingqi
See also the list of contributors who participated in this project.

## License
This project is licensed by the Diaz laboratory at MSKCC. Only Diaz group is authorized to use this pipeline.

## Acknowledgements
* Team DiazLab @ MSKCC

