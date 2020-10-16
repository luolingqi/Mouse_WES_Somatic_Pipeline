# Diaz Lab Mouse WES Somatic Variant Calling Pipeline (Primary Analysis)
### The legitimate pipeline was built on top of GATK/Mutect against the mm10 reference mouse resource bundle, optimized by serious mutation filters. 

&nbsp;&nbsp;&nbsp;&nbsp;
This is a version-controlled code repository for **Mouse Somatic Variant Calling Pipeline** in development by the Diaz group at MSKCC. The pipeline implementation is specific for execution by the members in the **Diaz group only!!!**


**Primary Lilac Locations for the Pipeline and associated resource bundles**

* Top-level directory of the pipeline: **_/home/luol2/lingqi_workspace/Mouse_WES_Pipeline_** w/sub-directories: 
  - Main pipeline implementation scripts: **_`/Primary`_**
  - Downstream R analysis for variant gain/loss, MAF distribution plots, Mutation signature, etc.: **_`/Seconcary`_**
  
* Primary Mouse Reference Genome (GRCm38/mm10): **_/home/luol2/lingqi_workspace/References/GRCm38_** w/contents:
  - Reference Genome files in diverse format: fasta, dict, bwamem index
  - The Agilent SureSelect exome enrichment target bed file
  - Germline SNP/INDEL recorded by various sources: Mouse Genome Project (MGP, Ensembl, etc.): **_`/Mus_musculus`_**
  
* ENSEMBL VEP resources for both Mouse and Human: **_/home/luol2/lingqi_workspace/vep_data_**

**Primary Scripts for automatic pipeline running:**
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
  
  
