# Diaz Lab Mouse WES Somatic Variant Calling Pipeline <br/> (BALB/c & C57BL/6)

This is a version-controlled code repository for **Mouse Somatic Variant Calling Pipeline** in development by the Diaz group at MSKCC. The pipeline implementation is specific for execution by the members in the **Diaz group only!!!**

The entire pipeline consists of two parts: **`Primary Analysis`** and **`Secondary Analysis`**.  
1. **Primary analysis** takes raw FASTQ files as inputs, follow the procedures described below to generate well-filtered & annotated mutation callset.  
2. **Secondary analysis** picks up the callset and produce multiple summary analyses in terms of both table and plot to cover the routine tasks like betwee-condition variant gain/loss comparison, variant number, VAF distribution by type and mutation signaure deconvolution. In more detail:  

![GitHub Logo](/images/Mouse_WES_Somatic_Mutation_Calling_Pipeline.png)

> The **`primary analysis`** was built on top of both **GATK Best Practice** for somatic calling and a high performance Tumor/Normal variant calling algorithm, **Mutect2**. The sample FASTQ files are aligned and called against the GRCh38 reference genome and its relevant resource bundle. The result set of variants are further manually filtered and VEP annotated to reach a clean, well function-annotated set. 
>
> The **`secondary analysis`** performs downstream summary analysis with the annotated callset. Before its operation, ***RESEARCHERS*** need to *MANUALLY* prepare 3 ***MANIFEST FILES*** to define the analytical scope for 1) `variant comparison by type`; 2) `variant allele frequency (VAF) plot`; 3) `mutation signature`. The template for these manifest files can be found in the tables shown below. A R script summarizes all the secondary analysis outputs and generates a markdown notebook with tables and plots of the final results. 

____________________________________
### **Primary Lilac Locations for the Pipeline and associated Resource Bundles**

* Top-level directory of the pipeline: **_/home/luol2/lingqi_workspace/Mouse_WES_Pipeline_** w/sub-directories: 
  - Main pipeline implementation scripts: **_`/Primary`_**
  - Downstream analysis & plotting for variant gain/loss, MAF distribution plots, Mutation signature, etc.: **_`/Seconcary`_**
  
* Primary Mouse Reference Genome (GRCm38/mm10): **_/home/luol2/lingqi_workspace/References/GRCm38_** w/contents:
  - Reference Genome files in diverse format: fasta, dict, bwamem index
  - The Agilent SureSelect exome enrichment target bed file
  - Germline SNP/INDEL recorded by various sources: Mouse Genome Project (MGP, Ensembl, etc.): **_`/Mus_musculus`_**
  
* ENSEMBL VEP resources for both Mouse and Human: **_/home/luol2/lingqi_workspace/vep_data_**

___________________________________
### **Primary Scripts & manifest files for Automatic Pipeline Running**
  * `step0_manifest_prep.sh` -- it creates the manifest files for secondary analysis
  
  * `step1_preprocessing_simple.sh` -- it takes fastq file and readgroup info as inputs, and runs the following as listed in the table below

  * `step2_Mutect2_VEP_simple.sh` -- it takes the duplicate-removed and BQSR-recalibrated file outputs from `step1_preprocessing_simple.sh`, and runs the following as listed in the table below

  * `step3_gain_loss_VEP_simple.sh` -- it takes a manifest sample comparison file (as shown in the template table below) and all the well annotated/filtered variant files from `step2_Mutect2_VEP_simple.sh`, output the gain/loss for all the pairwise comparisons (e.g. Parental vs Treatment)
  
  * `step4_MutSig_deconstruction.sh` -- it takes a manifest file (as shown in the template table below) recording vcf files for the individual samples and gain/loss comparisons, deconstruct and plot the mutation signatures
    
  * `items_to_compare.txt` -- The manifest file to define parental vs treatment comparisons, as shown in the template `Table 2`
  
  * `items_to_plot.txt` -- The manifest file to define what to plot for VAF distribution, as shown in the template `Table 3`
  
  * `items_to_plot_MutSig.txt` -- The manifest file to define what to plot for mutation signatures, as shown in the template `Table 4`

Step1: Data Quality Checking & Preprocessing  |  Step2: Variant Calling, filtering & Annotation | Step3: Variant Gain/Loss Comparison <br/> (e.g. Parental v.s. Treated) | Step4: Mutation Signature Deconvolution
-------------------------------------------   |  ---------------------------------------------- |  --------------------------------------------------------------- | ----------------------------------------
Quality checking of raw fastq files <br/> **(Fastqc - run_fastqc.sh)**  |  Somatic Variant Calling and filtration <br/> **(GATK Mutect2 - run_mutect2_and_Filter.sh)** | Variant Gain/Loss (Parental vs Treated) <br/> **(run_variants_gain_loss_treatment.vs.Parental.sh)** | Mutation Signature Deconvolution <br/> **(deconstructSigs - run_deconstructsigsmm10.sh)**
Adapter & low quality reads trimming <br/> **(Trimgalore - run_trim_galore.sh)** |  Removing Germline SNP/INDEL variants <br/> **(For BALB/cJ strain - run_remove_BALB_cJ_germline_snp_indels.sh)**
Trimmed fastq to uBAM format conversion <br/> **(required by GATK pipeline - run_fastq_to_uBAM.sh)**  |  ENSEMBL VEP variant annotation & type filtration <br/> **(missense, frameshit, nonsynonymous, etc. - run_VEP_annotation_BALB_cJ.sh)**
BWA MEM alignment to GRCm38/mm10 <br/> **(BWA MEM - run_bwa_mem.sh)**  |  Extra manual filtrations by quality <br/> **(AD, MBQ, MMQ, MPOS5, etc. - run_VEP_annotation_BALB_cJ.sh)**
Alignment quality metrics collection <br/> **(Qualimap - run_qualimap.sh)**  |  
Hybrid selection quality metrics collection <br/> **(GATK HsMetrics - run_CollectHsMetrics.sh)**  |  
Estimate and Apply MarkDuplicate and <br/> Base Quality Score recalibration <br/> **(GATK MarkDuplicate, BQSR  - run_markduplicate.sh)**  |  

**Table 2. Sample Comparison Manifest Template (`User Supplied`)**
**Parental Sample** | **Treated Sample**
------------------- | ------------------
Sample_SW480_P8W_IGO_10212_G_2 | Sample_SW480_P0W_IGO_10212_G_1
Sample_SW480_P8W_IGO_10212_G_2 | Sample_SW480_TMZ8W_IGO_10212_G_3
Sample_SW480_P8W_IGO_10212_G_2 | Sample_SW480_CDDP8W_IGO_10212_G_4
Sample_SW480_P8W_IGO_10212_G_2 | Sample_SW480_COMBO8W_IGO_10212_G_5

**Table 3. VAF Distribution Manifest Template (`Automatically Generated`)**
**cats** | **comparison** | **variant_type** | **files**
-------- | -------------- | ---------------- | ---------
Sample_SW480_CDDP8W_IGO_10212_G_4 |	loss | nonsynonymous | PITT_0522/Sample_SW480_CDDP8W_IGO_10212_G_4_vs_Sample_SW480_P8W_IGO_10212_G_2.VEP.ann.AF0.05_BQ20_MQ50.non_synonymous_variant_only.query
Sample_SW480_CDDP8W_IGO_10212_G_4 |	gain | nonsynonymous | PITT_0522/Sample_SW480_P8W_IGO_10212_G_2_vs_Sample_SW480_CDDP8W_IGO_10212_G_4.VEP.ann.AF0.05_BQ20_MQ50.non_synonymous_variant_only.query
Sample_SW480_CDDP8W_IGO_10212_G_4 | total |	nonsynonymous |	PITT_0522/Sample_SW480_CDDP8W_IGO_10212_G_4/Sample_SW480_CDDP8W_IGO_10212_G_4.aligned.duplicates_marked.recalibrated.targeted_sequencing.somvarius.tumor_only.sorted.rm_dbsnps.VEP.ann.AF0.05_BQ20_MQ50.non_synonymous_variant_only.query

**Table 4. Mutation Signature Manifest Template (`Automatically Generated`)**
**cats** | **comparison** | **files**
-------- | -------------- | ---------
Sample_SW480_COMBO8W_IGO_10212_G_5 | loss | PITT_0522/Sample_SW480_COMBO8W_IGO_10212_G_5_vs_Sample_SW480_P8W_IGO_10212_G_2.VEP.ann.AF0.05_BQ20_MQ50.vcf
Sample_SW480_TMZ8W_IGO_10212_G_3 | loss | PITT_0522/Sample_SW480_TMZ8W_IGO_10212_G_3_vs_Sample_SW480_P8W_IGO_10212_G_2.VEP.ann.AF0.05_BQ20_MQ50.vcf
Sample_SW480_CDDP8W_IGO_10212_G_4 | gain | PITT_0522/Sample_SW480_P8W_IGO_10212_G_2_vs_Sample_SW480_CDDP8W_IGO_10212_G_4.VEP.ann.AF0.05_BQ20_MQ50.vcf
Sample_SW480_COMBO8W_IGO_10212_G_5 | gain | PITT_0522/Sample_SW480_P8W_IGO_10212_G_2_vs_Sample_SW480_COMBO8W_IGO_10212_G_5.VEP.ann.AF0.05_BQ20_MQ50.vcf
Sample_SW480_CDDP8W_IGO_10212_G_4 | total | PITT_0522/Sample_SW480_CDDP8W_IGO_10212_G_4/Sample_SW480_CDDP8W_IGO_10212_G_4.aligned.duplicates_marked.recalibrated.targeted_sequencing.somvarius.tumor_only.sorted.rm_dbsnps.AF0.05_BQ20_MQ50.VEP.ann.vcf

__________________________________
### **Prerequisites for Running the Pipeline**<br/>

* The entire pipeline was built on MSK High Performance Computing (HPC) platform with all the individual building blocks/tools developed in worry-free encapsulated enviroment (Singularity). So, there is little dependency to the system we log in on Lilac, which means, **_anyone with an active Lilac account and basic skill of linux_** can easily run it without any headache of adjusting running environment/parameters.
* The input data structure needs to be organized as following, so that the pipeline can locate the pair-end fastq files in gz format in each sample folder.
```
DATA_PATH/
|-- PROJECT/
|   |-- SUBJECT/ # can be missed
|       |-- SAMPLE/
|           |-- *R1*fastq.gz
|           |-- *R2*fastq.gz
|       |-- NORMAL_SAMPLE/
|           |-- *R1*fastq.gz
|           |-- *R2*fastq.gz
```

__________________________
### **Main Pipeline Usage (Primary & Secondary Analysis!)**

The entire pipeline is split into the following four steps for the sake of efficient debugging and implementation. The operator will need to supervise the successful execution of each step by lauching subsequent one. There is one batch script linked to each of the 4 steps respectively as shown below.  

* step0_manifest_prep.sh
* step1_preprocessing.sh
* step2_Mutect2_VEP.sh
* step3_gain_loss_VEP_simple.sh
* step4_MutSig_deconstruction.sh

Please first copy the above shell scripts into a folder (e.g. data\_analysis) as your working directory to launch these 4 steps! Each step consists of multiple heavy-load tasks, which take time to accomplish. Please be patient to wait till one step ends successfully before kicking off next step. While I don't expect end users to debug any errors which interrupt the pipeline, the pipeline does log the job status and errors in a log file named like "nohup\_step\_\*.log". A note message "Mission Accomplished!" at the end of the log file indicates the success of the step. Make sure you see it before you go to next step.
  
  
```  
  *************************************************************************************
  **From the comparison file, prepare multiple manifest files for secondary analyses.**
  *************************************************************************************
  #########################
  # prepare manifest files 
  #########################
  .USAGE.
  `nohup sh step0_manifest_prep.sh DATA_PATH PROJECT SUBJECT COMPARISON 2>&1 >nohup_step0.log &`
  
  .OPTIONS.
  DATA_PATH  a root directory of the entire study, required.                                 e.g. /home/luol2/lingqi_workspace/Projects/Ben_Projects
  PROJECT    a project name, required.                                                       e.g. WES_mouse_Project_10212_E
  SUBJECT    a subject name if any.                                                          e.g. any name here, if no, just use '.'
  COMPARISON a file listing all pairwise comparisons (Parental vs Treated) , required.       e.g. Comparison_SW480.txt

  *************************************************************************************************************************************** 
  **You can write your own loop to run the following steps for multiple samples. The following USAGE example is for single-sample only.**
  *************************************************************************************************************************************** 
  #################
  # preprocessing 
  #################
  .USAGE.
  `nohup sh step1_preprocessing_simple.sh DATA_PATH PROJECT SUBJECT SAMPLE 2>&1 >nohup_step1_SAMPLE.log &`
  
  .OPTIONS.
  DATA_PATH  a root directory of the entire study, required.             e.g. /home/luol2/lingqi_workspace/Projects/Ben_Projects
  PROJECT    a project name, required.                                   e.g. WES_mouse_Project_10212_E
  SUBJECT    a subject name if any.                                      e.g. any name here, if no, just use '.'
  SAMPLE     a sample name, required.                                    e.g. Sample_CT26CDDP_M1_IGO_10212_E_13
  
  ###################################
  # Mutect2 calling & VEP annotation
  ###################################
  `nohup sh step2_Mutect2_VEP_simple.sh DATA_PATH PROJECT SUBJECT SAMPLE NORMAL_SAMPLE 2>&1 >nohup_step2_SAMPLE.log &`
  
  .OPTIONS.
  Same as the options in step 1, except the name of the normal sample (NORMAL_SAMPLE) is required here
  and the normal sample needs to be in the same SUBJECT folder as SAMPLE
  
  #########################################################################
  # Variant Gain/Loss analysis (e.g. Parental vs Treated) & VEP annotation 
  #########################################################################
  `nohup sh step3_gain_loss_VEP_simple.sh DATA_PATH PROJECT SUBJECT PARENTAL TREATED 2>&1 >nohup_step3_PARENTAL_vs_TREATED.log &`
  
  .OPTIONS.
  Same as the options in step 1, except that you need to specify both PARENTAL and TREATED sample for comparison. Loop through a manifest comparison file as Table 2 above to run the annotated gain/loss analysis as many as you want

  ################################################
  # Mutation Signature analysis (deconstructSigs)
  ################################################
  `nohup sh step4_MutSig_deconstruction.sh DATA_PATH PROJECT MANIFEST 2>&1 >nohup_step4_MutSig.log`
  
  .OPTIONS.
  Same as the option in step 1, except
  MANIFEST      a file recording the samples, variant types whose signatures are to be analyzed, together with file locations to the VCF files           e.g. data_analysis/items_to_plots_MutSig.txt Table 3 above shows an example

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

