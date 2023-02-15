# Integrated single cell analysis reveals co-evolution of malignant B cells and the tumor microenvironment in transformed follicular lymphoma

<br></br>
## Table of Contents
1. [Abstract](#abstract)
2. [Patient cohort](#patient-cohort)
3. [Code structure](#code-structure)
4. [Publication](#publication)
5. [Acknowledgement](#acknowledgement)

<br></br>
## Abstract
 * Histological transformation of follicular lymphoma to aggressive form occurs 2-3% per year with poor outcome. Divergent evolution and an altered tumour-microenvironment (TME) have been implicated during transformation. However, phenotypic consequences of this evolution and its implication in TME remain unknown. To address this, we performed single cell whole genome and whole transcriptome sequencing of paired pre/post transformation patient samples. We further performed scWTS of additional samples from patients without transformation. Our analysis revealed evolutionary dynamics of transformation at unprecedented resolution. Integration of scWGS and scWTS identified pathways upregulated during evolution. scWTS analysis revealed a shifting TME landscape, with an exhausted signature emerging during transformation. Using multi-color immunofluorescence we transferred these findings to a TME-based transformation biomarker, subsequently validated in independent pretreatment cohorts. Taken together, our results provide a comprehensive view of the combined genomic and phenotypic evolution of malignant cells during transformation, and shifting cross-talk between malignant cells and TME. 

<br></br>
## Patient cohort
 * 11 patients diagnosed with FL but no transformation. 
 * 11 patients diagnosed with FL and later transformed to DLBCL (two biopsies). 
 * 2 patients without FL for reactive lymph node biopsy.
   
<br></br>
## Code structure
1. Data preprocessing
** 10x BCL files was processed with [Cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) (v3.0.2) and then processed with the preprocess_scRNA.R script to generate a filtered SingleCellExperiment object. BCR data was processed with [Cellranger vdj](https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/using/vdj) to generate clonotype information. DLP data was processed with [single cell pipeline](https://single-cell-pipeline.readthedocs.io/en/latest/) to generate phylogenetic trees and clones.
2. Differential expression
** DE was performed with the de.R script.
3. [Clonealign](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1645-z)
** Clonealign was performed using the scripts in clonealign folder.
4. IHC data analysis was performed using the scripts in IHC_analysis folder.
5. Survival analysis was performed using the survival.R script.

<br></br>
## Publication
 <I>[Integrated single cell analysis reveals co-evolution of malignant B cells and the tumor microenvironment in transformed follicular lymphoma](https://www.biorxiv.org/content/10.1101/2022.11.17.516951v1)</I>

<br></br>
## Acknowledgement
 * Roth lab, University of British Columbia.