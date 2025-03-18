# DRUG-seq U2OS MoA Box  

[![GitHub:Novartis_DRUG--seq](https://img.shields.io/badge/GitHub-Novartis_DRUG--seq-6699FF.svg)](https://github.com/Novartis/DRUG-seq)
[![Zenodo:DOI](https://img.shields.io/badge/Zenodo-Novartis_DRUG--seq_U2OS_MoABox_Dataset-6699FF.svg)](https://dx.doi.org/10.5281/zenodo.14291446)

## Authorship
Andrea Hadjikyriacou, Chian Yang, Martin Henault, Robin Ge, Leandra Mansur, Alicia Lindeman, Carsten Russ, Steffen Renner, Marc Hild, Jeremy Jenkins, Caroline Gubser-Keller, Jingyao Li, Daniel J. Ho, Marilisa Neri, Frederic D. Sigoillot, Robert Ihry

## Background
DRUG-seq is an open-source high throughput transcriptomics platform. In order to facilitate the development of new computational analysis methods, that require large-scale unbiased transcriptomic data, we are making the DRUG-seq MoA Box dataset publicly available to drive innovation. Please reference the dataset using Ye et. al., 2018 [[1]](#1), and the associated GitHub and Zenodo links provided. We look forward to hearing about the new methods development and hoping to see the future publications that follow!

We profiled a large set of small molecules (N = 4,343), including the previously released Mechanism of Action (MoA) Box (Canham et al., 2020 [[2]](#2) and [Novartis MoA Box GitHub](https://github.com/Novartis/MoaBox)) in the Osteosarcoma U2-OS cell line in a 4-dose response (0.01, 0.1, 1 and 10 microMolar) at one treatment timepoint (24.0hrs), using the high-throughput transcriptomics profiling method DRUG-seq (Ye et al., 2018 [[1]](#1), Li et al., 2022 [[3]](#3), and [Novartis DRUG-seq GitHub](https://github.com/Novartis/DRUG-seq) ). The profiles cover an additional 640 public compound structures part of the evolving Novartis MoA Box collection.

Here we share the following metadata and data:
- plate well-level gene UMI counts
- contrast-level for each compound, dose, timepoint and batch: gene fold change and associated statistics.
- reduced feature space derivatives of the UMI counts (PCA, UMAP) results (and associated pre-filtered and normalized count matrix data files) are provided before/after limma RemoveBatchEffect batch correction.

## Dataset Content Statistics and processing information
### 1. Content description and medata files
  The dataset consists of 52 batches of 3 replicate 384-well plates (156 x 384-well plates), resulting in a total of 156*384 = 59,904 Next Generation Sequencing samples (uniquely identified by biosample_id/external_biosample_id columns in provided metadata) and associated transcriptome profiles. The subset of unused wells (N = 3276) was annotated with the filler well_type=EMPTY and water (H2O) identifier (cmpd_sample_id: EC-27-RY89). There were 56,628 wells=samples effectively used in the data processing post UMI counts generation, representing 4,343 unique compounds tested in 4-dose response and triplicate (one replicate per plate across 3 replicate plates in a batch), except for a few positive controls (Homoharringtonine cmpd_sample_id: EA-18-FP00; BTdCPU cmpd_sample_id: SE-15-AV21; cmpd_sample_id: BD-11-DV28 ) shared between all batches and plates that were represented in triplicate per plate (9 wells per batch) at 10 microMolar and the Reference Control (RC) DMSO wells (cmpd_sample_id: CB-43-EP73) represented in column 23 (16 wells per plate). The positive controls were selected with the expectation they can have a large impact on cell proliferation and transcriptome profile.

  a. Metadata files (3 files) for 59,904 DRUG-seq samples and MoABox public compounds annotation: 
   ***./DRUGseq_U2OS_MoABox_plate_wells_metadata_public.txt*** # contains plate and well-level metadata for the 59,904 DRUG-seq samples including plate and well identifiers, compound identifier (cmpd_sample_id), dose, timepoint. The cmpd_sample_id column allows mapping to compounds metadata in the next file. 
   ***./MoABox_compounds_metadata.txt*** # contains 4,345 compounds (4,343 compounds + DMSO + Water) with annotation of structure (inchi_key and SMILES) and mechanism of action (MoA) when available.
   ***./MoABox_compounds_target_pairs_public.txt*** # contains row-separated compounds (inchi_key) / target (gene, symbol) pairs annotation when available. 
  *For reference*: 
  - Water (H2O; well_type = EMPTY): cmpd_sample_id: EC-27-RY89; inchi_key: 'XLYOFNOQVPJJNP-UHFFFAOYSA-N'; SMILES: 'O' 
  - Dimethyl Sulfoxide (DMSO; well_type=RC): cmpd_sample_id: CB-43-EP73; inchi_key: 'IAZDPXIOMUYVGZ-UHFFFAOYSA-N'; SMILES: 'CS(C)=O' 
  - Wells (NGS samples) are uniquely identified by metadata columns including ‘biosample_id’/‘external_biosample_id’ 
  - Compound samples are uniquely identified by metadata columns including ‘cmpd_sample_id’ 

  *Example of a DRUG-seq 384-well plate content format (384-well plate VH02012944)*
   <img src="./images/example_384well_plate_map.png" width="50%" style="display: block; margin-left: 0;">

  - Tested samples (well_type = ‘SA’) are represented as circles 
  - The series of 16 Reference Control (well_type = ‘RC’) DMSO-treated wells are represented by filled squares in column 23 
  - Non-effective wells (no treatment) are represented by empty squares (well_type = EMPTY) and annotated as the water cmpd_sample_id. 
  - SA wells in column 24 consist of 3 replicates of 10 microMollar for each of 3 positive controls (BD-11-DV28, BTdCPU and Homoharringtonine). 

  b. Genes annotation files 
    ***./drugseq_ensembl_v98_annotation_and_entrez_mapping.RData*** # contains 3 items, drugseq_ensg_v98 has DRUG-seq data gene.ID to Ensembl (Version 98) ENSG gene identifiers mapping and two items have full join and inner join mapping to ENTREZ gene ID and associated gene information (drugseq_ensg_v98_entrez_mapping, drugseq_ensg_v98_entrez_mapping_notNULL) 


### 2. Data processing and data files

The dataset was processed with the Novartis DRUG-seq data processing/analysis pipeline, following the method described in Li et al., 2022 [[3]](#3) with the exception that after the TRUE NULL first step, two robust Reference Control (RC) DMSO wells were selected per plate, hence a total of 6 RC wells per batch, and used in fold change calculations final step. In order to generate the PCA and UMAP views, samples for contrasts defined as ‘active’ over TRUE NULL-defined background activity threshold (global threshold across 52 batches, 5-percentile of DMSO versus DMSO contrasts evaluated: >17 differentially expressed genes).  

  a. plate_well-level gene UMI counts are provided for all 59,904 wells in the file (well_type: SA, RC and EMPTY). Please note the dataset was saved with gzip compression to minimize the file size. We recommend loading the object in R and saving it locally as e.g. uncompressed RData file to allow faster loading for subsequent processing steps. 
   ***./Exp_gzip.RData***

  b. Sample groups and comparisons metadata files 
   ***./robust_RC_ReferenceControl_DMSO_wells.txt*** # set of Reference Control (RC) DMSO wells selected as robust after the TRUE NULL step1 and used in fold change calculations per batch 
   ***./comparisons_metadata.RData*** # contain 2 objects, comparison_group_member_info_public shows the samples groups setup to feed into the comparisons and comparisons_info_public shows metadata for the comparisons. 

  c. Fold Change calculations result files 
   ***./de.RData*** # data.frame aggregating all 17,731contrasts across 52 batches 

### 3. PCA and UMAP dimensionality reduction and resulting files: 

  a. UMI count matrix filtering, before/after batch correction 
   ***./mat.filtered_before_RemoveBatchEffect.RData*** # UMI count matrix after filtering samples (samples for contrasts with a number of differentially expressed genes above the threshold of 17; TRUE NULL 95%ile threshold across batches and select RC control wells ) and genes per the published method, before limma::RemoveBatchEffect was applied on plate_barcodes 
   ***./mat.filtered_after_RemoveBatchEffect.RData*** # UMI count matrix after filtering samples (samples for contrasts with a number of differentially expressed genes above the threshold of 17; TRUE NULL 95%ile threshold across batches and select RC control wells ) and genes per the published method, after limma::RemoveBatchEffect was applied on plate_barcodes 

  b. PCA result files (before/after batch correction) 
   ***./pca_filtered_before_RemoveBatchEffect.RData*** 
   ***./pca_filtered_after_RemoveBatchEffect.RData***    

  c. UMAP results files (2D and 3D, before and after RemoveBatchEffect) 
   ***./umap_2D_coordinates_filtered_before_and_after_RemoveBatchEffect.txt*** 
   ***./umap_3D_coordinates_filtered_before_and_after_RemoveBatchEffect.txt*** 

## Files content excerpts: 
  ***./DRUGseq_U2OS_MoABox_plate_wells_metadata_public.txt*** 
   <img src="./images/DRUGseq_U2OS_MoABox_plate_wells_metadata_public.png" width="90%" style="display: block; margin-left: 0;">

  ***./MoABox_compounds_metadata.txt***
   <img src="./images/MoABox_compounds_metadata.png" width="80%" style="display: block; margin-left: 0;">

  ***./MoABox_compounds_target_pairs_public.txt***
   <img src="./images/MoABox_compounds_target_pairs_public.png" width="70%" style="display: block; margin-left: 0;">

  ***./robust_RC_ReferenceControl_DMSO_wells.txt***
   <img src="./images/robust_RC_ReferenceControl_DMSO_wells.png" width="100%" style="display: block; margin-left: 0;">

  ***./drugseq_ensembl_v98_annotation_and_entrez_mapping.RData*** (two objects) 
   <img src="./images/drugseq_ensembl_v98_annotation_and_entrez_mapping.png" width="50%" style="display: block; margin-left: 0;">

  ***./Exp_gzip.RData***
   Structure of the Exp R object:
    <p style="font-size: 9px;">
    List of 52 batches named ‘1’ to ‘52’, each a list of 3 replicate plates (VH020…), each a list of 2 objects UMI.counts and Annotation 
    List of 52 
    $ 13:List of 3 
      ..$ VH02010536:List of 2 
      .. ..$ UMI.counts: num [1:59594, 1:384] 16 4 0 0 0 0 0 7 28 68 ... 
      .. .. ..- attr(*, "dimnames")=List of 2 
      .. .. .. ..$ : chr [1:59594] "EDC3,grch38_15" "ARHGEF10L,grch38_1" "MTCO3P39,grch38_4" "RBMY2KP,grch38_Y" ... 
      .. .. .. ..$ : chr [1:384] "VH02010536_AACACCTAGT" "VH02010536_AATTGCGATG" "VH02010536_TTGGTCAGTA" "VH02010536_GTTCATTGCC" ... 
      .. ..$ Annotation:'data.frame':	384 obs. of  18 variables: 
      .. .. ..$ batch_id             : int [1:384] 13 13 13 13 13 13 13 13 13 13 ... 
      .. .. ..$ plate_barcode        : chr [1:384] "VH02010536" "VH02010536" "VH02010536" "VH02010536" ... 
      .. .. ..$ plate_index          : chr [1:384] "GCGCTCTA" "GCGCTCTA" "GCGCTCTA" "GCGCTCTA" ... 
      .. .. ..$ well_id              : chr [1:384] "E02" "E08" "L19" "K12" ... 
      .. .. ..$ plate_replicate      : int [1:384] 1 1 1 1 1 1 1 1 1 1 ... 
      .. .. ..$ well_index           : chr [1:384] "AACCGGCGTA" "ATAACGCCTC" "GTTCCGGTGA" "CCTTGTATTC" ... 
      .. .. ..$ col                  : int [1:384] 2 8 19 12 4 11 13 22 14 10 ... 
      .. .. ..$ row                  : int [1:384] 5 5 12 11 12 8 10 12 4 2 ... 
      .. .. ..$ biosample_id         : int [1:384] 1997694 1997700 1997879 1997848 1997864 1997775 1997825 1997882 1997682 1997630 ... 
      .. .. ..$ external_biosample_id: chr [1:384] "CE-17-BF79" "GE-15-QX77" "ID-37-NE71" "WB-39-UD77" ... 
      .. .. ..$ cmpd_sample_id       : chr [1:384] "FC-56-ZH16" "CA-39-ZJ18" "BA-44-WG29" "ED-32-JY57" ... 
      .. .. ..$ well_type            : chr [1:384] "SA" "SA" "SA" "SA" ... 
      .. .. ..$ cell_line_name       : chr [1:384] "U-2 OS" "U-2 OS" "U-2 OS" "U-2 OS" ... 
      .. .. ..$ cell_line_ncn        : chr [1:384] "FH55-48QE" "FH55-48QE" "FH55-48QE" "FH55-48QE" ... 
      .. .. ..$ concentration        : num [1:384] 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 ... 
      .. .. ..$ unit                 : chr [1:384] "uM" "uM" "uM" "uM" ... 
      .. .. ..$ hours_post_treatment : chr [1:384] "24.0" "24.0" "24.0" "24.0" ... 
      .. .. ..$ Sample               : chr [1:384] "FC-56-ZH16_0.01uM_24.0hr_13" "CA-39-ZJ18_0.01uM_24.0hr_13" "BA-44-WG29_0.01uM_24.0hr_13" "ED-32-JY57_0.01uM_24.0hr_13" ...       
    </p>

  ***./mat.filtered_before_RemoveBatchEffect.RData*** 
   <img src="./images/mat.filtered_before_RemoveBatchEffect.png" width="60%" style="display: block; margin-left: 0;">

  ***./mat.filtered_after_RemoveBatchEffect.RData*** 
   <img src="./images/mat.filtered_after_RemoveBatchEffect.png" width="60%" style="display: block; margin-left: 0;">

  ***./pca_filtered_before_RemoveBatchEffect.RData*** 
   <img src="./images/pca_filtered_before_RemoveBatchEffect.png" width="40%" style="display: block; margin-left: 0;">

  ***./pca_filtered_after_RemoveBatchEffect.RData***    
   <img src="./images/pca_filtered_after_RemoveBatchEffect.png" width="40%" style="display: block; margin-left: 0;">

  ***./umap_2D_coordinates_filtered_before_and_after_RemoveBatchEffect.txt*** 
  ***./umap_3D_coordinates_filtered_before_and_after_RemoveBatchEffect.txt*** 
  ![umap_3D_coordinates_filtered_before_and_after_RemoveBatchEffect](./images/umap_3D_coordinates_filtered_before_and_after_RemoveBatchEffect.png)
  ***./comparisons_metadata.RData*** (2 objects)
  - comparison_group_member_info_public 
   <img src="./images/comparison_group_member_info_public.png" width="90%" style="display: block; margin-left: 0;">
  - comparisons_info_public 
   <img src="./images/comparisons_info_public.png" width="90%" style="display: block; margin-left: 0;">
  (comparison_group_id1/2 map to comparison_group_id in comparison_group_member_info_public)
  ***./de.RData*** 
  <img src="./images/de.png" width="90%" style="display: block; margin-left: 0;">
 
### References
<a id="1">[1]</a>  Ye C, Ho DJ, Neri M, Yang C, Kulkarni T, Randhawa R, Henault M, Mostacci N, Farmer P, Renner S, Ihry R, Mansur L, Keller CG, McAllister G, Hild M, Jenkins J, Kaykas A. 
DRUG-seq for miniaturized high-throughput transcriptome profiling in drug discovery. Nat Commun. 2018 Oct 17;9(1):4307. PMID: 30333485; PMCID: PMC6192987. 
[![DOI:10.1101/2021.01.08.425840](https://img.shields.io/badge/DOI-10.1038/s41467--018--06500--x-B31B1B.svg)](https://doi.org/10.1038/s41467-018-06500-x)

<a id="2">[2]</a>  Canham SM, Wang Y, Cornett A, Auld DS, Baeschlin DK, Patoor M, Skaanderup PR, Honda A, Llamas L, Wendel G, Mapa FA, Aspesi P Jr, Labbé-Giguère N, Gamber GG, Palacios DS, Schuffenhauer A, Deng Z, Nigsch F, Frederiksen M, Bushell SM, Rothman D, Jain RK, Hemmerle H, Briner K, Porter JA, Tallarico JA, Jenkins JL. Systematic Chemogenetic Library Assembly. Cell Chem Biol. 2020 Sep 17;27(9):1124-1129. doi: 10.1016/j.chembiol.2020.07.004. Epub 2020 Jul 23. PMID: 32707038. 
[![DOI:10.1016/j.chembiol.2020.07.004](https://img.shields.io/badge/DOI-10.1016/j.chembiol.2020.07.004-B31B1B.svg)](https://doi.org/10.1016/j.chembiol.2020.07.004)

<a id="3">[3]</a> Li J, Ho DJ, Henault M, Yang C, Neri M, Ge R, Renner S, Mansur L, Lindeman A, Kelly B, Tumkaya T, Ke X, Soler-Llavina G, Shanker G, Russ C, Hild M, Gubser Keller C, Jenkins JL, Worringer KA, Sigoillot FD, Ihry RJ. DRUG-seq Provides Unbiased Biological Activity Readouts for Neuroscience Drug Discovery. ACS Chem Biol. 2022 Jun 17;17(6):1401-1414. doi: 10.1021/acschembio.1c00920. Epub 2022 May 4. PMID: 35508359; PMCID: PMC9207813. 
[![DOI:10.1021/acschembio.1c00920](https://img.shields.io/badge/DOI-10.1021/acschembio.1c00920-B31B1B.svg)](https://doi.org/10.1021/acschembio.1c00920)

### GitHub links 
[![GitHub:Novartis_MoaBox](https://img.shields.io/badge/GitHub-Novartis_MoaBox-6699FF.svg)](https://github.com/Novartis/MoaBox)
[![GitHub:Novartis_DRUG--seq](https://img.shields.io/badge/GitHub-Novartis_DRUG--seq-6699FF.svg)](https://github.com/Novartis/DRUG-seq)
[![GitHub:Novartis_DRUG--seq_U2OS_MoABox_Dataset](https://img.shields.io/badge/GitHub-Novartis_DRUG--seq_U2OS_MoABox_Dataset-6699FF.svg)](https://github.com/Novartis/DRUG-seq/tree/main/data/Novartis_drugseq_U2OS_MoABox)
### Zenodo link: 
[![Zenodo:DOI](https://img.shields.io/badge/Zenodo-Novartis_DRUG--seq_U2OS_MoABox_Dataset-6699FF.svg)](https://dx.doi.org/10.5281/zenodo.14291446)