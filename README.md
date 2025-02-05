# SeQuoIA


---
## Paper information

**Title:** SeQuoIA: a novel single-cell BCR sequencing analysis pipeline for the identification of BCR-driven selection mechanisms 

**Authors:** Eglantine HECTOR &  Pierre MILPIED 

Aix-Marseille Université, INSERM, CNRS, Centre d'Immunologie de Marseille-Luminy, Marseille (CIML), France. 


**Quick summary:**
SeQuoIA is a snakemake pipeline performing single cell BCR repertoire processing, clonotype assignement and phylogeny reconstuction. 
Additionally, a model of SHM is performed and allows the inference of selection pressure at the codon or sequence level. 
The associated paper describes the methodology, the benchmark of phylogeny reconstruction and selection quantification, as well as new biological results obtained with the application of the pipeline. 

**Correspondance:** hector@ciml.univ-mrs.fr (until june 2025), milpied@ciml.univ-mrs.fr

**Citation:** 
If you find SeQuoIA useful for your publication, please cite:


---
## Contributors (code)

* Eglantine HECTOR (hector@ciml.univ-mrs.fr)

---
## Usage and warnings 

The user can provide his own scBCRseq inputs and rerun an analysis from the paper (Sprumont et al.).
In the future, we plan to provide large cohort human datasets with precomputed selection scores for interactive data exploration in R. 
For input format, the user is referred to the paper supplementary table 2.

SeQuoIA is a pipeline and not a software, and may be not user friendly for non bioinformaticians. 
We also warn potential users that computational time may be long (2 days to 3 weeks), depending on dataset size and complexity and computational resources. 
The BCRseq community is welcome to reuse the concepts introduced in the paper or pieces of the code, although they are not optimal in terms of resource usage.  

---
## Depedencies

In order to prepare the environment for analysis execution, it is required to:
- Install Singularity and Docker
- Install Snakemake

Below you will find detailed instruction for each of these steps.



### Snakemake


You need to install Snakemake to run the complete analysis workflow. Use your preferred solution : https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

### Singularity  
You need to install Singularity v2.6 on your system to run the complete analysis. Follow the instructions here : https://sylabs.io/guides/2.6/admin-guide/
Several docker environments stored in a .sif format will be called by snakemake. 

### Docker images 

The environment for custom analyses in an interactive Rstudio environent can be compiled from a docker file in the 02_Container folder. 
Go to YOUR_WORKING_DIRECTORY/SeQuoIA_Repertoire_Analysis/02_Container/Data_Analysis and run:

```
docker build -t rstudio_sequoia_data_analysis  -f rstudio_SeQuoIA_data_analysis.dockerfile .
docker run --name rstudio_sequoia_data_analysis -d -p 8787:8787 -v /mnt:/mnt rstudio_sequoia_data_analysis
```

---
## SeQuoIA Installation
To use SeQuoIA, clone the github repository and set the WORKING_DIR environment variable in the snakekefile to the path to this folder. 
Once done, you should obtain the following subfolder structure, each of them containing several files.


```
    SeQuoIA_Repertoire_Analysis
    ├── 00_Rawdata
    │   ├── 01_VDJ_BCR_contigs
    │   └── 02_scRNAseq
    ├── 01_Reference
    │   ├── Genomes_Alignements
    │   ├── IMGT_Coordinates
    │   └── Substitution_Models
    ├── 02_Container
    │   ├── BCR_preprocessing
    │   ├── Trees
    │   ├── IGoR
    │   ├── Selection_R
    │   └── Data_Analysis
    ├── 03_Script
    │   ├── 01_BCRseq_processing
    │   ├── 02_Phylogeny
    │   ├── 03_IgoR_VDJ_optional
    │   └── 04_Analysis
    ├── 04_Workflow
    │   └── config_file.yaml    
    ├── 05_Output


```

The scripts are listed and described in Supplementary Table 1 from the paper.
The output folder is empty at this stage and will be automatically generated by SeQuoIA. 

---


## Reproducible example

Provided config file parameters are set to run a reproducible example. The user can regenerate QC figure 3E from the paper. 
Otherwise, load your data as explained in tutorial section. 


## Tutorial 


### 1) Loading input data 

Typical inputs feature seurat metadata object and a cellranger VDJ outputs. 
The input format is detailed in the paper supplementary tables and can be easily created from any single cell technology. 

These inputs are to be loaded in the 00_Rawdata folder as follows:

```
SeQuoIA_Repertoire_Analysis
    ├── 00_Rawdata
    │   ├── 01_VDJ_BCR_contigs
    │   │	└──YOUR_DATASET_NAME
    │   │		└──SUBDATASET
    │   │			├──PREFIX.fasta
    │   │			└──PREFIX_annotations.csv
    │   └── 02_scRNAseq
    │   	└──YOUR_DATASET_NAME
    │			└──YOUR_METADATA_TABLE.txt/csv/


```

NB: metadata file name should contain "metadata" pattern (case insensitive)

### 2) Parameters configuration 

Set you dataset, subdataset, input prefix (typically all_contig or filtered_contig as in the outputs of CellrangerVDJ but anything could work) in the config files, as in the folder and filenames.
Fill in ressource and organism in the config file (04_Workflow). 
The donor column is optional. If matching a column name in metadata, clonotype assingment will be run in a donor-based approach. 
Your will find examples and precisions on the different options as commentaries in the config file 



### 3) Preprocessing & phylogeny reconstruction 

This step is mandatory for subsequent analyses. 
We remind the user that the directory variable is to be set in all snakefiles before running the pipeline. 
The first snakemake script can be run with following command line: 

```
cd YOUR_PATH/SeQuoIA_Repertoire_Analysis
snakemake --use-singularity  --singularity-args "-B /mnt/Volume:/mnt/Volume -B /run/user:/run/user" -s snakefile_BCRseq_annotation_inference_phylogeny.yaml  -j 1 --latency-wait 30 --keep-incomplete 
```
Additional options can be included in the command line. For more information, refer to snakemake dedicated page: 
https://snakemake.readthedocs.io/en/stable/executing/cli.html

The preprocessing step will generate:
* preprocessed and filtered sequence tables in the AIRR format, to be found in YOUR_WORKING_DIRECTORY/SeQuoIA_Repertoire_Analysis/05_Output/01_BCRseq_processing/01_Generating_airr_BCR_filtered_files/YOUR_DATASET/YOUR_SUBDATASET
* Clonotype assignement html reports and tables, to be found in YOUR_WORKING_DIRECTORY/SeQuoIA_Repertoire_Analysis/05_Output/01_BCRseq_processing/02_Clonotype_assignment/YOUR_DATASET/YOUR_SUBDATASET
* NCA/germline sequences split by clonotype, and a summary table of each sequence featuring AIRR data, clonotype name and reconstructed ancestral sequencesto be found in YOUR_WORKING_DIRECTORY/SeQuoIA_Repertoire_Analysis/05_Output/01_BCRseq_processing/03_UCA_MRCA_inference_annotation/YOUR_DATASET/YOUR_SUBDATASET

Reconstructed trees are to be found in the following directory: YOUR_WORKING_DIRECTORY/SeQuoIA_Repertoire_Analysis/05_Output/02_Phylogeny/YOUR_DATASET/YOUR_SUBDATASET
Subfolders are created for each clonotype, with subfolder names featuring unique CloneID, clone size and VDJ composition to ease data exploration. 
In each subfolder, you will find: 
* html report summarizing phylogenetic reconstruction as well as tree statistics 
* fasta files of BCR sequences, including observed nodes and intermediary sequences 
* a summary table of selected trees in newick format and their quality scores 

### 4) Optional IgoR application for clonal selection 

This step is optional and calls another package (IGOr) for VDJ combination analysis. 
As in previous section. Change the -s snakefile variable with the corresponding file name "snakefile_BCRseq_VDJ_combinations.yaml". 
A dataframe summarizing VDJ combination probabilities for each clonotype will be generated in the following folder: YOUR_WORKING_DIRECTORY/SeQuoIA_Repertoire_Analysis/05_Output/04_Analysis/04_VDJ_usage/YOUR_DATASET/YOUR_SUBDATASET

### 5) Mutation Modeling + score calculation 

See "Preprocessing & phylogeny reconstruction" section. The corresponding snakefile name corresponds to "snakefile_BCRseq_selection_scores.yaml" in the main directory. 


The main output is generated in the following location: YOUR_WORKING_DIRECTORY/SeQuoIA_Repertoire_Analysis/05_Output/04_Analysis/04_VDJ_usage/YOUR_DATASET/YOUR_SUBDATASET/All_Samples_Merged_Scores.txt
In the same folder, you will find an interactive html report summarizing clonotype composition, local selection scores summary, and selection score distribution across categories defined in the optional variable METADATA_COLUMNS of the config file. These analyses are basic and performed in unfiltered data. Custom analysis is thus highly recommended. 

Other outputs can be found in the following locations: 
* YOUR_WORKING_DIRECTORY/SeQuoIA_Repertoire_Analysis/05_Output/04_Analysis/03_Selection/YOUR_DATASET/YOUR_SUBDATASET : summary table of theoretical probabilities and QCs
* YOUR_WORKING_DIRECTORY/SeQuoIA_Repertoire_Analysis/05_Output/04_Analysis/03_Selection/03_Selection_Scores_Local/YOUR_DATASET/YOUR_SUBDATASET : html report with intraclonal scores, codon selection scores plots 


### 7) Downstream analyses


Custom analyses can be performed with the following outputs:
* All_Samples_Merged_Scores.txt: this table summarises Clonotype, Tree relationships, selection scores and other BCR features (Isotype, mutation load...) for each cell (rows)
* Tree folder: created automatically and located in the following directory YOUR_WORKING_DIRECTORY/SeQuoIA_Repertoire_Analysis/05_Output/02_Phylogeny/YOUR_DATASET_NAME/YOUR_SUBDATASET. 
* Pgens from optional IGOr PGen_VDJ_summary.txt: this output from the optional IGOr module summarises Pgen per clonotype (rows)

The scripts used for the generation of paper analyses and figures can be found in the 06_Paper_Figures folder that is not part of the pipeline but provided as a template for custom analyses. 

Future developments may include a R library for data exploration and the provision of human datasets for selection mechanism exploration. Please be patient! 


### Piping SeQuoIA modules 

All steps (except user custom downstream analyses) can be piped into a single command line, thanks to the AllPartsSnakemake.bash script provided in the main directory. 

```
./allPartsSnakemake.bash > SeQuoIA_full_log.txt 
```

Uncomment optional steps if necessary. 


