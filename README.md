# ![QIAGEN logo](https://www.qiagen.com/sfc/images/qiagen-logo.png) 
**Welcome QDI team!**

Following up on my application for QDI's [Senior Bioinformatics Scientist](https://www.qiagen.com/us/about-us/careers/jobs/details?jobId=17775&jobTitle=Senior%20Bioinformatics%20Scientist%2C%20QDI), Jacquey asked me to provide code demonstrating expertise with Python and R in an omics context. Jacquey also asked that the code be relevant to a pipeline I've written, or benchmarking/testing methods. I think these are _very_ sensible things to ask. Unfortunately, nearly all of the bioinformatics code I've written recently belongs to my former employer. So instead, I've decided to write a small pipeline to address a question that has been in the back of my mind since I came across this MAQC paper relating to [single cell mRNA-seq benchmarking](#https://pubmed.ncbi.nlm.nih.gov/33349700/).

You may be familiar with computation tools used to estimate the proportion of cell types that contributed to a bulk mRNA-seq sample. These tools usually rely on data enriched for a particular cell type or identification of cell type specific gene expression patterns from single cell mRNA-seq data. There are a number of tools available to deconvolute bulk mRNA-seq data, and I've been wondering which methods work the best and if the ratio of cells present in a mixture affects the accuracy of any deconvolution approach. Existing literature on [this topic](https://www.nature.com/articles/s41467-020-19015-1) is valuable, but relies on creating 'pseudo-bulk' mRNA-seq samples and then comparing predictions made using the single cell data that created the pseudo-bulk data.

I would like to know how deconvolution methods perform when single cell and bulk libraries are prepared independently, and if choice of single cell technologies impacts the results. I hope that by looking over this workflow you will get a sense for how I approach problems and if I might be a good fit for your organization. If you would rather get straight to the Python and R, please checkout this [Python](./code/python_stuff/python_cmd_tools.ipynb) and this [R](./code/R_stuff/run_deconvolution.R) code.

# Table of Contents
**[Welcome!](#qiagen-logo)**<br>
**[Introduction](#introduction)**<br>
**[Preparation](#preparation)**<br>
**[Usage Instructions](#usage-instructions)**<br>

# Introduction
The purpose of this code is to download bulk and single cell mRNA-seq data for particular cell lines, use the bulk data sets to create mixtures of cell types as test cases, and then see how the chosen deconvolution tools perform with single cell mRNA-seq data from the same cell lines as test cases. The cell lines being used here are HCC1395 and HCC1395BL. HCC1395 and HCC1395BL are both derived from the same individual. HCC1395 is a breast cancer cell line, and HCC1395BL is a normal B lymphocyte line.

Data corresponding to these lines is available on NCBI's [Short Read Archive](https://www.ncbi.nlm.nih.gov/sra), so the first steps will be to download and process that data. Next, test cases will be created such that each test case has about the same number of total reads, and has known fractions of reads from  HCC1395 and HCC1395L lines. Finally, the results from the deconvolution tools can be obtained, and the results plotted. 

# Preparation

## System requirements
This pipeline requires the workflow language [NextFlow](https://www.nextflow.io/) and the containerization tool [Docker](https://www.docker.com/) to run. I won't provide installation instructions since each tools has extensive documentation. <br>

## Docker Images
Most of this pipeline's Docker images are available on [Docker Hub](https://hub.docker.com/). However, I placed the Python environment and deconvolution tools in the same Docker image, which will need to be built. <br>
This image can be built by running the following code from this directory:

## Environmental Variables
The NextFlow workflow needs it user to supply three file paths corresponding to directories on you local system:
- `fastqDir` - the directory where the raw files are to be downloaded to
- `outputDir` - the directory where results will written to
- `workDir` - the directory where temporary files will be written

# Usage Instructions

To perform each step, you can use the following code.

## Download the desired files and their metadata
```
nextflow run nf-core/fetchngs \
    --input ./SRR_Acc_List.txt \
    -profile docker \
    --force_sratools_download \ # I find this setting more reliable than the default FTP approach
    --outdir $outputDir 
```

## Create cell type mixture samples
```
# needs a docker container with compatible dir structure 
jupyter nbconvert --execute ./code/python_stuff/create_convolution.ipynb
```

## pre-process the bulk mRNA-seq data
```
nextflow run nf-core/rnaseq \
    --input ./bulk_samplesheet.csv \
    --outdir $MAQC_DATA \
    --save_merged_fastq true \
    --with_umi false \
    --skip_bbsplit true \
    --max_memory 80.GB \
    --max_cpus 8 \
    --skip_bbsplit true \
    --genome GRCh38 \
    --save_align_intermeds true \
    --skip_preseq true \
    -profile docker \
    --skip_rseqc
```

## Process the scRNA-seq data
```
nextflow run nf-core/scrnaseq \
    -profile docker \
    --outdir ${MAQC_DATA}/scRNA \
    --input ./scRNA_samplesheet.csv \
    --aligner star \
    --genome GRCh38 \
    --max_memory 80GB \
    --max_cpus 8 \
    --protocol 10XV2

```

# Summary/Results
You can have a look at the summary report this pipeline generates by opening [example_report.html](./example_report.html)

# To-Do
- shrink the QIAGEN logo
- add a License that makes sense