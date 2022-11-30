# MIMOSA

__MIMOSA__, or __MWAS Imputing Methylome Obliging Summary-level mQTLs and Associated LD matrices__, is a set of models that substantially improve the prediction accuracy of DNA methylation and downstream MWAS power through the use of a large, summary-level mQTL dataset provided by the Genetics of DNA Methylation Consortium (GoDMC).  Here, we provide a tutorial to download the models and run MWAS with them.  

### Citation

If you use the MIMOSA models, please cite

> Melton, H. J., Zhang, Z., Wu, L., & Wu, C. (2022). MIMOSA: A resource consisting of improved methylome imputation
models increases power to identify CpG site-phenotype associations. Under Review.

## Step 1: Download DNAm prediction models

The MIMOSA models are available on Zenodo (DOI: 10.5281/zenodo.7325055).  The files included are _cg********.rds_, which each consist of a set of weights for SNPs used in predicting DNAm at the CpG site cg*********.  Once you load one of these files into R, you'll have a list with elements 1) a dataframe with (among other things) mQTL p-value, rsID, chromosome, position (of SNP), a1, a2, CpG site, and position of CpG site; 2) a set of DNAm prediction weights corresponding to the SNPs in 1); 3) the prediction accuracy (R^2) of the model for the CpG site in the test dataset (which comes from the Framingham Heart Study).  Once you have the weights, you're halfway to running your MWAS.

## Step 2: Prepare GWAS Summary Statistics

The GWAS summary statistics you use for your MWAS will need to be in the correct format to work with MIMOSA.  Please follow these conventions:

* Split the raw data into 22 subfiles, one for each chromosome.  If the trait you're considering is asthma (AST), then name these subfiles AST-1.sumstats, ..., AST-22.sumstats.
* Ensure the header of the subfiles are in ALL CAPS.  At minimum, you should have columns for SNP (where this is rsIDs), CHR, A1, A2, Z.
* You can manually process the summary statistics, or you can use the provided APSS.R function

### APSS.R

_APSS.R_ is an interactive R function that helps you easily process and shape GWAS summary statistics.  There are 3 main arguments: 1) _directory.working_ is simply your current working directory; 2) _filename_ is the name of the GWAS summary statistics file to process; 3) _BIG_ is the maximum number of GBs that will be __fully__ loaded in and worked with (default argument is 2).  If the size of the GWAS summary statistics file is larger than _BIG_, _APSS.R_ will perform an exploratory read of the data first.  This will hopefully shorten runtime and handle GWAS summary statistics files larger than 10 GBs.

The code provided below to create DNAm prediction models relies on specific, processed data.  We cannot provide all of the data due to privacy concerns, but we can provide a list of the used datasets.  Please reach out if you need assistance preparing data.


## Step 3: Conduct MWAS



The main code used to train DNAm prediction models is _mainbody_cpp_rsid_precise_justCreateWeightsFinal.R_.  You'll only need to specify the _name\_batch_ argument, which is the desired output name.  This code produces a set of 900 potential imputation models at each CpG site, of which we select the best in step 3.  

### Parallelizability

_mainbody_cpp_rsid_precise_justCreateWeightsFinal.R_ is parallel-friendly.  There is a small bit of code that guaranteees mutual exclusion for all subjobs, meaning you can run the code on as many parallel instances as needed and it will automatically detect unfinished jobs.  Note that it is currently set up to run on a slurm-based computing cluster, though this can be changed with minimal effort.  Simply edit lines 4-5 where the _job.id_ variable is defined.

### Alignment

| Data | Usage | Where to Find |
| --- | --- | --- |
| git status | List all new or modified files | fdfs |

_mainbody_cpp_rsid_precise_justCreateWeightsFinal.R_ assumes that you want to align mQTL files by rsID, so you'll need to have rsIDs in your mQTL files.  




