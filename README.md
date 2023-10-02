# MIMOSA

__MIMOSA__, or __MWAS Imputing Methylome Obliging Summary-level mQTLs and Associated LD matrices__, is a set of models that substantially improve the prediction accuracy of DNA methylation and downstream MWAS power through the use of a large, summary-level mQTL dataset provided by the Genetics of DNA Methylation Consortium (GoDMC).  Here, we provide a tutorial to download the models and run MWAS with them.  

### Citation

If you use the MIMOSA models, please cite

> Melton, H. J., Zhang, Z., Deng, H., Wu, L., & Wu, C. (2023). MIMOSA: A resource consisting of improved methylome imputation
models increases power to identify CpG site-phenotype associations. 

## Step 1: Download DNAm prediction models

The MIMOSA models for DNAm prediction in whole blood are available on zenodo.org [here](https://zenodo.org/record/7325055).  The files included are _cg********.rds_, which each consists of a set of five weights for SNPs used in predicting DNAm at the CpG site cg*********.  Once you load one of these files into R, you'll have a list of five lists, with each corresponding to each of the penalized regression methods used to build the models (ElNet, MNet, SCAD, MCP, and LASSO).  Each of these lists will contain elements 1) TRUE/FALSE if this method produced a satisfactory DNAm prediction model; 2) a dataframe with (among other things) mQTL p-value, rsID, chromosome, position (of SNP), a1, a2, and CpG site; 2) a set of DNAm prediction weights corresponding to the SNPs in 1); 3) the prediction accuracy (R^2) of the model for the CpG site in the test dataset (which comes from the Framingham Heart Study).  Once you have the weights, you're halfway to running your MWAS.

## Step 2: Prepare GWAS Summary Statistics

The GWAS summary statistics you use for your MWAS will need to be in the correct format to work with MIMOSA.  Please follow these conventions:

* Split the raw data into 22 subfiles, one for each chromosome.  If the trait you're considering is asthma (AST), then name these subfiles AST-1.sumstats, ..., AST-22.sumstats.
* Ensure the header of the subfiles are in ALL CAPS.  At minimum, you should have columns for SNP (where this is rsIDs), CHR, A1, A2, Z.
* You can manually process the summary statistics, or you can use the provided APSS.R function

### APSS.R

_APSS.R_ is an interactive R function that helps you easily process and shape GWAS summary statistics.  There are 3 main arguments: 1) _directory.working_ is simply your current working directory. 2) _filename_ is the name of the GWAS summary statistics file to process. 3) _BIG_ is the maximum number of GBs that will be __fully__ loaded in and worked with (default argument is 2).  If the size of the GWAS summary statistics file is larger than _BIG_, _APSS.R_ will perform an exploratory read of the data first.  This will hopefully shorten runtime and handle GWAS summary statistics files larger than 10 GBs.


## Step 3: Conduct MWAS

The script used to run MWAS is _MIMOSA-MWAS.R_.  There are five arguments for this script: 1) _path.ref_ is the pathway to your LD reference panel (we used the 1000 Genomes Project and downloaded from [here](https://www.internationalgenome.org/data)).  2) _trait_ is the name of the trait of interest.  This argument should be the same as the name of your GWAS summary statistics files (i.e. AST, using the example from earlier).  3) _path.trait_ is the pathway to your GWAS summary statistics.  4) _path.out_ is the pathway to wherever you'd like to save your results.  5) _path.weight_ is the pathway to the directory where you stored the MIMOSA models.

_MIMOSA-MWAS.R_ relies on a few other key scripts that should be in your working directory: _ACAT.R_ and _PatchUp.R_.  These are available to download from our Github page.  _MIMOSA-MWAS.R_ also requires the file _IlluminaHumanMethylation450kanno.rds_ (downloaded from Illumina's website) in your working directory to match positions to the CpG sites.

Note that _MIMOSA-MWAS.R_ is currently set up to run across 300 parallel instances on a slurm cluster.  If you need to change this, you'll need to edit lines 1-2 that determine which of the 300 instances is running.  You'll probably want to still use the variable _id.job_ to store which instance is running, since it is referenced later in the code.  You would also need to edit the for loop on line 64 depending on how many jobs you're running.

A slurm submission script, _MIMOSA-MWAS.sh_, is provided for ease of use.

### Output Format

| Column | Name | Description |
| --- | --- | --- |
| 1 | runtime | Runtime |
| 2 | CpG | Name of CpG site |
| 3 | chromosome | Chromosome |
| 4 | model_best | Whichever of the 5 penalized regression models performed best at DNAm prediction on the test data |
| 5 | r2_best | R^2 for best DNAm prediction model on test data |
| 6-10 | p_Method | p-value for MWAS for each method |
| 11-15 | z_Method | z-score for MWAS for each method |
| 16 | p_Union | p-value of MWAS from best performing model (in terms of R^2 on testing data) |
| 17 | p_ACAT | p-value from ACAT on all satistfactory models |
| 18 | CpG_pos | Position of CpG site |

## Disclaimer

The DNAm prediction models and software are provided “as is”, without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose and noninfringement. in no event shall the authors or copyright holders be liable for any claim, damages or other liability, whether in an action of contract, tort or otherwise, arising from, out of or in connection with the models or the use or other dealings in the models.

## License

Maintainer: Hunter Melton (hjm19d@fsu.edu)

[MIT](http://opensource.org/licenses/MIT)

Copyright (c) 2013-present, Hunter Melton (hjm19d@fsu.edu), Zichen Zhang (zz17@fsu.edu), Chong Wu (cwu18@mdanderson.org)
