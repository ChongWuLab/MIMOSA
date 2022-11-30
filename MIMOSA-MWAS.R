array  <- Sys.getenv("SLURM_ARRAY_TASK_ID")
id.job <- as.numeric(array)

# Packages
library(BEDMatrix)
suppressMessages(library(data.table))
suppressMessages(library(ddpcr))
suppressMessages(library(dplyr))
library(optparse)

# Options
args=(commandArgs(TRUE))
for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
}
cat(path.ref," ",trait," ",path.out," ", path.trait, path.weight, "\n")


# User defined functions

# PatchUp
PatchUp <- function(M) {
    M <- apply(M, 2, function(x) {
        x[is.na(x)] <- mean(x, na.rm = TRUE)
        return(x)
    })

    return(M)
}


out <- data.frame(
    CpG         = character(),
    chromosome  = numeric(),
    r2_test     = numeric(),
    p           = numeric(),
    z           = numeric(),
    runtime     = numeric(),
    pos         = numeric(),
    stringsAsFactors = FALSE
)

files <- dir(path.weight)


# We have DNAm prediction models for 123,430 CpG sites
# Split into 412 pieces, so 300 arrays
for (cpg.index in (1 + (id.job - 1) * 412):(id.job * 412)) {

    if ( cpg.index > length(files)){break}
    cat(cpg.index, "\n")

    # Start tracking runtime
    time.start <- proc.time()[3]

    # The vector to store all the updates in this iteration
    update <- rep(NA, )

    # TODO: 3. Use your directory 
    # Load weight
    weightList <- readRDS(paste0(path.weight, files[cpg.index]))
    snps <- weightList[[1]]
    out.weight <- weightList[[2]]
    
    cpg.file <- files[cpg.index]
    cpg.name <- gsub(".rds","",cpg.file)
    

    update[3] <- weightList[[3]]
    if (update[3] <= 0.005) {
       cat("maxR2 <= 0.005 in test data.  Skipped", "\n \n")
       next
    }

    # Get the chromosome
    chr <- snps$SNPChr[1]


    # Read the summary statistics file
    if (file.exists(paste0(path.trait , trait, "-", chr, ".sumstats"))) {
       ss <- paste0(path.trait , trait, "-", chr, ".sumstats") %>% fread() %>% as.data.frame()
    } else {
      cat("No GWAS Summary Statistics file found.  Skipped", "\n \n")	   
      next
    }
    names(ss) <- colnames(ss) %>% tolower()
    ss["ID"] <- paste0(ss$chr, "_", ss$pos)

    # Load the reference panel
    quiet(
        bim.ref <- as.data.frame(fread(paste0(path.ref, chr, ".bim")))
    )
    quiet(
        genotype.ref <- BEDMatrix(paste0(path.ref, chr), simple_names = TRUE)
    )

    names(bim.ref)[c(2, 5, 6)] <- c("SNP" ,"a1", "a2")

    
    # Create new identifier for reference panel
    bim.ref["ID"] <- paste0(bim.ref$V1, "_", bim.ref$V4)

    # Find the common snps in all three data sets
    list.common <- intersect(bim.ref$SNP, snps$SNP) %>% intersect(., ss$snp) ###just match on rsID
    
    # Skip if no common snps found
    if (length(list.common) == 0) {
       cat("No common SNPs found.  Skipped.", "\n \n")
       next
    }

    # Trim genotype.ref
    genotype.temp <- genotype.ref[, bim.ref$SNP %in% list.common]
    genotype.temp <- as.matrix(genotype.temp)
    bim.temp      <- bim.ref[bim.ref$SNP %in% list.common, ]

    # Fix the NAs in reference panel
    if (sum(is.na(genotype.temp)) != 0) {
        genotype.temp <- PatchUp(genotype.temp)
    }

    # Trim ss
    ss.temp <- ss[ss$snp %in% list.common, ]

    # Trim snps and wgt.matrix
    index.temp <- snps$SNP %in% list.common

    snps       <- snps[index.temp, ]
    out.weight <- out.weight[index.temp,] 
    
    rm(index.temp)

    # Re-arrange datasets
    m.1 <- match(ss.temp$snp, bim.temp$SNP)
    m.2 <- match(ss.temp$snp, snps$SNP)

    bim.temp      <- bim.temp[m.1, ]
    genotype.temp <- genotype.temp[, m.1]
    genotype.temp <- as.matrix(genotype.temp)

    snps       <- snps[m.2, ]
    out.weight <- out.weight[m.2]

    # Align the mismatched alleles
    problem.1 <- ss.temp$a1 != bim.temp$a1
    genotype.temp[, problem.1] <- 2 - genotype.temp[, problem.1]
    genotype.temp <- as.matrix(genotype.temp)

    problem.2 <- ss.temp$a1 != snps$a1

    # Compute LD matrix
    genotype.temp <- scale(genotype.temp)
    genotype.temp <- as.matrix(genotype.temp)
    matrix.LD  <- t(genotype.temp) %*% genotype.temp / (nrow(genotype.temp) - 1)

    # Catch when there is only one row in wgt.matrix
    if ("numeric" %in% class(out.weight)) {
        out.weight2 <- out.weight %>% as.matrix() %>% t() %>% as.data.frame()
    }

    # Catch any NA weights
    if (sum(is.na(out.weight) > 0)){
       index.NA <- is.na(out.weight)
       out.weight[index.NA] <- 0
       rm(index.NA)
    }

    # Catch potential over-fitting
    if (max(abs(out.weight)) >= 10){
       out.weight <- rep(0, length(out.weight))
    }   

    # Association test 
    weights <- out.weight

    # Skip if weight is a zero vector
    if (sum(weights) == 0) {
        update[5] <- NA
	cat("All weights in list.common are 0.  Skipped.", "\n, \n")
	next
    }

    # Keep the non-zero components of weights vector
    keep <- (weights != 0)
    #print(sum(keep)
    cat(sum(keep), " nonzero weights.", "\n")
    weights <- weights[keep]
    
    # Compute MWAS z-score, r2, and p-value
    z.mwas  <- as.numeric(weights %*% ss.temp$z[keep])
    r2.mwas <- as.numeric(weights %*% matrix.LD[keep, keep] %*% weights)

    update[4]  <- 2 * (pnorm(abs(z.mwas / sqrt(r2.mwas)), lower.tail = F))
    update[5] <- z.mwas
    

    #################
    # Output format #
    #################
    # 1.     CpG site name
    # 2.     Chromosome
    # 3.     R^2 on testing data 
    # 4.     P-value of MWAS
    # 5.     Z-score of MWAS
    # 6.     Runtime
    # 7.     CpG site position
    
    # Update
    update[1]  <- snps$CpG[1]
    update[2]  <- chr
    
    # Update 7
    update[7] <- snps$CpGPos[1]

    # Stop tracking runtime
    time.end <- proc.time()[3]
    update[6] <- time.end - time.start

    out[nrow(out) + 1, ] <- update
    cat(cpg.index, " successful.", "\n \n")
}

# Write the result
dir.create(paste0(path.out, trait))
write.table(out, file = paste0(path.out, trait, "/", trait, "-", id.job, ".txt"), row.names = FALSE, quote = FALSE)


# Preparing your summary statistics
# Please follow my naming patterns!
# For instance, the raw data file for Asthma was named "AST" and was split into 22 smaller files.
# The smaller files were named AST-1.sumstats, ..., AST-22.sumstats.
# For each smaller files, make sure that your header is ALL CAPS (to avoid potential problems, please follow the way that I organized my ss files).
# You can manually process your raw ss files. (Read this useful blog: https://huwenboshi.github.io/data%20management/2017/11/23/tips-for-formatting-gwas-summary-stats.html)
# To save time, you can try the tool "APSS".

##########################
# An example run of APSS #
##########################
# source("APSS.R")
# APSS(
#    directory.working = "./"
#    filename = "LDL",
#    auto = FALSE,
#    do.return = FALSE,
#    BIG = 3
#    )
