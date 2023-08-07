slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
job.id <- as.numeric(slurm_arrayid)

# Packages
library(BEDMatrix)
suppressMessages(library(data.table))
suppressMessages(library(ddpcr))
suppressMessages(library(dplyr))
library(optparse)


#path.ref <- opt$path.ref
#trait    <- opt$trait
#path.out <- opt$path.out

args=(commandArgs(TRUE))
for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
}
cat(path.ref," ",trait," ",path.out," ", path.trait, " ", path.weight,"\n")

#examples
#path.ref='resources/1000-genomes/1000G.EUR.ALLSNP.QC.CHR'
#trait='HCH'
#path.trait='GWAS-summary-statistics/'
#path.out='Results/'
#path.weight='MIMOSA-models/

# ACAT
source("ACAT.R")

# PatchUp
source("PatchUp.R")

# Object to store the output
out <- data.frame(
    runtime     = character(),
    CpG         = character(),
    chromosome  = numeric(),
    model_best  = character(),
    r2_best     = numeric(),
    p_ElNet     = numeric(),
    p_MNet      = numeric(),
    p_SCAD      = numeric(),
    p_MCP       = numeric(),
    p_LASSO     = numeric(),
    z_ElNet     = numeric(),
    z_MNet      = numeric(),
    z_SCAD      = numeric(),
    z_MCP       = numeric(),
    z_LASSO     = numeric(),
    p_Union     = numeric(),
    p_ACAT      = numeric(),
    CpG_pos    = numeric(),
    stringsAsFactors = FALSE
)

# Main iteration

files <- dir(path.weight) #123365  300 arrays of 415 CpG sites
cpg.names <- gsub(".rds","", files)


for (cpg.index in  (1+(415*(job.id-1))):(min(job.id*415, 123365))) {
    # Start tracking runtime
    time.start <- proc.time()[3]

    # The vector to store all the updates in this iteration
    update <- rep(NA, 18)

    # Load weight
    w.check <- try(readRDS(paste0("../WeightsMWASPrepped/MWASPreppedAllrmFAVORnew/", files[cpg.index])))
    if("try-error" %in% class(w.check)){next}
    weight.list <- readRDS(paste0(path.weight, files[cpg.index]))
    names(weight.list) <- c("ElNet", "MNet", "SCAD", "MCP", "LASSO")
    #the weight lists have five sub-lists, one for each of ElNet, MNet, SCAD, MCP, LASSO
    #each sub-list has five elements: T/F for if the weight exists for that method, mQTL summary stats, actual weights, R^2 on test data, lambda

    R2sVec <- c()
    for (i in 1:5){
    	if(!(length(weight.list[[i]][[4]]) == 0)){  #catch if R^2 test is empty but not NA
       	   R2sVec <- append(R2sVec, weight.list[[i]][[4]])
    	} else {
    	   weight.list[[i]][[1]] <- FALSE
    	   R2sVec <- append(R2sVec, NA)
    	}
     
    }
    
    update[5] <- max(R2sVec[1:5], na.rm = TRUE)
    if (update[5] <= 0.005) {
        next
    }

    # Get the chromosome
    snpsTrue <- c()
    for (i in 1:5){
      if (weight.list[[i]][[1]]){
        snpsTrue <- append(snpsTrue, i)
      }
    }
    snps.index <- min(snpsTrue)
    snps <- weight.list[[snps.index]][[2]]
    chr <- snps$SNPChr[1]

    # Read the summary statistics file
    if (file.exists(paste0(path.trait, trait, "-", chr, ".sumstats"))) {
        ss <- paste0(path.trait, trait, "-", chr, ".sumstats") %>% fread() %>% as.data.frame()
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

    # Find the best model
    bestIdx <- which(R2sVec == update[5])
    bestIdx <- bestIdx[1]
    model.best <- names(weight.list)[bestIdx]
    
    # Create new identifier for reference panel
    bim.ref["ID"] <- paste0(bim.ref$V1, "_", bim.ref$V4)


    
    # Find the common snps in all three data sets
    # Iterate by method
    for (i in 1:5){
    	gc()
  	if (!weight.list[[i]][[1]]){
	     next
	  }   
	    snps <- weight.list[[i]][[2]]
    	list.common <- intersect(bim.ref$SNP, snps$SNP) %>% intersect(., ss$snp)

        # Skip if no (or one) common snps found
        if (length(list.common) <= 1) {
            next
        }

        # Trim genotype.ref
        genotype.temp <- genotype.ref[, bim.ref$SNP %in% list.common]
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
	      out.weight <- weight.list[[i]][[3]][index.temp]
	      if ("numeric" %in% class(out.weight)) {
           out.weight <- t(as.matrix(out.weight))
      	}

        rm(index.temp)

        # Re-arrange every data sets
        m.1 <- match(ss.temp$snp, bim.temp$SNP)
        m.2 <- match(ss.temp$snp, snps$SNP)

        bim.temp      <- bim.temp[m.1, ]
    	  genotype.temp <- genotype.temp[, m.1]

      	snps       <- snps[m.2, ]
      	out.weight <- out.weight[,m.2 ]

      	# Align the mismatched alleles - ref panel (bim.temp) and snps (i.e. weights) are already aligned
       	problem <- ss.temp$a1 != snps$a1
      	ss.temp$a1 <- snps$a1
      	ss.temp$a2 <- snps$a2
      	ss.temp$beta[problem] <- -1 * ss.temp$beta[problem]
      	ss.temp$z[problem] <- -1 * ss.temp$z[problem]
    	
      
      
    	
        # Compute LD matrix
        genotype.temp <- scale(genotype.temp)
      	matrix.LD  <- t(genotype.temp) %*% genotype.temp / (nrow(genotype.temp) - 1)

      # Catch: When there is only one row in wgt.matrix
    	if ("numeric" %in% class(out.weight)) {
           out.weight <- out.weight %>% as.matrix() %>% t() %>% as.data.frame()
    	}

	# Catch: Over-fitting and NAs that may have snuck through
      	for (j in 1:ncol(out.weight)){
      	  if (is.na(out.weight[1,j])){
      	    out.weight[1,j] <- 0
      	  }
      	}
        if (max(abs(out.weight)) >= 10) {
            out.weight <- rep(0, ncol(out.weight))
        }

        # Settings
        weights <- out.weight

        # Skip if weight is a zero vector
        if (sum(weights) == 0) {
            update[5 + i] <- NA  #if weights add to 0, set pval to NA
            next
        }

        # Keep the non-zero components of weights vector
        keep <- (weights != 0)
        weights <- weights[keep]

        # Compute TWAS z-score, r2, and p-value
        z.mwas  <- as.numeric(weights %*% ss.temp$z[keep])
        r2.mwas <- as.numeric(weights %*% matrix.LD[keep, keep] %*% weights)

        update[5 + i]  <- as.numeric(2 * (pnorm(abs(z.mwas / sqrt(r2.mwas)), lower.tail = F)))
        update[10 + i] <- as.numeric(z.mwas)
    }
    
    # Union of all models
    update[16] <- update[bestIdx + 5]

    # ACAT on all models
    check.na   <- !is.na(update[6:10])
    check.sign <- R2sVec > 0
    check.final	<- check.na & check.sign
    if (sum(check.final) == 0) {
	     update[17] <- NA
    } else {
       update[17] <- ACAT(update[6:10][check.final], R2sVec[check.final] / sum(R2sVec[check.final]))
    }


    # Stop tracking runtime
    time.end <- proc.time()[3]

    #################
    # Output format #
    #################
    # 1.     Runtime
    # 2.     CpG Site
    # 3.     Chromosome
    # 4.     Best model
    # 5.     Best model's R^2 on testing data
    # 6-10.  P-value of MWAS from weights constructed by (ElNet, MNet, SCAD, MCP, LASSO)
    # 11-15. Z-score of MWAS from weights constructed by (ElNet, MNet, SCAD, MCP, LASSO)
    # 16.    P-value of MWAS from best performing model (in terms of R^2 on testing data)
    # 17.    P-value from ACAT on all models
    # 18.    CpG site position
    
    # Update
    update[1]  <- time.end - time.start
    update[2]  <- snps$CpG[1]
    update[3]  <- chr
    update[4]  <- model.best
    
    annotation <- readRDS("IlluminaHumanMethylation450kanno.rds")
    update[18] <- annotation[update[2],'pos']
    

    out[nrow(out) + 1, ] <- update
    cat(cpg.index, "\n")
}

# Write the result
dir.create(paste0(path.out, trait))
write.table(out, file = paste0(path.out, trait, "/",trait, "-", job.id), row.names = FALSE, quote = FALSE)
