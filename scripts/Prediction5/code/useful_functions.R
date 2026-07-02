

#' Get cross-validation folds
#'
#' @param index length of total set to select from, numeric vector of length 1
#' @param nb.folds number of fold to divide, default is 5
#' @param seed (optional) for reproducibility, set a seed.
#'
#' @return list of testing set indexes for each fold.

getFolds <- function(index, nb.folds=5, seed=NULL){
  stopifnot(nb.folds <= index,
            ifelse(! is.null(seed), is.numeric(seed), TRUE))
  if(!is.null(seed)){set.seed(seed, kind="Mersenne-Twister")}
  split <- split(sample(1:index), rep(1:nb.folds, length=index))
  names(split) <- paste0("Fold", seq(1:nb.folds))
  return(split)
}


#' Format and curate vcf file
#'
#' @param vcf.p2f path to the vcf file
#' @param matrix.gt matrix of genotypes, must contains CHROM and POS columns, default is NULL
#' @param p2f.export.vcf path to export curated vcf file, default is NULL (no export), useful for further Beagle imputation.
#' @param IDnum logical, remove the prefix idxx. from the genotype name, default is FALSE
#' @param mrk.info data frame with marker information pulled from the vcf, with columns CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, default is NULL. Needed if `vcf.p2f` is not provided and `p2f.export.vcf` is not null.
#' @param corresp.geno.name table of correspondence between genotype names, with the first column being the name to be matched and the second one being the updated name.
#' @param remove.chrUkn logical, remove markers from unknown chromosome, default is TRUE
#' @param check.mrk.dups logical, check for duplicated markers, default is TRUE
#' @param check.geno.dups logical, check for duplicated genotypes, default is TRUE
#' @param remove.nonPolyMrk logical, remove non-segregating markers, default is TRUE
#' @param tresh.heterozygous numeric, threshold of the maximum proportion of heterozygous genotypes, default is NULL
#' @param thresh.NA.ind numeric, threshold of the maximum proportion missing values for individuals, default is 0.5
#' @param thresh.NA.mrk numeric, threshold of the maximum proportion of missing values for markers, default is 0.5
#' @param format.names logical, use a custom function to transform genotype names, default is FALSE
#' @param concordance.function logical, use a custom function to match genotype names, default is FALSE
#' @param corresp.geno.list named list of correspondence between genotype names, with correct spelling as name and possible synonym as character vector, default is NULL.
#' @param imputation character, impute missing values with kNNI or Beagle or don't impute, default is NULL (no imputation)
#' @param p2f.beagle path to export file for Beagle imputation, default is NULL
#' @param thresh.MAF numeric, threshold of minor allele frequency, default is NULL
#' @param verbose numeric, level of verbosity, default is 1
#'
#' @return data frame with 0/1/2 values from the curated vcf file, with genotypes in row and markers in columns.
#' @author Prediction5 Team
#' @seealso [format_names()], [getGenoTas_to_DF()], [concordance_match_name()]
#' @importFrom vcfR read.vcfR extract.gt getFIX
#' @importFrom methods new
#' @importFrom dplyr arrange filter mutate select distinct
#' @importFrom tibble as_tibble

#' @export

format_curate_vcf <- function(vcf.p2f=NULL,
                              matrix.gt=NULL,
                              mrk.info=NULL,
                              corresp.geno.name=NULL,
                              p2f.export.vcf=NULL,
                              IDnum=FALSE,
                              remove.chrUkn=TRUE,
                              check.mrk.dups=TRUE,
                              remove.nonPolyMrk=TRUE,
                              tresh.heterozygous=NULL,
                              check.geno.dups=TRUE,
                              thresh.NA.ind=0.5,
                              thresh.NA.mrk=0.5,
                              format.names=FALSE,
                              concordance.function=FALSE,
                              corresp.geno.list=NULL,
                              imputation=NULL,
                              p2f.beagle=NULL,
                              thresh.MAF=NULL,
                              verbose=1){

  ### Verifications
  stopifnot(thresh.NA.mrk >= 0, thresh.NA.mrk <= 1,
            thresh.NA.ind >=0, thresh.NA.ind <= 1,
            is.logical(check.geno.dups), is.logical(check.mrk.dups),
            is.logical(remove.chrUkn), is.logical(format.names),
            is.logical(IDnum),is.logical(remove.nonPolyMrk),
            ## check if vcf.p2f or matrix.gt is provided
            xor(is.null(vcf.p2f), is.null(matrix.gt)),
            imputation %in% c(NULL,"kNNI","Beagle"),
            length(imputation) <=1)
  #vcf.file <- as.data.frame(vcf.file)


  ### 1. Load VCF file
  if(!is.null(vcf.p2f)){
    stopifnot(file.exists(vcf.p2f))

    if(verbose > 0) print("Load the vcf file...")
    ## load vcf file and convert to gt
    vcf <- vcfR::read.vcfR(vcf.p2f, verbose = FALSE)
    gt <- as.data.frame(vcfR::extract.gt(vcf, element = "GT", as.numeric=FALSE, IDtoRowNames = FALSE))
    meta <- vcf@meta
    mrk.info <- tibble::as_tibble(vcfR::getFIX(vcf))
    gt$ID <- mrk.info$ID
    colsGT <- colnames(gt)[!colnames(gt) %in% c("ID")]
    mrk.info$POS <- as.numeric(mrk.info$POS)
    mrk.info <- dplyr::arrange(mrk.info, CHROM,POS)
    ## combine marker information with genotype
    vcf.file <- cbind(mrk.info[,c("CHROM","POS")],gt[match(mrk.info$ID,gt$ID),colsGT])
    vcf.file$POS <- as.numeric(vcf.file$POS)

    rm(vcf, gt)

  } else {
    stopifnot(!is.null(matrix.gt), all(c("CHROM","POS") %in% colnames(matrix.gt)[c(1,2)]))
    vcf.file <- matrix.gt
    meta <- NULL
  }

  if(verbose >0){
    print(paste0("Initial dimensions: ", nrow(vcf.file)," markers and ",
                 (ncol(vcf.file) - 2), " genotypes"))
  }

  ## remove the prefix idxx. from the genotype name
  if(IDnum){
    colnames(vcf.file) <- gsub(x=colnames(vcf.file),pattern="^id[0-9]+\\.",replacement="")
  }
  ## replace patterns "./." by NA
  vcf.file[,3:ncol(vcf.file)] <- as.data.frame(apply(vcf.file[,3:ncol(vcf.file)],2, function(x)
    gsub(pattern="\\.",x=x,replacement=NA)))
  ### replace multi-allelic marker by mono-allelic, i.e., "2/2"" by "1/1"
  # vcf.file[,3:ncol(vcf.file)] <- as.data.frame(apply(vcf.file[,3:ncol(vcf.file)],2, function(x)
  #   gsub(pattern="2/2",x=x,replacement="1/1")))
  vcf.file[,3:ncol(vcf.file)] <- as.data.frame(apply(vcf.file[,3:ncol(vcf.file)],
                                                     2, function(x){
                                                       x=gsub(pattern="2/2",x=x,replacement="1/1")
                                                       x=gsub(pattern="0/2",x=x,replacement="0/1")
                                                       x=gsub(pattern="2/0",x=x,replacement="1/0")
                                                     }))
  ## set position as numeric
  vcf.file$POS <- as.numeric(stringr::str_trim(vcf.file$POS))

  ## update genotype names
  if(!is.null(corresp.geno.name)){
    stopifnot(ncol(corresp.geno.name) >=2)
    if(verbose >0){
      print("Update genotype names..")
    }

    colnames(vcf.file) <- plyr::mapvalues(colnames(vcf.file),
                                          from=corresp.geno.name[,1],
                                          to=corresp.geno.name[,2], warn_missing=FALSE)
  }

  ## handle unknown chromosome, remove markers if few
  if(remove.chrUkn){
    chr.unk <- c("CHRUnknown","U","ChrUnknown","Unknown","Un","ChrUn","chrU","ChrU")
    if(verbose >0){
      print(paste0("Removed ",nrow(vcf.file[vcf.file$CHROM %in% chr.unk,]),
                   " markers from unknown chromosome"))
    }
    vcf.file <- vcf.file[!vcf.file$CHROM %in% chr.unk,]
  }

  ## check for duplicated markers
  # extract map information
  vcf.file$chr_pos <- paste0(vcf.file$CHROM, "_",vcf.file$POS)
  if(check.mrk.dups){
    dup_snp <- unique(vcf.file$chr_pos[duplicated(vcf.file$chr_pos)])
    if(verbose >0){
      print(paste0("Found ",length(dup_snp)," duplicated markers"))
    }
    ## for duplicated loci, keep the one with the lowest missing value
    idx2rem <- c()
    for(mrk in dup_snp){
      ## get row id for duplicated markers
      idx.mrkdup <- grep(mrk,vcf.file$chr_pos)
      ## detect rows with less missing values
      idx2rem <- c(idx2rem,idx.mrkdup[-which.min(rowSums(is.na(vcf.file[vcf.file$chr_pos == mrk,])))])
    }
    ## remove duplicated columns with more missing values
    if(length(idx2rem)>0){
      vcf.file <- vcf.file[-idx2rem,]
    }
  }
  rownames(vcf.file) <- vcf.file$chr_pos
  ## if there is no chr before chr_pos, add it
  if(all(grepl("[1-7]|U",substr(rownames(vcf.file),1,1)))){
    rownames(vcf.file) <- paste0("chr",rownames(vcf.file))
  }
  cols2rem <- c("CHROM","POS", "chr_pos")
  vcf.file[,cols2rem]<- NULL

  ## 2. remove duplicated genotypes

  ## format genotype name with custom function format_names
  if(format.names){
    colnames(vcf.file) <- format_names(colnames(vcf.file))$corrected
  } else if(any(grep(" ", colnames(vcf.file))) & !is.null(p2f.export.vcf)){
    ## need to remove white spaces for synbreed export, replace it by underscore
    colnames(vcf.file) <- gsub(x=colnames(vcf.file),pattern=" ",replacement = "_")
  }
  if(concordance.function){
    stopifnot(!is.null(corresp.geno.list))
    conc <- concordance_match_name(match.name=colnames(vcf.file), corresp.list=corresp.geno.list)
    colnames(vcf.file) <- conc$new.names
  }

  if(check.geno.dups){
    ## classic detection of duplicated genotypes
    dup.genos1 <- unique(colnames(vcf.file)[duplicated(colnames(vcf.file))])

    ### detect duplicated genotypes with .1/.2/... appended at the end of the name
    dup.genos2 <- grep(pattern="\\.[1-9]{1,2}$", colnames(vcf.file), value=TRUE)
    dup.genos2 <- unique(stringr::str_remove(dup.genos2, pattern="\\.[1-9]{1,2}$"))
    dup.genos <- unique(c(dup.genos1,dup.genos2))
    if(verbose >0){
      print(paste0("Found ",length(dup.genos)," duplicated genotypes"))
    }
    ## for duplicated genotypes, fill if missing data and keep the one with less missing data
    for(gen in dup.genos){
      idx.genodup <- grep(paste0("^",gen,"[.1-9]{0,2}"),colnames(vcf.file))
      #colnames(vcf.file.sel)[idx.genodup]
      if(length(idx.genodup) > 1){
        tmp <- vcf.file[,idx.genodup]

        ## find the column with the least NA values
        idxLessNa <- which.min(colSums(is.na(tmp)))
        ## fill NA values with data from other columns if available
        idxna <- which(is.na(tmp[,idxLessNa]))
        for(i in 1:ncol(tmp)){
          tmp[idxna,idxLessNa] <- dplyr::coalesce(tmp[idxna,idxLessNa],tmp[idxna,i])
        }
        ## suppress duplicated columns, keep the one with less NA values and filled
        vcf.file[,idx.genodup[idxLessNa]] <- tmp[,idxLessNa]
        colnames(vcf.file)[idx.genodup[idxLessNa]] <- gen
        idx2rem <- setdiff(idx.genodup,idx.genodup[idxLessNa])
        vcf.file = vcf.file[,-idx2rem]
      }
    }
    stopifnot(length(unique(colnames(vcf.file))) == ncol(vcf.file))
  }


  ## 3. Curate: remove SNPs/ genos with too many missing data

  ### remove markers with too many NA
  snp.miss <- apply(vcf.file, 1, function(x) sum(is.na(x)))
  snp2rem <- names(snp.miss)[snp.miss > thresh.NA.mrk*ncol(vcf.file)]
  if(verbose > 0){
    print(paste0("Removed ",length(snp2rem)," markers with too many missing data"))
  }
  vcf.file <- vcf.file[!rownames(vcf.file) %in% snp2rem,]
  if(nrow(vcf.file) == 0){
    stop("No markers left after removing markers with too many missing data")
  }

  ### remove individuals with too many NA
  gen.miss <- apply(vcf.file, 2, function(x) sum(is.na(x)))
  geno2rem <- names(gen.miss)[gen.miss > thresh.NA.ind*nrow(vcf.file)]
  if(verbose > 0){
    print(paste0("Removed ",length(geno2rem)," individuals with too many missing data"))
  }
  vcf.file <- vcf.file[,!colnames(vcf.file) %in% geno2rem]
  if(ncol(vcf.file) < 4){
    stop("No genotypes left after removing genotypes with too many missing data")
  }

  ## Remove non-polymorphic markers
  if(remove.nonPolyMrk){
    nb.allele <- apply(vcf.file, 1, function(x) length(table(t(x))))
    mrk2rem <- names(nb.allele)[nb.allele <= 1]
    if(verbose > 0){
      print(paste0("Removed ",length(mrk2rem)," non-polymorphic markers"))
    }
    if(length(mrk2rem) >0) vcf.file <- vcf.file[!rownames(vcf.file) %in% mrk2rem,]
    rm(mrk2rem)
  }

  ## Remove markers with too many heterozygous genotypes
  if(!is.null(tresh.heterozygous)){
    hetero <- apply(vcf.file, 2, function(x) {sum(x %in% c("0/1","1/0","0|1","1|0"), na.rm = TRUE)})
    geno2rem <- names(hetero)[hetero > tresh.heterozygous * nrow(vcf.file)]
    if(verbose > 0){
      print(paste0("Removed ",length(geno2rem)," genotypes with too many heterozygous markers"))
    }
    vcf.file <- vcf.file[,!colnames(vcf.file) %in% geno2rem,]
  }

  rownames(vcf.file) <- gsub(x = rownames(vcf.file),pattern="Chr",replacement = "chr")
  ## recode numerically: replace 0/0=0, 0/1 or 1/0=1 and 1/1=2
  rowN <- rownames(vcf.file)

  vcf.num <- vcf.file
  vcf.num <-
    apply(vcf.num, 2,
          function(x) {
            as.numeric(stringr::str_replace_all(gsub("[[:punct:]]", "", x),
                                                c("00"="0","01"="1","10"="1","11"="2")))
            })
  vcf.num <- as.data.frame(vcf.num)
  rownames(vcf.num) <- rowN

  ### 4. (optional) export to a formatted vcf file for Beagle imputation

  if(!is.null(p2f.export.vcf)){
    ## if no extension or vcf extension provided, extend to vcf.gz extension
    if(!grepl("vcf.gz$",p2f.export.vcf)){
      p2f.export.vcf <- paste0(gsub(".vcf","",p2f.export.vcf),".vcf.gz")
    }

    ## recreate vcfR object
    vcf2exp <- methods::new(Class = "vcfR")
    ### meta element
    if(!is.null(meta)){
      vcf2exp@meta <- meta
    } else {
      vcf2exp@meta <- c("##fileformat=VCFv4.2",
                        "##reference=xxx",
                        "##comment=Generated from format_curate_vcf function",
                        #"##filedate= 20065",
                        "##source='write.vcf of vcfR'",
                        "##FORMAT=<ID=GT,Number=1,Type=String,Description='Genotype'>")
    }

    ### fix element: add marker information
    #### define a default mrk.info data frame if not provided
    if(is.null(vcf.p2f) & is.null(mrk.info)){
      chr_pos <- rownames(vcf.num)
      mrk.info <- data.frame(CHROM=unlist(lapply(strsplit(chr_pos,"_"),function(x) x[1])),
                             POS=unlist(lapply(strsplit(chr_pos,"_"),function(x) x[2])),
                             ID=chr_pos,
                             REF="A",ALT="T",
                             QUAL=".",FILTER="PASS",INFO=".",
                             FORMAT="GT")
      mrk.info$CHROM <- gsub("Chr","chr",mrk.info$CHROM)

    } else {
      ## use the mrk.info from the vcf file, subset selected markers
      mrk.info$ID <- paste0(gsub("Chr","chr",mrk.info$CHROM),"_",mrk.info$POS)
      ## if marker name start with chromosome number, add chr
      if(all(grepl("^[1-9]|U",substr(mrk.info$ID,1,1)))){
        mrk.info$ID <- paste0("chr",mrk.info$ID)
      }
      mrk.info <- mrk.info[match(rownames(vcf.num), mrk.info$ID),]
    }
    if(!all(c("INFO","FORMAT") %in% colnames(mrk.info))){
      mrk.info <- cbind(mrk.info,
                        #QUAL=".",FILTER="PASS",
                        INFO=".",FORMAT="GT")
      mrk.info <- mrk.info[,c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")]
    }
    mrk.info <- mrk.info[match(rowN, mrk.info$ID),]
    stopifnot(nrow(mrk.info) == nrow(vcf.num))
    mrk.info <- as.matrix(mrk.info)
    mrk.info[,"POS"] <- gsub(" ","",mrk.info[,"POS"])
    vcf2exp@fix <- mrk.info

    ### gt element, coded in 0/0, 0/1, 1/1
    vcf2exp@gt <- as.matrix(vcf.file)

    ### export
    vcfR::write.vcf(x=vcf2exp, file=p2f.export.vcf, mask=FALSE, APPEND=FALSE)

  } ## end if !is.null(p2f.export.vcf)


  ### Imputation of missing markers with kNNI
  if(length(imputation) >0 && imputation == "kNNI"){
    if(is.null(p2f.export.vcf)){
      stop("Need to export vcf file to perform imputation")
    }
    ## load rTASSEL and rJava packages
    if(!requireNamespace("rTASSEL", quietly = TRUE)){
      stop("Please install rTASSEL package to perform kNNI imputation")
    }
    if(!requireNamespace("rJava", quietly = TRUE)){
      stop("Please install rJava package to perform kNNI imputation")
    }

    tasGenoVCF <- rTASSEL::readGenotypeTableFromPath(path = p2f.export.vcf)
    ## use KNNI as imputation method, using default parameters
    tasGenoImp <- rTASSEL::imputeLDKNNi(tasObj=tasGenoVCF, highLDSSites = 30,
                                        knnTaxa = 10, maxDistance = 1e+07)
    ## use an internal function to convert this object to a data frame
    vcf.file <- getGenoTas_to_DF(tasGenoImp)

    ## format it to export it as a vcf file
    vcf2exp@gt <- as.matrix(vcf.file)
    # gp$geno <- t(vcf.file[rownames(gp$map),rownames(gp$geno)])
    # stopifnot(identical(colnames(gp$geno),rownames(gp$map)))

    ## detect and warn if missing values left

    if(verbose > 0) {
      ### proportion of missing value per column
      print("Summary of missing data percentage per marker: ")
      print(summary(colMeans(is.na(gp$geno)) *100))

      ### proportion of missing value per line
      print("Summary of missing data percentage per genotype: ")
      print(summary(rowMeans(is.na(gp$geno)) *100))

    }
    # # export VCF
    # synbreed::write.vcf(gp,file=p2f.export.vcf)
    ### export
    vcfR::write.vcf(x=vcf2exp, file=p2f.export.vcf, mask=FALSE, APPEND=FALSE)

  } else if(length(imputation) >0 && imputation == "Beagle"){
    if(is.null(p2f.export.vcf)){
      stop("Need to export vcf file to perform imputation")
    }
    if(!file.exists(p2f.beagle)){
      stop("Beagle jar file not found")
    }
    ## run Beagle imputation
    p2f.imp <- gsub(".vcf.gz","_imputed",p2f.export.vcf)
    system(paste0("java -Xmx2g -jar ",p2f.beagle," gt=",p2f.export.vcf,
                  " out=",p2f.imp,
                  " nthreads=4 ne=100000"))

    ## load imputed vcf
    vcf.imp <- vcfR::read.vcfR(paste0(p2f.imp,".vcf.gz"), verbose = FALSE)
    gt.imp <- vcfR::extract.gt(vcf.imp, element = "GT", as.numeric=FALSE, IDtoRowNames = FALSE)
    mrk.info.imp <- tibble::as_tibble(vcfR::getFIX(vcf.imp))
    vcf.file <- cbind(mrk.info.imp[,c("CHROM","POS")],gt.imp)
    vcf.file$POS <- as.numeric(vcf.file$POS)
    vcf.file <- dplyr::arrange(vcf.file, CHROM,POS)
    rownames(vcf.file) <- paste0(vcf.file$CHROM, "_",vcf.file$POS)
    cols2rem <- c("CHROM","POS")
    vcf.file[,cols2rem]<- NULL
    rowN <- rownames(vcf.file)
    ## recode numerically: replace 0/0=0, 0/1 or 1/0=1 and 1/1=2
    vcf.num <- apply(vcf.file, 2, function(x) {
      as.numeric(stringr::str_replace_all(gsub("[[:punct:]]", "", x),
                                          c("00"="0","01"="1","10"="1","11"="2")))})
    vcf.num <- as.data.frame(vcf.num)
    rownames(vcf.num) <- rowN
  } ## end if Beagle imput


  ## Remove markers with MAF below threshold after imputation
  if(!is.null(thresh.MAF)){
    thresh.MAF <- as.numeric(thresh.MAF)
    if(thresh.MAF < 0 | thresh.MAF > 0.5){
      stop("MAF should be between 0 and 0.5")
    }

    maf <- apply(vcf.num, 2, function(x) {
      t = table(x) ; t = t[names(t) %in% c(0,2)]
      if(length(t) < 2) maf = 0
      else maf = sort(t)[1]/length(x)
      return(maf)
    })

    mrk2rem <- names(maf)[maf < thresh.MAF]

    if(verbose > 0){
      print(paste0("Removed ",length(mrk2rem)," markers with MAF below ",thresh.MAF))
    }
    vcf.num <- vcf.num[!rownames(vcf.num) %in% mrk2rem,]
  }


  if(verbose >0){
    print(paste0("Final dimensions: ",ncol(vcf.num), " genotypes and ", nrow(vcf.num)," markers"))
  }
  return(t(vcf.num)) ## output in dimension geno x marker
}



## from R metan package
# Function to make HTML tables
#' Print nice table using DT
#'
#' @param table data frame to print
#' @param rownames logical, if TRUE, print row names, default is FALSE
#' @param digits integer, number of digits to print for numeric columns, default is 3
#' @param ... additional arguments to pass to `DT::datatable`
#'
#' @returns a DT datatable object
#' @export
#' @examples
#' print_table(mtcars)
#' @importFrom DT datatable formatSignif

print_table <- function(table, rownames = FALSE, digits = 3, ...){
  df <- DT::datatable(table, rownames = rownames, extensions = 'Buttons',
                      options = list(scrollX = TRUE,
                                     dom = '<<t>Bp>',
                                     buttons = c('copy','csv', 'excel', 'pdf', 'print')), ...)
  num_cols <- c(as.numeric(which(sapply(table, class) == "numeric")))
  if(length(num_cols) > 0){
    DT::formatSignif(df, columns = num_cols, digits = digits)
  } else{
    df
  }
}




#' Compute genomic prediction different methods
#'
#' @param geno genomic data with genotypes in row (GID in rownames) and marker in columns. Values should be column centered and scaled.
#' @param pheno phenotypic data with genotypes in row (in GID column) and traits in columns. Phenotypic value should be corrected for year and location effects beforehand.
#' @param traits character vector of trait names
#' @param GP.method character vector of length one of genomic prediction methods to use. Must be one of "rrBLUP", "GBLUP", "RKHS", "RKHS-KA", "RandomForest", "BayesA", "BayesB" or "LASSO".
#' @param nreps number of repetitions for cross-validation, default is 10
#' @param nfolds number of folds for cross-validation, default is 10
#' @param h bandwith parameter for RKHS, default is 1. If multiple values are provided, method will be RKHS Kernel Averaging.
#' @param nIter number of iterations for RKHS, default is 8000
#' @param burnIn number of burn-in iterations for RKHS, default is 1500
#' @param ntree number of trees for RandomForest, default is 100
#' @param nb.mtry number of mtry for RandomForest, default is 10
#' @param p2d.temp path to directory to export temporary genomic prediction results, default is NULL (could cause error in parallelization if NULL).
#' @param nb.cores number of cores to parallelize the computation, default is 1 (no parallelization)
#' @param p2f.stats path to file to export genomic prediction results, default is NULL
#'
#' @return a list with the following elements: `obspred` with observed vs. predicted genotypic values, and `gp.stats` with genomic prediction statistics
#' @author Prediction5 Team
#' @seealso [getFolds()]
#' @importFrom BGLR BGLR
#' @importFrom parallel makeCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach %dopar%
#' @importFrom rrBLUP kinship.BLUP kin.blup
#' @importFrom caret train predict.train
#' @importFrom randomForest randomForest
#' @importFrom dplyr bind_rows
#' @importFrom stats na.omit cor
#' @importFrom utils write.csv write.table
#' @importFrom glmnet cv.glmnet
#' @importFrom purrr map_dfr
#' @export

compute_GP_methods <- function(geno, pheno, traits, GP.method, nreps=10,
                               nfolds=10, h=1, nb.mtry=10, nIter=8000,burnIn=1500,
                               ntree=100,p2d.temp=NULL,
                               nb.cores=1, p2f.stats=NULL){


  stopifnot(GP.method %in% c("rrBLUP","GBLUP","RKHS","RKHS-KA","RandomForest",
                             "BayesA","BayesB","BayesC","LASSO"),
            all(is.numeric(h)), length(GP.method) == 1,
            "GID" %in% colnames(pheno))
  if(!is.null(p2d.temp) && !dir.exists(p2d.temp) &
     GP.method %in% c("RKHS","RKHS-KA","BayesA","BayesB","BayesC")){
    dir.create(p2d.temp, recursive = TRUE)
  }

  ## format inputs
  geno <- as.matrix(geno)
  pheno <- as.data.frame(pheno)
  int.geno <- intersect(rownames(geno),pheno$GID)
  geno <- geno[int.geno,]
  pheno <- pheno[match(int.geno, pheno$GID),c("GID",traits)]

  # prepare outputs
  ## create a data frame to store accuracy
  gp.stats <- expand.grid(trait=traits, rep=seq(nreps), corP=NA,
                          GP.method=GP.method, nb.genos.TS=NA,
                          nb.genos.VS=NA, nb.snps=ncol(geno))
  ## list to store observed / predicted breeding values
  pred.list <- list()


  ## for loop across traits
  for(tr in traits){
    write(paste0("Working on trait ",tr), stderr())
    ## format geno and pheno for the trait, after removing missing values
    pheno_tr <- stats::na.omit(pheno[,c("GID",tr)])
    inds <- intersect(pheno_tr$GID, rownames(geno))
    pheno_tr <- pheno_tr[match(inds, pheno_tr$GID),]
    geno_tr <- geno[inds,]

    if(length(inds) < 50){
      print(paste("Too few genotypes for ",tr,"trait"))
    } else{
      ## vector of phenotypic values
      y <- pheno_tr[,tr] ; names(y) <- pheno_tr$GID

      ## genomic relationship matrix
      G <- as.matrix(crossprod(geno_tr)/ncol(geno_tr))


      ## create cluster
      cl <- parallel::makeCluster(type="SOCK",spec=nb.cores)
      doSNOW::registerDoSNOW(cl)

      for(r in 1:nreps){
        write(r, stderr())
        ## get cross-validation folds
        folds <- getFolds(nrow(pheno_tr), nb.folds=nfolds, seed=r+57425)


        if(GP.method %in% "rrBLUP"){
          requireNamespace("rrBLUP")
          ## parallel computation of GP for rrBLUP
          out <- foreach::foreach(f=1:nfolds, .packages="rrBLUP") %dopar%{
            ## estimate marker effects on training set
            res <- rrBLUP::kinship.BLUP(y=y[-folds[[f]]],
                                        G.train=geno_tr[-folds[[f]],],
                                        G.pred=geno_tr[folds[[f]],], K.method="RR")$g.pred
            return(res)
          }

        } else if(GP.method %in% "GBLUP"){
          requireNamespace("rrBLUP")
          ## parallel computation of GP for GBLUP
          out <- foreach::foreach(f=1:nfolds, .packages="rrBLUP") %dopar%{
            ## compute the genomic relationship matrix
            K <- tcrossprod(geno_tr)/ncol(geno_tr)

            ## function inputs
            data <- pheno_tr
            ## set phenotypes to NA for the training set
            data[[tr]][folds[[f]]] <- NA
            ## estimate marker effects on training set
            res <- rrBLUP::kin.blup(data=data,geno="GID",pheno=tr,K=K)
            return(res$pred[folds[[f]]])
          }


        } else if(GP.method %in% "RKHS" & length(h) ==1){
          requireNamespace("BGLR")
          ## define distance matrix for RKHS
          D <- as.matrix(dist(geno_tr,method="euclidean"))^2
          D <- D/mean(D)
          ### compute kernel
          K <- exp(-h*D)

          ## parallel processing
          out <- foreach::foreach(f=1:nfolds,.errorhandling = "pass",.packages="BGLR") %dopar% {
            ## estimate marker effects on training set
            yNA <- y
            yNA[folds[[f]]] <- NA
            if(!is.null(p2d.temp)){
              p2f.temp <- paste0(p2d.temp,"/",GP.method,"_",tr,
                                 "_",r,"_",f)
            }
            fit <- BGLR::BGLR(y=yNA,ETA=list(list(K=K,model='RKHS')),
                              nIter=nIter,burnIn=burnIn,
                              saveAt=p2f.temp,verbose=FALSE)
            return(as.data.frame(fit$yHat[folds[[f]]]))
          }

        } else if(GP.method %in% "RKHS-KA" | length(h) >1){
          requireNamespace("BGLR")
          ## define distance matrix for RKHS
          D <- as.matrix(dist(geno_tr,method="euclidean"))^2
          D <- D/mean(D)
          KList <- list()
          for(i in 1:length(h)){
            KList[[i]]<-list(K=exp(-h[i]*D),model='RKHS')
          }

          ## parallel processing
          out <- foreach::foreach(f=1:nfolds,.errorhandling = "pass",.packages="BGLR") %dopar% {
            ## estimate marker effects on training set
            yNA <- y
            yNA[folds[[f]]] <- NA
            if(!is.null(p2d.temp)){
              p2f.temp <- paste0(p2d.temp,"/",GP.method,"_",
                                 tr,"_",r,"_",f)
            }
            fit <- BGLR::BGLR(y=yNA,ETA=KList,
                              nIter=nIter,burnIn=burnIn,
                              saveAt=p2f.temp,verbose=FALSE)
            return(as.data.frame(fit$yHat[folds[[f]]]))
          }


        } else if(GP.method %in% "RandomForest"){
          requireNamespace("caret")

          ## optimize mtry = number of randomly selected variables at each split
          #tunegridrf <- expand.grid(.mtry=seq(1,ncol(geno_tr)/3,length.out=nb.mtry))
          tunegridrf <- expand.grid(.mtry=c(100,500,1000,2000,5000))
          ## parallel processing
          out <- foreach::foreach(f=1:nfolds,.errorhandling="pass",.packages="caret") %dopar% {
            ## estimate marker effects on training set

            fit <- caret::train(y=y[-folds[[f]]],
                                x=geno_tr[-folds[[f]],],
                                method = "rf",tuneGrid = tunegridrf, ntree=ntree)
            return(as.data.frame(caret::predict.train(fit,geno_tr[folds[[f]],])))
          }
        }else if(GP.method %in% "LASSO"){
          requireNamespace("glmnet")
          out <- foreach::foreach(f=1:nfolds, .packages="glmnet") %dopar% {
            ### Fit the model with glmnet package
            cv <- glmnet::cv.glmnet(y=y[-folds[[f]]],x=geno_tr[-folds[[f]],], alpha=1)
            fit <- glmnet::glmnet(y=y[-folds[[f]]],x=geno_tr[-folds[[f]],],
                                  alpha=1,lambda=cv$lambda.min)
            return(geno_tr[folds[[f]],] %*% as.matrix(fit$beta))
          }

        } else if(GP.method %in% "BayesA"){
          requireNamespace("BGLR")
          out <- foreach::foreach(f=1:nfolds,.packages="BGLR") %dopar% {

            if(!is.null(p2d.temp)){
              p2f.temp <- paste0(p2d.temp,"/",GP.method,"_",tr,"_",r,"_",f)
            }
            ## estimate marker effects on training set
            fit <- BGLR::BGLR(y=y[-folds[[f]]],
                              ETA=list(list(X=geno_tr[-folds[[f]],],model='BayesA')),
                              nIter=nIter,burnIn=burnIn,
                              saveAt=p2f.temp,verbose=FALSE)
            return(geno_tr[folds[[f]],] %*% as.matrix(fit$ETA[[1]]$b))
          }

        } else if(GP.method %in% "BayesB"){
          requireNamespace("BGLR")
          out <- foreach::foreach(f=1:nfolds,.packages="BGLR") %dopar% {
            if(!is.null(p2d.temp)){
              p2f.temp <- paste0(p2d.temp,"/",GP.method,"_",tr,"_",r,"_",f)
            }
            ## estimate marker effects on training set
            fit <- BGLR::BGLR(y=y[-folds[[f]]],
                              ETA=list(list(X=geno_tr[-folds[[f]],],model='BayesB')),
                              nIter=nIter,burnIn=burnIn,
                              saveAt=p2f.temp,verbose=FALSE)
            ## return predicted value on validation set for this fold
            return(geno_tr[folds[[f]],] %*% as.matrix(fit$ETA[[1]]$b))
          }
        }else if(GP.method %in% "BayesC"){
          requireNamespace("BGLR")
          out <- foreach::foreach(f=1:nfolds,.packages="BGLR") %dopar% {
            if(!is.null(p2d.temp)){
              p2f.temp <- paste0(p2d.temp,"/",GP.method,"_",tr,"_",r,"_",f)
            }
            ## estimate marker effects on training set
            fit <- BGLR::BGLR(y=y[-folds[[f]]],
                              ETA=list(list(X=geno_tr[-folds[[f]],],model='BayesC')),
                              nIter=nIter,burnIn=burnIn,
                              saveAt=p2f.temp,verbose=FALSE)
            ## return predicted value on validation set for this fold
            return(geno_tr[folds[[f]],] %*% as.matrix(fit$ETA[[1]]$b))
          }

        } ## end if method


        ## combine all predicted values across folds
        ### (too few individuals to calculate predictive ability for each fold)
        pred.obs.all <- purrr::map_dfr(out,as.data.frame)
        pred.obs.all <- data.frame(GID=rownames(pred.obs.all),
                                   obs=pheno[match(rownames(pred.obs.all),pheno$GID),tr],
                                   pred=pred.obs.all[,1],
                                   trait=tr,rep=r, GP.method=GP.method)
        pred.list[[tr]][[r]] <- pred.obs.all


        ## estimnate predictive ability as Pearson's correlation between observed and predicted
        idx <- which(gp.stats$trait %in% tr & gp.stats$rep %in% r)
        gp.stats$corP[idx] <- round(cor(pred.obs.all$obs,pred.obs.all$pred),3)
        gp.stats$nb.genos.TS[idx] <- length(y)-mean(sapply(folds, length))
        gp.stats$nb.genos.VS[idx] <- mean(sapply(folds, length))

      } ## end for rep

      parallel::stopCluster(cl)
    } ## end if ind>50

  } ## end for trait

  ## if file path provided, export results.
  if(!is.null(p2f.stats)){
    ext <- tools::file_ext(p2f.stats)
    if(ext %in% c("csv")){
      utils::write.csv(gp.stats, file=p2f.stats, row.names=FALSE)
    } else if(ext %in% c("txt","tsv")){
      utils::write.table(file=p2f.stats, gp.stats, sep="\t", row.names=FALSE)
    } else if(ext %in% "rds"){
      saveRDS(gp.stats, file=p2f.stats)
    }

  }

  return(list(obspred=pred.list, gp.stats=gp.stats))

}


#' Run genomic prediction on all genotypes, output predicted values
#'
#' @param geno genomic data with genotypes in row (GID in rownames) and marker in columns. Values should be column centered and scaled.
#' @param pheno phenotypic data with genotypes in row (in GID column) and traits in columns
#' @param traits character vector of trait names, must correspond to `pheno` column names
#' @param GP.method character vector of genomic prediction methods to use, must be one of "rrBLUP", "RKHS", "BayesA", "BayesB"
#' @param runCV logical, if TRUE, run cross-validation on the common genotypes in `pheno` and `geno`, default is FALSE
#' @param testSetGID GID of the test set genotypes, if NULL, will provide predicted genotypic values for all genotypes in `geno`
#' @param nreps number of repetitions for cross-validation, default is 10
#' @param nfolds number of folds for cross-validation, default is 10
#' @param h bandwith parameter for RKHS, default is 1.
#' @param nb.mtry number of randomly selected variables at each split for RandomForest, default is 10
#' @param nIter number of iterations for RKHS, default is 8000
#' @param burnIn number of burn-in iterations for RKHS, default is 1500
#' @param ntree number of trees for RandomForest, default is 100
#' @param p2d.temp path to directory to export temporary genomic prediction results, default is NULL (could cause error in parallelization if NULL).
#' @param nb.cores number of cores to parallelize the computation, default is 1 (no parallelization)
#' @param p2f.stats.cv path to file to export cross-validation genomic prediction results, default is NULL
#' @param p2f.pred path to file to export genotypic values from genomic prediction, default is NULL
#' @param verbose integer, level of verbosity, default is 1
#'
#' @return a data frame with the predicted values for all genotypes and traits
#' @seealso [compute_GP_methods()]
#' @author Prediction5 Team
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach %dopar% foreach
#' @importFrom rrBLUP kinship.BLUP
#' @importFrom tidyr pivot_wider
#' @importFrom tidyselect all_of
#' @export

compute_GP_allGeno <- function(geno, pheno, traits, GP.method,
                               runCV=FALSE, testSetGID=NULL,
                               nreps=10,nfolds=10, h=1, nb.mtry=10,
                               nIter=8000,burnIn=1500,
                               ntree=100,p2d.temp=NULL,
                               nb.cores=1, p2f.stats.cv=NULL,
                               p2f.pred=NULL, verbose=1){

  ## initial verifications
  stopifnot(GP.method %in% c("rrBLUP","GBLUP","RKHS","BayesA","BayesB","BayesC"),
            #"LASSO","RKHS-KA","RandomForest"),
            all(is.numeric(h)), length(GP.method) == 1,
            "GID" %in% colnames(pheno),
            all(traits %in% colnames(pheno)))
  if(!is.null(p2d.temp) && !dir.exists(p2d.temp) &
     GP.method %in% c("RKHS","RKHS-KA","BayesA","BayesB","BayesC")){
    dir.create(p2d.temp, recursive = TRUE)
  }

  ## formatting
  geno <- as.matrix(geno)
  ntraits <- length(traits)
  pheno.traits <- pheno[,c("GID",traits)]

  ## cross-validation
  if(runCV){
    cmonGID <- intersect(rownames(geno), pheno$GID)
    if(verbose > 0){
      print(paste0("Compute genomic prediction in cross-validation for ",
                   length(cmonGID)," common genotypes"))
    }
    ## use custom function
    gp.cv <- compute_GP_methods(
      geno = geno[cmonGID,],
      pheno = pheno[match(cmonGID, pheno$GID),],
      traits = traits,
      nfolds = nfolds,
      nreps = nreps,
      GP.method = GP.method,
      h = h,
      nIter = nIter,
      burnIn = burnIn,
      p2d.temp=p2d.temp,
      nb.cores=nb.cores,
      p2f.stats=p2f.stats.cv)
  }


  ## create cluster
  cl <- parallel::makeCluster(type="SOCK",spec=nb.cores)
  doSNOW::registerDoSNOW(cl)
  ## run genomic prediction on all genotypes

  if(GP.method %in% "rrBLUP"){
    ## parallel computation of GP for rrBLUP
    out <- foreach::foreach(f=1:ntraits,.errorhandling = "pass",
                            .combine = "rbind",.packages = "rrBLUP") %dopar% {

                              ## estimate marker effects on all available data
                              yNA <- merge(data.frame(GID=c(rownames(geno))),
                                           pheno.traits[,c("GID",traits[f])], by="GID", all.x=TRUE)
                              yGID <- yNA$GID
                              ## for test set geno if provided, set to NA (remove from training set)
                              if(!is.null(testSetGID)){
                                yNA[yNA$GID %in% testSetGID,tr] <- NA
                              }
                              fit <- rrBLUP::kinship.BLUP(y=yNA[[traits[f]]],
                                                          G.train=geno,
                                                          G.pred=geno, K.method="RR")$g.pred
                              ## output predicted values
                              tmp <- data.frame(GID=yGID, GP.method=GP.method,trait=traits[f],yHat=fit)
                              return(tmp)
                            }
  } else if(GP.method %in% "GBLUP"){
    requireNamespace("rrBLUP")
    ## parallel computation of GP for GBLUP
    out <- foreach::foreach(f=1:ntraits, .packages="rrBLUP",.errorhandling = "pass",
                            .combine = "rbind") %dopar%{
                              ## estimate marker effects on all available data
                              yNA <- merge(data.frame(GID=c(rownames(geno))),
                                           pheno.traits[,c("GID",traits[f])], by="GID", all.x=TRUE)
                              yGID <- yNA$GID

                              ## compute the genomic relationship matrix
                              K <- tcrossprod(geno[yGID,])/ncol(geno)

                              ## estimate marker effects on training set
                              res <- rrBLUP::kin.blup(data=yNA,geno="GID",pheno=traits[f],K=K)
                              tmp <- data.frame(GID=yGID, GP.method=GP.method,
                                                trait=traits[f],yHat=res$pred)
                              return(tmp)
                            }

  }else if(GP.method == "RKHS"){
    D <- as.matrix(dist(geno,method="euclidean"))^2
    D <- D/mean(D, na.rm=TRUE)
    ### compute kernel
    K <- exp(-h*D)
    rm(D)

    if(is.null(p2d.temp)){
      print(warning("No temporary directory provided, this may cause issues."))
    }

    ## parallel processing over traits
    out <- foreach::foreach(
      f=1:ntraits,.errorhandling = "pass",
      .combine = "rbind",.packages = "BGLR") %dopar% {
        yNA <- merge(data.frame(GID=c(rownames(geno))),
                     pheno.traits[,c("GID",traits[f])],
                     by="GID", all.x=TRUE)
        yNA <- yNA[match(rownames(K),yNA$GID),]
        ## if testSetGID provided, set to NA (remove from training set)
        if(!is.null(testSetGID)){
          yNA[yNA$GID %in% testSetGID,traits[f]] <- NA
        }
        stopifnot(identical(yNA$GID, rownames(geno)))

        p2f.temp <- paste0(p2d.temp,"/RKHS_PredAll_testSet_",traits[f])
        ## fit model
        fit <- BGLR::BGLR(y=yNA[[traits[f]]],
                          ETA=list(list(K=K,model='RKHS')),
                          nIter=nIter,burnIn=burnIn,
                          saveAt=p2f.temp,verbose=FALSE)
        ## output predicted values
        tmp <- data.frame(GID=yNA$GID, GP.method=GP.method,
                          trait=traits[f],yHat=fit$yHat)
        return(tmp)
      }

  } else if(GP.method %in% c("BayesA","BayesB","BayesC")){
    if(is.null(p2d.temp)){
      print(warning("No temporary directory provided, this may cause issues."))
    }
    ## parallel processing over traits
    out <- foreach::foreach(
      f=1:ntraits,.errorhandling = "pass",
      .combine = "rbind",.packages = "BGLR") %dopar% {
        yNA <- merge(data.frame(GID=c(rownames(geno))),
                     pheno.traits[,c("GID",traits[f])], by="GID", all.x=TRUE)

        yNA <- yNA[match(rownames(K),yNA$GID),]
        ## if testSetGID provided, set to NA (remove from training set)
        if(!is.null(testSetGID)){
          yNA[yNA$GID %in% testSetGID,traits[f]] <- NA
        }
        stopifnot(identical(yNA$GID, rownames(geno)))

        p2f.temp <- paste0(p2d.temp,"/",GP.method,"_PredAll_testSet_",traits[f])
        ## fit model
        fit <- BGLR::BGLR(y=yNA[[traits[f]]],
                          ETA=list(list(X=geno,model=GP.method)),
                          nIter=nIter,burnIn=burnIn,
                          saveAt=p2f.temp,verbose=FALSE)
        ## output predicted values
        tmp <- data.frame(GID=yNA$GID, GP.method=GP.method,
                          trait=traits[f],yHat=fit$yHat)
        return(tmp)
      }
  }

  parallel::stopCluster(cl)
  ## set the predicted values in a wide format, with one column per trait
  all.pred <- out %>%
    dplyr::select(GID,GP.method,trait,yHat) %>%
    tidyr::pivot_wider(names_from="trait", values_from="yHat")

  if(!is.null(testSetGID)){
    ## if testSetID provided, keep only the genotypes in the test set
    all.pred <- all.pred[all.pred$GID %in% testSetGID,]
  }


  ## export predicted values
  if(!is.null(p2f.pred)){
    ## if file path provided, export results.
    if(!is.null(p2f.stats)){
      ext <- tools::file_ext(p2f.pred)
      if(ext %in% c("csv")){
        utils::write.csv(all.pred, file=p2f.pred, row.names=FALSE)
      } else if(ext %in% c("txt","tsv")){
        utils::write.table(file=p2f.pred, all.pred, sep="\t",
                           row.names=FALSE)
      } else if(ext %in% "rds"){
        saveRDS(all.pred, file=p2f.pred)
      }
    }

  } ## end if p2f.pred

  if(runCV){
    return(list(obspred=gp.cv$obspred, gp.stats=gp.cv$gp.stats,
                all.pred=all.pred))
  } else {
    return(all.pred) ## doesn't work
  }
}

