#!/usr/bin/Rscript

##' @name qMAPModules
##' @rdname qMAPModules
##' @title qMAP -quantitative mapping for differential fragmentation of small non-coding RNAs
##' @details
##' For more information, see the user manual:
##' \code{browseURL(system.file("docs/qMAP_v1.0.0_User_Manual.pdf", package = "qMAPModules"))}
##'
##' After installing the package qMAPModules, call the library as: library(qMAPModules).
##'
##' There are two modules or main functions in this package : qMAP and qMAP_mh. qMAP has two model/methods: Model1 and Model2
##'
##' Ask help for the description and help menus of each of these functions :
##' \cr ?qMAP.single_family
##' \cr ?qMAP.MH.single_family
##'
##' If want to run examples, it is necessary to set working directory to the base folder of the installed package qMAPModules. Such as : setwd("/home/..../R/.../qMAPModules/").
##'
##' @param matrix_file A user-provided read-count matrix where the first and second columns are labeled as "Sequence" and "Annotation", respectively. The remaining columns contain read counts from different samples under study.
##' @param cl A sample classification vector allowing only two groups: 1 for controls, 2 for cases, and -1 for the samples to be excluded from the analysis.
##' @param parental_rna_filename A FASTA file containing sequences of parental RNA families description.
##' @param sncrna_family This parameter specifies a single sncRNA family name to search for in the user-provided parental RNA FASTA file. If a match is found, differential fragmentation analysis will be performed for that family. Use any single family name, such as : “mature-tRNA-Ala-AGC”, “tRNA-Ala-AGC”, “mature-mt_tRNA-His-GTG”, “mt_tRNA-His-GTG”, “5S-rRNA”, “16S-rRNA”
##' @param max.mismatch This parameter specifies the maximum number of allowed mismatches (default is 1) when mapping sequence reads onto the FASTA sequence of the parental RNA.
##' @return with qMAP.single_family, it will return, for each queried parental RNA (sncrna_family), an output data-frame with  mean_diff (mean difference in coverage between two groups), p (the P-value indicating statistical significance ), and species_num (number of contributing sncRNA species mapped to the parental RNA).
##' @return with qMAP.MH.single_family, it will return, for each queried parental RNA (sncrna_family), an output data-frame with sequence (the sncRNA species sequence), or (the odds ratio), p (the P-value inferred from the Mantel-Haenszel statistic), and start (the start position on the parental RNA where the sncRNA species is mapped).
#'
#' @importFrom Biostrings matchPattern start
#' @importFrom dplyr %>% arrange
#' @importFrom rlang .data
#' @importFrom seqinr read.fasta
#' @importFrom stats mantelhaen.test lm t.test coef pt
#' @importFrom usethis use_package
#' @importFrom utils read.delim combn write.table
NULL
#'
#' Run qMAP on sncRNA matrix file (file-based input)
#'
#' This is a variant for 'qMAP.single_family( )' that reads the input matrix from a file.
#'
#' @param input_file Path to a tab-delimited read-count matrix. The first two columns should be "Sequence" and "Annotation", and the rest are counts.
#' @param ... Additional arguments passed to `qMAP.single_family()`: `cl`, `sncrna_family`, `parental_rna_filename`, etc.
#'
#' @return A data frame with results from the qMAP model
#' @export
qMAP.single_family.1 <- function(input_file, ...) {
  qMAP.single_family(function() read.delim(input_file), ...)
}
#' Run qMAP on sncRNA data frame (in-memory)
#'
#' This is a variant for 'qMAP.single_family( )' where the matrix is already in memory as a data frame.
#'
#' @param input_dataframe A data frame with "Sequence", "Annotation", and read-count columns.
#' @param ... Additional arguments passed to `qMAP.single_family()`: `cl`, `sncrna_family`, `parental_rna_filename`, etc.
#'
#' @return A data frame with results from the qMAP model
#' @export
qMAP.single_family.2 <- function(input_dataframe, ...) {
  qMAP.single_family(function() input_dataframe, ...)
}
#' qMAP.single_family - identification of differentially fragmented parental RNAs between groups through quantitative mapping of sncRNAs
#'
#' @param matrix_file A user-provided read-count matrix where the first and second columns are labeled as "Sequence" and "Annotation", respectively. The remaining columns contain read counts from different samples under study.
#' @param cl A sample classification vector allowing only two groups: 1 for controls, 2 for cases, and -1 for the samples to be excluded from the analysis.
#' @param parental_rna_filename A FASTA file containing sequences of parental RNA families description.
#' @param sncrna_family This parameter specifies a single sncRNA family name to search for in the user-provided parental RNA FASTA file. If a match is found, differential fragmentation analysis will be performed for that family. Use any single family name, such as : “mature-tRNA-Ala-AGC”, “tRNA-Ala-AGC”, “mature-mt_tRNA-His-GTG”, “mt_tRNA-His-GTG”, “5S-rRNA”, “16S-rRNA”
#' @param max.mismatch This parameter specifies the maximum number of allowed mismatches (default is 1) when mapping sequence reads onto the FASTA sequence of the parental RNA.
#' @param method Specify either “1” or “2” to select which model to run for identifying differentially fragmented parental RNAs. Use “2” to perform analysis with Model 2 (method = 2) or use “1” (default) to run Model 1 (method = 1).
#' @param min_count This parameter sets the minimum allowed mean read-count across samples for the sncRNA family to be considered in the analysis. The default value is 10, meaning only families with mean read-counts greater than or equal to 10 across samples will be included.
#' @param shuffling_round This parameter specifies the minimum number of sample shuffling rounds performed during the statistical model to identify differential fragmentation. The default value is 100 for Model 1 (method =1),  and 1000 for Model 2 (method =2)).
#' @param output  This parameter specifies the name of the output file, By default, if no name is provided,  the output file will be saved in the working directory with a name like “qMAP_****_output.txt", where * represents various  parameters used in the run. If '"temp"', the function creates a temporary output file (used for examples/testing only).
#' @return with qMAP, it will return, for each queried parental RNA (sncrna_family), an output data-frame with  mean_diff (mean difference in coverage between two groups), p (the P-value indicating statistical significance ), and species_num (number of contributing sncRNA species mapped to the parental RNA).
#' @note In example runs, `output = "temp"` is used to avoid writing files during checks. In real analyses, you may provide a file path or let the function auto-name the output (see manual).
#'
#' @examples
#' # Note:  After installation, one can find the example files in "../qMAPModules/extdata/"
#' # To run on a specific input, provide the full path as: qMAPModules(a = "/path/to/your_matrix.txt")
#' sample_matrix_1 <- system.file("extdata", "sample_matrix1.txt", package = "qMAPModules")
#' sample_matrix_2 <- system.file("extdata", "sample_matrix2.txt", package = "qMAPModules")
#' sample_matrix_3 <- system.file("extdata", "sample_matrix3.txt", package = "qMAPModules")
#' mouse_gtsrna_list <- system.file("extdata", "mmu_gtsrna.txt", package = "qMAPModules")
#' mouse_rRNAs <- system.file("extdata", "mouse_rRNAs.fa", package = "qMAPModules")
#' mouse_tRNAs <- system.file("extdata", "mm10-tRNAs_CCA.fa", package = "qMAPModules")
#' # Run differential fragmentation analysis for 5S-rRNA family using  corresponding FASTA sequence
#' # in the user provided parental rna file (mouse_rRNAs.fa) by mapping sequencing reads in the query
#' # read-count matrix file (sample_matrix1.txt) comprising read counts for 5 control (cl = 1) and
#' # 4 case (cl =2) samples. It will run this for method 1 (default)  with allowed maximum mismatch 1
#' # (default) and with 100 shuffling rounds (default, for method 1), only if the mean total read-
#' # counts across the samples for the mapped reads (against query sncRNA family) > 10 (default).
#' qMAP.single_family.1(input_file = sample_matrix_1, cl=c(1,1,1,1,1,2,2,2,2),
#'                     sncrna_family="5S-rRNA", parental_rna_filename= mouse_rRNAs, output = "temp")
#'
#' # Run differential fragmentation analysis in a similar manner as for example 1 above.
#'qMAP.single_family.1(sample_matrix_1, cl=c(1,1,1,1,1,2,2,2,2), sncrna_family="5S-rRNA",
#'			parental_rna_filename=mouse_rRNAs, method = "1", max.mismatch =1,
#'			shuffling_round = 100, min_count = 10, output = "temp")
#'
#' # Run differential fragmentation analysis for 16S-rRNA family using  corresponding FASTA sequence
#' # in the user provided parental rna file (mouse_rRNAs.fa) by mapping sequence reads in the query
#' # read-count matrix file (sample_matrix2.txt) comprising read counts for 3 control (first three
#' # sample columns in the matrix) (cl = 1) and 3 case (last three sample columns in the matrix)
#' # (cl = 2) samples. Here, 3 samples excluded (cl = -1). It will run this for method 1 (default) with
#' # allowed maximum mismatch 1 (default) and with 100 shuffling rounds (default, for method 1), only
#' # if the mean total read-counts across the samples for the mapped reads (against query sncRNA
#' # family) > 10 (default). The result output will be saved to the specified file.
#' qMAP.single_family.1(sample_matrix_2, cl=c(1,1,1,-1,-1,-1,2,2,2), sncrna_family="16S-rRNA",
#' 			parental_rna_filename = mouse_rRNAs,  output = "temp")
#'
#' # Run differential fragmentation analysis for tRNA-Gly-CCC family by applying method 2 using
#' # corresponding FASTA sequence in the user provided parental rna file (mm10-tRNAs_CCA.fa) by mapping
#' # sequence reads in the query read-count matrix file (sample_matrix1.txt) comprising read counts for
#' # 5 control (cl = 1) and 4 case (cl = 2) samples. It will run this with allowed maximum mismatch 1
#' # (default) and with 1000 shuffling rounds (default, for method 2), only if the mean total read-
#' # counts across the samples for the mapped reads (against query sncRNA family) > 10 (default).
#' qMAP.single_family.1(sample_matrix_1, cl=c(1,1,1,1,1,2,2,2,2), sncrna_family="tRNA-Gly-CCC",
#' 			parental_rna_filename = mouse_tRNAs, method = "2", output = "temp")
#'
#' # Run differential fragmentation analysis for tRNA-Gly-CCC family by applying method 2 using
#' # corresponding FASTA sequence in the user provided parental rna file (mm10-tRNAs_CCA.fa) by mapping
#' # sequence reads in the query read-count matrix file (sample_matrix1.txt) comprising read counts for
#' # 5 control (cl = 1) and 4 case (cl = 2) samples. It will run this with allowed maximum mismatch 1
#' # (default) and with 1000 shuffling rounds (default, for method 2), only if the mean total read-
#' # counts across the samples for the mapped reads (against query sncRNA family) > 2.
#' qMAP.single_family.1(sample_matrix_1, cl=c(1,1,1,1,1,2,2,2,2), sncrna_family="tRNA-Gly-CCC",
#' 		parental_rna_filename = mouse_tRNAs, min_count = 2, method = "2", output = "temp")
#'
#' # Run differential fragmentation analysis for tRNA-Val-AAC family by applying method 2 using
#' # corresponding FASTA sequence in the user provided parental rna file (mm10-tRNAs_CCA.fa) by
#' # mapping sequence reads (with mismatch zero) in the query read-count matrix file
#' # (sample_matrix3.txt) comprising read counts for 4 Contioin1 (cl = 1), 4 Condition2 (cl = 2)
#' # samples, and rest 8 for other conditions. It will run this with 1000 shuffling rounds (default,
#' # for method 2), only if the mean total read-counts across the samples for the mapped reads (against
#' # query sncRNA family) > 10 (default).
#' qMAP.single_family.1(sample_matrix_3, cl=c(1,1,1,1,2,2,2,2,-1,-1,-1,-1,-1,-1,-1,-1),
#' 			sncrna_family="tRNA-Val-AAC", parental_rna_filename = mouse_tRNAs,
#' 			max.mismatch = 0, method = "2", output = "temp")
#'
#' # Running with query read-count matrix as a data frame as input
#' # First read and save the read-count matrix file as a dataframe:
#' a <- read.delim(sample_matrix_1, header = TRUE)
#' # Then call qMAP function to run differential fragmentation analysis for tRNA-His-GTG family by
#' # applying method 2 using corresponding FASTA sequence in the user provided parental rna file
#' # (mm10-tRNAs_CCA.fa) by mapping sequence reads (with mismatch 1) in the query read-count matrix
#' # file (a, sample_matrix1.txt) comprising read counts for 5 control (cl = 1) and 4 case (cl =2)
#' # samples. It will run this with 1000 shuffling rounds (default, for method 2), only if the mean
#' # total read-counts across the samples for the mapped reads (against query sncRNA family) > 10
#' # (default).
#' qMAP.single_family.2(a, cl=c(1,1,1,1,1,2,2,2,2), sncrna_family="tRNA-His-GTG",
#' 		parental_rna_filename= mouse_tRNAs, max.mismatch = 1, method = "2", output = "temp")
#'
#' # Run differential fragmentation analysis for each of the RNA families (listed in a file; here it
#' # is mmu_gtsrna.txt) using  corresponding FASTA sequences in the user provided parental rna file
#' # (mm10-tRNAs_CCA.fa) by mapping sequence reads in the query read-count matrix file
#' # (sample_matrix2.txt) comprising read counts for 3 control (first three sample columns in the
#' # matrix) (cl = 1) and 3 case (next three sample columns in the matrix) (cl = 2) samples. Here,
#' # last 3 samples excluded (cl = -1). It will run this for method 1 with allowed maximum mismatch 1
#' # (default) and with 100 shuffling rounds, only if the mean total read-counts across the samples
#' # for the mapped reads (against query sncRNA family) > 10 (default). The result output will be
#' # saved to the specified file, qmap_multifamily_results.txt.
#' # To run qMAP for multiple families, it is suggested to run with the second variant
#' # (qMAP.single_family.2) in a for loop as follows:
#' # First read and save the read-count matrix file as a dataframe
#' a = read.delim(sample_matrix_2)
#' cl = c(1,1,1,2,2,2,-1,-1,-1)
#' result = c()
#' # Read sncRNA families from list in the file
#' sncrna_family_list = read.table(mouse_gtsrna_list)[,1]
#' # Then for each family in the list, call qMAP function
#' for (i in 1:length(sncrna_family_list)) {
#' 	result = rbind(result, qMAP.single_family.2(a, cl=cl, sncrna_family=sncrna_family_list[i],
#' 	parental_rna_filename = mouse_tRNAs, method=1, shuffling_round=100, output="temp"))
#' 	}
#' # Creating a safe temporary output file exclusively for run examples (see manual)
#' out_file <- tempfile("qmap_multifamily_results", fileext = ".txt")
#' write.table(result, file = out_file, quote=FALSE, row.names=FALSE, sep="\t")
#'
#' @export

### Function to run qMAP for a single sncRNA family, calling Model 1 (default, method = 1), or Model 2 (if, method =2)
### cl: sample classification vector, only two groups are allowed (1 for controls, 2 for cases, and -1 for the samples to be excluded)
qMAP.single_family <- function(matrix_file, cl, sncrna_family, parental_rna_filename, max.mismatch = getOption("max.mismatch", 1), min_count = getOption("min_count", 10), shuffling_round =  NULL, method = 1, output = NULL)
{

  #### qmap function body starts ####
  args = commandArgs(trailingOnly=TRUE)

  ########## Part1: defining different functions to call ##########

  options(warn=1)


  # Function to check if all values in the matrix are whole numbers
  check_whole_numbers <- function(exp) {

    # Check for non-numeric values
    if (!all(sapply(exp, is.numeric))) {
      stop("All columns in the read-count matrix must be numeric. There are alphabets/characters in the matrix.")
    }
    ### Check for non-finite
    if (!all(sapply(exp, is.finite))) {
      stop("The read-count matrix contains non-finite values (NA, Inf, or NaN).")
    }

    ### Check for non-integer values
    if (!all(sapply(exp, is.integer))) {
      message("Warning: There are some non-integer values in the read-count matrix (e.g., 0.5, 1.01). Only whole numbers should be used.")
    }
    return(TRUE)
  }

  ## Function to calculate coverage pattern (sncRNA mapping) along parental RNA
  calc.coverage <- function(expr, start, len, large_rna_length)
  {
    coverage = rep(0, large_rna_length)
    for (i in 1:length(expr))
    {
      if (is.na(start[i]) || is.na(len[i])) next  ## skip invalid entries

      for (j in start[i]:(start[i]+len[i]-1))
      {
        if (j > 0 && j <= large_rna_length) {  ## avoid index out of bound
          coverage[j] = coverage[j] + expr[i]
        }
      }
    }
    if (sum(coverage) == 0) return(coverage)
    return(coverage / sum(coverage))
  }

  ## Function to calculate t-statistics between two coverage
  calc.t_statistic <- function(coverage1, coverage2)
  {
    t_list = rep(NA, nrow(coverage1))
    for (i in 1:nrow(coverage1))
    {
      t_list[i] = abs(t.test(coverage1[i,], coverage2[i,])$statistic)
    }
    return(t_list)
  }

  ## Function to run Model 1, a statistical model to identify differential fragmentation (reshuffling + parametric test)
  diff_frag.model1 <- function(coverage, cl, shuffling_round = 100)
  {
    coverage1 = coverage[,cl==1]
    coverage2 = coverage[,cl==2]
    t_list = calc.t_statistic(coverage1, coverage2)
    n = length(cl)
    n1 = length(cl[cl==1])
    n2 = length(cl[cl==2])
    rna_len = nrow(coverage)

    if (choose(n, n1) < shuffling_round)
    {
      cl_shuffled_matrix = matrix(NA, ncol=choose(n, n1), nrow=length(cl))
      combination = combn(1:n, n1)
      for (i in 1:ncol(combination))
      {
        cl_shuffled_matrix[combination[,i], i] = 1
        cl_shuffled_matrix[-combination[,i], i] = 2
      }
    } else {
      cl_shuffled_matrix = matrix(NA, ncol=shuffling_round, nrow=length(cl))
      for (i in 1:ncol(cl_shuffled_matrix))
      {
        cl_shuffled_matrix[,i] = cl[sample(n)]
      }
    }

    t_shuffled_list = rep(NA, rna_len * ncol(cl_shuffled_matrix))
    for (i in 1:ncol(cl_shuffled_matrix))
    {
      cl_shuffled = cl_shuffled_matrix[,i]
      coverage1_shuffled = coverage[,cl_shuffled==1]
      coverage2_shuffled = coverage[,cl_shuffled==2]
      t_shuffled_list[((i-1)*rna_len+1):(i*rna_len)] = calc.t_statistic(coverage1_shuffled, coverage2_shuffled)
    }

    data = data.frame(c(t_list, t_shuffled_list), rep(1:rna_len, ncol(cl_shuffled_matrix)+1), c(rep(2, rna_len), rep(1, rna_len*ncol(cl_shuffled_matrix))))
    colnames(data) = c("diff", "pos", "group")

    mod = lm(data$diff ~ factor(data$pos) + data$group)
    lm_summary = summary(mod)
    p = pt(coef(lm_summary)["data$group", 3], mod$df, lower.tail = FALSE)
    mean_diff = sum(abs(rowMeans(coverage1) - rowMeans(coverage2)))
    return(c(mean_diff, p))
  }

  ## Function to run Model 2, a statistical model to identify differential fragmentation (reshuffling) (non-parametric)
  diff_frag.model2 <- function(coverage, cl, shuffling_round = 1000)
  {
    coverage1 = coverage[,cl==1]
    coverage2 = coverage[,cl==2]
    mean_diff = sum(abs(rowMeans(coverage1) - rowMeans(coverage2)))

    n = length(cl)
    mean_diff_shuffled = rep(NA, shuffling_round)
    for (i in 1:shuffling_round)
    {
      cl_shuffled = cl[sample(n)]
      coverage1_shuffled = coverage[,cl_shuffled==1]
      coverage2_shuffled = coverage[,cl_shuffled==2]
      mean_diff_shuffled[i] = sum(abs(rowMeans(coverage1_shuffled) - rowMeans(coverage2_shuffled)))
    }
    p = length(mean_diff_shuffled[mean_diff_shuffled>=mean_diff]) / shuffling_round
    return(c(mean_diff, p))
  }

########## Part2: Processing of input matrix and different parameters to be used in the function ##########

  ### trim sncRNA family name and look for the corresponding fasta sequences in the user provided parental rna file
  sncrna_family = sub("mature-", "", sncrna_family)
  s = unlist(read.fasta(parental_rna_filename, seqtype="DNA", as.string=T, forceDNAtolower=F))
  if (grepl("rRNA", sncrna_family)) s = s[names(s) == sncrna_family]
  else s = s[grepl(sncrna_family, names(s))]
  uni_s = unique(s)
  z = match(uni_s, s)
  parental_rna_list = s[z]

  if (length(s) == 0)
  {
    print(paste0("The specified sncRNA family - ", sncrna_family,  ", does not exist in the provided parental RNA file."))
    return(NULL)
  }


  # If matrix_file is a function, call it
  if (is.function(matrix_file)) {
    a <- matrix_file()
  } else if (is.data.frame(matrix_file) || is.matrix(matrix_file)) {
    a <- matrix_file
  } else {
    stop("`matrix_file` must be a function or a data frame/matrix.")
  }

  ### Processing of input matrix and different parameters to be used in the function ###
  ### a: a user-provided read-count matrix with first and second column as "Sequence" and "Annotation", respectively. Rest of the columns having read counts from different samples under study. Sequence must be unique.###
  #a = matrix_file()
  anno = a[,1:2]  ### Separated the columns with unique Sequence (or ID) and their Annotation information
  e = a[,3:ncol(a)]  ### Separated the columns with read-count values

  check_whole_numbers(e)

  ### Confirming the sample size in cl option matched with that in the matrix file
  cl_sum = length(cl[cl==1]) + length(cl[cl==2]) + length(cl[cl==-1])
  if (ncol(e) != cl_sum)
  {
    print(paste0("Error: Number of samples (",   ncol(e), ") in the input matrix does not match the number of samples specified in the cl option (", cl_sum, ")."))
    return(NULL)
  }

  rownames(e) = a[,1]
  e = e[, cl==1 | cl==2]
  cl = cl[cl==1 | cl==2]
  if (length(cl[cl==1]) < 2 | length(cl[cl==2]) < 2) ### Confirming the sample size of both group is more than two
  {
    print("Error: sample size of at least one group is less than two.")
    return(NULL)
  }

  ### exclude the sncRNA species (rows in the matrix) with zero read-counts
  anno = anno[rowMeans(e) > 0,]
  e = e[rowMeans(e) > 0,]

  ### Processing the input matrix for query sncRNA family
  if (grepl("rRNA", sncrna_family))
  {
    is.sncrna_family = (anno$Annotation==sncrna_family)
  } else
  {
    if (grepl("tRNA", sncrna_family)) is.sncrna_family = grepl(paste("mature-", sncrna_family, sep=""), anno$Annotation)
    else is.sncrna_family = grepl(sncrna_family, anno$Annotation)
  }
  e_tmp = e[is.sncrna_family,]   ### Extracting read-count matrix values for the species with query sncRNA family annotation

  if (nrow(e_tmp) == 0)
  {
    print(paste0("There are no read fragments in the input matrix file annotated as - ", sncrna_family  , "."))
    return(NULL)
  }

  if (mean(colSums(e_tmp)) < min_count) ### proceed only for the sncRNA family with mean read-counts across samples > min_count (default 10)
  {
    print(paste0("Read fragments in the input matrix annotated as - ", sncrna_family  , ", have an average total read count per sample  < ", min_count, ". Try with less-stringent min_count threshold."))
    return(NULL)
  }

  ### Function to find the match position of the small sequence within the long sequence
  find_position_with_mismatches <- function(long_seq, small_seq) {
    match_result <- matchPattern(small_seq, long_seq, max.mismatch = max.mismatch)
    ### If there are matches, return the first match position (1-based index)
    if (length(match_result) > 0) {
      return(start(match_result)[1])  ### Return the first position
     } else {
      return(NA)  ### If no match is found
    }
  }

  mean_diff = rep(NA, length(parental_rna_list))
  p = rep(NA, length(parental_rna_list))
  species_num = rep(NA, length(parental_rna_list))
  ### running through the sequence length of parental rna
  for (i in 1:length(parental_rna_list))
  {
    start = mapply(find_position_with_mismatches, parental_rna_list[i], rownames(e_tmp))
    e_i = e_tmp[!is.na(start),]
    start = start[!is.na(start)]
    len = nchar(rownames(e_i))

    if (nrow(e_i) == 0)
    {
      print(paste0("No sncRNA species from the query matrix were mapped to the parental RNA '", names(parental_rna_list[i]), "' allowing ",  max.mismatch, " mismatch."))
      next
    }

    if (mean(colSums(e_i)) < min_count)
    {
      print(paste0("The mean total read-count per sample for reads mapped on ", names(parental_rna_list[i]), " < ", min_count, ". Try with less-stringent min_count threshold."))
      next
    }

    coverage = c()
    if (nrow(e_i) > 0) {
      for (j in 1:ncol(e_i))
      {
        coverage = cbind(coverage, calc.coverage(e_i[,j], start, len, nchar(parental_rna_list[i])))
      }
      colnames(coverage) = colnames(e_i)
      rownames(coverage) = 1:nchar(parental_rna_list[i])

      coverage[is.nan(coverage)] <- 0

      ### Check if all rows in coverage have the same value within each row
      all_rows_constant <- all(apply(coverage, 1, function(row) length(unique(row)) == 1))

      ### Proceed only if NOT all rows are constant
      if (!all_rows_constant) {

        if (is.null(shuffling_round)) {
          if (method == 1) {
            shuffling_round <- getOption("shuffling_round", 100)
          } else if (method == 2) {
            shuffling_round <- 1000
          } else {
            return(message("Error: method must be either 1 or 2"))
          }
        }

        if (method==2) tmp_result = diff_frag.model2(coverage, cl, shuffling_round = shuffling_round) ### running Model 2 (method 2) to identify differential fragmentation
        else tmp_result = diff_frag.model1(coverage, cl, shuffling_round = shuffling_round) ### else Model 1 (method 1) will be run, by default
        mean_diff[i] = tmp_result[1]
        p[i] = tmp_result[2]
        species_num[i] = nrow(e_i)
      }
    }
  }

  result = data.frame(names(parental_rna_list), mean_diff, p, species_num)
  colnames(result) = c("parental_rna", "mean_diff", "p", "species_num")
  result <- result[!is.na(mean_diff),]
  if (nrow(result) == 0){
    print(paste0("No sncRNA species sequences from the query matrix were mapped on parental RNA '", sncrna_family, "' with the chosen thresholds (mismatch <=", max.mismatch, ", min_count >=", min_count, ".)"))
    return(NULL) }

  ### save to file
  if (identical(output, "temp")) {
    output <- tempfile(fileext = ".txt")
    message("No output filename provided. Saving to temporary file: ", output)
  }else if (is.null(output)) {
    output <- file.path(getwd(), paste0("qMAP_", method, "_", sncrna_family, "_mismatch_", max.mismatch,  "_shuffling_round_", shuffling_round, "_output.txt"))
    message("No output filename provided. Saving to: ", output)
  } else {
    # Normalize path to absolute to avoid confusion
    output <- normalizePath(output, mustWork = FALSE)
  }
  # Try writing the file
  tryCatch({
  write.table(result, file = output, sep = "\t", row.names = FALSE, quote = FALSE)
  message("Result saved to: ", output)
  }, error = function(e) {
    warning("Failed to save output file. Error: ", conditionMessage(e))
  })
  cat("Files written to:", dirname(output), "\n\n")
  return(result)
}

########## qMAP function ends ##########


#' Run qMAP.MH on sncRNA matrix file (file-based input)
#'
#' This is a variant for 'qMAP.MH.single_family( )' that reads the input matrix from a file.
#'
#' @param input_file Path to a tab-delimited read-count matrix. The first two columns should be "Sequence" and "Annotation", and the rest are counts.
#' @param ... Additional arguments passed to `qMAP.MH.single_family()`: `cl`, `sncrna_family`, `parental_rna_filename`, etc.
#'
#' @return A data frame with results from the qMAP.MH module
#' @export
qMAP.MH.single_family.1 <- function(input_file, ...) {
  qMAP.MH.single_family(function() read.delim(input_file), ...)
}
#' Run qMAP.MH on sncRNA data frame (in-memory)
#'
#' This is a variant for 'qMAP.MH.single_family( )' where the matrix is already in memory as a data frame.
#'
#' @param input_dataframe A data frame with "Sequence", "Annotation", and read-count columns.
#' @param ... Additional arguments passed to `qMAP.MH.single_family()`: `cl`, `sncrna_family`, `parental_rna_filename`, etc.
#'
#' @return A data frame with results from the qMAP.MH module
#' @export
qMAP.MH.single_family.2 <- function(input_dataframe, ...) {
  qMAP.MH.single_family(function() input_dataframe, ...)
}

#' qMAP.MH.single_family - to identify individual sncRNA species that contribute to the observed differential fragmentation of a given parental RNA across groups
#'
#' @param matrix_file A user-provided read-count matrix where the first and second columns are labeled as "Sequence" and "Annotation", respectively. The remaining columns contain read counts from different samples under study.
#' @param cl A sample classification vector allowing only two groups: 1 for controls, 2 for cases, and -1 for the samples to be excluded from the analysis.
#' @param parental_rna_filename A FASTA file containing sequences of parental RNA families description.
#' @param sncrna_family This parameter specifies a single sncRNA family name to search for in the user-provided parental RNA FASTA file. If a match is found, differential fragmentation analysis will be performed for that family. Use any single family name, such as : “mature-tRNA-Ala-AGC”, “tRNA-Ala-AGC”, “mature-mt_tRNA-His-GTG”, “mt_tRNA-His-GTG”, “5S-rRNA”, “16S-rRNA”
#' @param max.mismatch This parameter specifies the maximum number of allowed mismatches (default is 1) when mapping sequence reads onto the FASTA sequence of the parental RNA.
#' @param minimum_mean_count This parameter sets the minimum allowed mean read count (default: 5) across samples for the individual sncRNA species. Only sncRNA species (i.e., rows in the read-count matrix)  with mean read counts greater than or equal to minimum_mean_count will be included in the analysis.
#' @param output  This parameter specifies the name of the output file, By default, if no name is provided,  the output file will be saved in the working directory with a name like “qMAP_MH_output_X.txt", where X represents query sncrna_family parameter used in the run. If '"temp"', the function creates a temporary output file (used for examples/testing only).
#' @return It will return, for each queried parental RNA (sncrna_family), an output data-frame with sequence (the sncRNA species sequence), or (the odds ratio), p (the P-value inferred from the Mantel-Haenszel statistic), and start (the start position on the parental RNA where the sncRNA species is mapped).
#' @note In example runs, `output = "temp"` is used to avoid writing files during checks. In real analyses, you may provide a file path or let the function auto-name the output (see manual).
#'
#' @examples
#' # Note:  After installation, one can find the example files in "../qMAPModules/extdata/"
#' # To run on a specific input, provide the full path as: qMAPModules(a = "/path/to/your_matrix.txt")
#' sample_matrix_1 <- system.file("extdata", "sample_matrix1.txt", package = "qMAPModules")
#' sample_matrix_2 <- system.file("extdata", "sample_matrix2.txt", package = "qMAPModules")
#' sample_matrix_3 <- system.file("extdata", "sample_matrix3.txt", package = "qMAPModules")
#' mouse_gtsrna_list <- system.file("extdata", "mmu_gtsrna.txt", package = "qMAPModules")
#' mouse_rRNAs <- system.file("extdata", "mouse_rRNAs.fa", package = "qMAPModules")
#' mouse_tRNAs <- system.file("extdata", "mm10-tRNAs_CCA.fa", package = "qMAPModules")
#' # Run the identification of sncRNA species contributing to the differential fragmentation of 5S-rRNA
#' # family using corresponding FASTA sequence in the user provided parental rna file (mouse_rRNAs.fa)
#' # by mapping sequence reads in the query read-count matrix file (sample_matrix1.txt) comprising read
#' # counts for 5 control (cl=1) and 4 case (cl=2) samples, with allowed maximum mismatch 1 (default)
#' # and considering only the reads with the row mean read-counts across the samples > 5 (default).
#' qMAP.MH.single_family.1(input_file = sample_matrix_1, cl=c(1,1,1,1,1,2,2,2,2),
#'                     sncrna_family="5S-rRNA", parental_rna_filename= mouse_rRNAs, output = "temp")
#'
#' # Run the identification of sncRNA species contributing to the differential fragmentation of
#' # tRNA-Gly-CCC family using corresponding FASTA sequence in the user provided parental rna file
#' # (mm10-tRNAs_CCA.fa) by mapping sequence reads in the query read-count matrix file
#' # (sample_matrix1.txt) comprising read counts for 3 control (first three sample columns in the
#' # matrix) (cl = 1) and 3 case (last three sample columns in the matrix) (cl = 2) samples, with
#' # allowed maximum mismatch 1 (default) and considering only the reads with the row mean read-counts
#' # across the samples > 10. Here, 3 samples are excluded (cl = -1).
#' qMAP.MH.single_family.1(sample_matrix_2, cl=c(1,1,1,-1,-1,-1,2,2,2), sncrna_family="tRNA-Gly-CCC",
#'                   minimum_mean_count  = 10,  parental_rna_filename=mouse_tRNAs, output = "temp")
#'
#' # Run the identification of sncRNA species contributing to the differential fragmentation of
#' # tRNA-Val-AAC family using corresponding FASTA sequence in the user provided parental rna file
#' # (mm10-tRNAs_CCA.fa) by mapping sequence reads in the query read-count matrix file
#' # (sample_matrix3.txt) comprising read counts for 4 Contioin1 (cl = 1), 4 Condition2 (cl =2)
#' # samples, and rest 8 for other conditions, with allowed maximum mismatch 0 and considering only
#' # the reads with the row mean read-counts across the samples > 10. The result output will be saved
#' # to the specified file.
#' qMAP.MH.single_family.1(sample_matrix_3, cl=c(1,1,1,1,2,2,2,2,-1,-1,-1,-1,-1,-1,-1,-1),
#'                   sncrna_family="tRNA-Val-AAC", minimum_mean_count  = 10, max.mismatch = 0,
#'                   parental_rna_filename = mouse_tRNAs, output = "temp")
#'
#' # Running with query read-count matrix as a data frame as input.
#' # First read and save the read-count matrix file as a dataframe:
#' a <- read.delim(sample_matrix_1, header = TRUE)
#' # Then call qMAP_mh function to run the identification of sncRNA species contributing to the
#' # differential fragmentation of tRNA-His-GTG family using corresponding FASTA sequence in the user
#' # provided parental rna file (mm10-tRNAs_CCA.fa) by mapping sequence reads in the query read-count
#' # matrix file (sample_matrix1.txt) comprising read counts for 5 control (cl = 1) and 4 case (cl =2)
#' # samples, with allowed maximum mismatch 1 (default) and considering only the reads with the row
#' # mean read-counts across the samples > 5 (default).
#' qMAP.MH.single_family.2(a, cl=c(1,1,1,1,1,2,2,2,2), sncrna_family="tRNA-His-GTG",
#'			 	parental_rna_filename= mouse_tRNAs, output ="temp")
#'
#' # Run the identification of sncRNA species contributing to the differential fragmentation of each of
#' # the RNA families (listed in a file; here it is mmu_gtsrna.txt) using corresponding FASTA sequences
#' # in the user provided parental rna file (mm10-tRNAs_CCA.fa) by mapping sequence reads in the query
#' # read-count matrix file (sample_matrix2.txt) comprising read counts for 3 control (first three
#' # sample columns in the matrix) (cl = 1) and 3 case (next three sample columns in the matrix) (cl=2)
#' # samples. Here, last 3 samples excluded (cl = -1). It will run with allowed maximum mismatch 1
#' # (default) and considering only the reads with the row mean read-counts across the samples > 5
#' # (default). The result output will be saved to the specified file, qmap.mh_multifamily_results.txt.
#' # To run qMAP_mh for multiple families, it is suggested to run with second variant
#' # (qMAP.MH.single_family.2) in a loop as follow:
#' # First read and save the read-count matrix file as a dataframe:
#' a = read.delim(sample_matrix_2)
#' cl = c(1,1,1,2,2,2,-1,-1,-1)
#' result = c()
#' # Read sncRNA families from list in the file:
#' sncrna_family_list = read.table(mouse_gtsrna_list)[,1]
#' # Then for each family in the list, call qMAP_mh function
#' for (i in 1:length(sncrna_family_list)) {
#'     result = rbind(result, qMAP.MH.single_family.2(a, cl=cl, sncrna_family=sncrna_family_list[i],
#'     parental_rna_filename = mouse_tRNAs, output = "temp"))
#'   }
#' # Creating a safe temporary output file exclusively for run examples (see manual)
#' out_file <- tempfile("qmap.mh_multifamily_results", fileext = ".txt")
#' write.table(result, file = out_file, quote=FALSE, row.names=FALSE, sep="\t")
#'
#' @export

### function to perform Mantel-Haenszel test for a single sncRNA family
### cl: sample classification vector: only two groups are allowed (1 for controls, 2 for cases, and -1 for the samples to be excluded)
qMAP.MH.single_family <- function(matrix_file, cl, sncrna_family, parental_rna_filename, max.mismatch = getOption("max.mismatch", 1), minimum_mean_count = getOption("minimum_mean_count", 5), output = NULL)
{
  #### qmap function body starts ####
  args = commandArgs(trailingOnly=TRUE)

  ########## Part1: defining different functions to call ##########

  options(warn=1)

  # Function to check if all values in the matrix are whole numbers
  check_whole_numbers <- function(exp) {

    # Check for non-numeric values
    if (!all(sapply(exp, is.numeric))) {
      stop("All columns in the read-count matrix must be numeric. There are alphabets/characters in the matrix.")
    }
    ### Check for non-finite
    if (!all(sapply(exp, is.finite))) {
      stop("The read-count matrix contains non-finite values (NA, Inf, or NaN).")
    }

    ### Check for non-integer values
    if (!all(sapply(exp, is.integer))) {
      message("Warning: There are some non-integer values in the read-count matrix (e.g., 0.5, 1.01). Only whole numbers should be used.")
    }
    return(TRUE)
  }

  ### Function to call Mantel-Haenszel test
  calc.MH <- function(rna_count, cl)
  {
    total = colSums(rna_count)
    or = rep(NA, nrow(rna_count))
    p = rep(NA, nrow(rna_count))
    for (i in 1:nrow(rna_count))
    {
      count_matrix = rbind(rna_count[i,], total - rna_count[i,])
      count_matrix1 = count_matrix[,cl==1]
      count_matrix2 = count_matrix[,cl==2]

      data = c()
      for (j in 1:ncol(count_matrix1))
      {
        for (k in 1:ncol(count_matrix2))
        {
          data = append(data, c(count_matrix2[,k], count_matrix1[,j]))
        }
      }
      mh_data = array(data, dim=c(2, 2, length(data)/4))
      tmp = mantelhaen.test(mh_data)
      or[i] = tmp$estimate
      p[i] = tmp$p.value
    }
    result = data.frame(rownames(rna_count), or, p)
    return(result)
  }

########## Part2: Processing of input matrix and different parameters to be used in the function ##########

    ### trim sncRNA family name and look for the corresponding fasta sequences in the user provided parental rna file
  sncrna_family = sub("mature-", "", sncrna_family)
  s = unlist(read.fasta(parental_rna_filename, seqtype="DNA", as.string=T, forceDNAtolower=F))
  if (grepl("rRNA", sncrna_family)) s = s[names(s) == sncrna_family]
  else s = s[grepl(sncrna_family, names(s))]
  uni_s = unique(s)
  z = match(uni_s, s)
  parental_rna_list = s[z]

  if (length(s) == 0)
  {
    print(paste0("The specified sncRNA family - ", sncrna_family,  ", does not exist in the provided parental RNA file."))
    return(NULL)
  }

  ### Processing of input matrix and different parameters to be used in the function ###
  ### a: a user-provided read-count matrix with first and second column as "Sequence" and "Annotation", respectively. Rest of the columns having read counts from different samples under study. Sequence must be unique.###
  a = matrix_file()
  anno = a[,1:2]  ### Separate the columns with unique Sequence and its Annotation information
  e = a[,3:ncol(a)]  ### Separate the columns with read-count values

  check_whole_numbers(e)

  ### Confirming the sample size in cl option matched with that in the input matrix file
  cl_sum = length(cl[cl==1]) + length(cl[cl==2]) + length(cl[cl==-1])
  if (ncol(e) != cl_sum)
  {
    print(paste0("Error: Number of samples (",   ncol(e), ") in the input matrix does not match the number of samples specified in the cl option (", cl_sum, ")."))
    return(NULL)
  }

  rownames(e) = a[,1]
  e = e[, cl==1 | cl==2]
  cl = cl[cl==1 | cl==2]
  if (length(cl[cl==1]) < 2 | length(cl[cl==2]) < 2)  ### Confirming the sample size of both group is not less than two
  {
    print("Error: sample size of at least one group is less than two.")
    return(NULL)
  }

  ### exclude the sncRNA species (rows in the matrix) with zero read-counts
  anno = anno[rowMeans(e) > 0,]
  e = e[rowMeans(e) > 0,]


  ### Processing the input matrix for query sncRNA family
  if (grepl("rRNA", sncrna_family))
  {
    is.sncrna_family = (anno$Annotation==sncrna_family)
  } else
  {
    if (grepl("tRNA", sncrna_family)) is.sncrna_family = grepl(paste("mature-", sncrna_family, sep=""), anno$Annotation)
    else is.sncrna_family = grepl(sncrna_family, anno$Annotation)
  }
  e_tmp = e[is.sncrna_family,]   ### Extracting read-count matrix values for the species with query sncRNA family annotation

  if (nrow(e_tmp) == 0)
  {
    print(paste0("There are no read fragments in the input matrix file annotated as - ", sncrna_family  , "."))
    return(NULL)
  }

  ### exclude the sncRNA species (rows in the input matrix) with mean count less than minimum_mean_count
  anno = anno[rowMeans(e_tmp) >= minimum_mean_count,]
  e_tmp = e_tmp[rowMeans(e_tmp) >= minimum_mean_count,]

  if (nrow(e_tmp) < 2 )
  {
    print(paste0("There are fewer than 2 sncRNA species annotated as - ", sncrna_family, ", and having mean count >= ", minimum_mean_count, "."))
    return(NULL)
  }

  ### add one to each values in the input matrix to avoid any mathematical calculation involving zeros
  e_tmp <- e_tmp + 1

  ### Find the match position of the query sequence within the parent rna
  ### Function to find the match position of the small sequence within the long sequence
  find_position_with_mismatches <- function(long_seq, small_seq) {
    match_result <- matchPattern(small_seq, long_seq, max.mismatch = max.mismatch)
    ### If there are matches, return the first match position (1-based index)
    if (length(match_result) > 0) {
      return(start(match_result)[1])  ### Return the first position
    } else {
      return(NA)  ### If no match is found
    }
  }

  ### Running function find_position_with_mismatches over the length of parent rna
  result = c()
  for (i in 1:length(parental_rna_list))
  {
    start = mapply(find_position_with_mismatches, parental_rna_list[i], rownames(e_tmp))
    e_i = e_tmp[!is.na(start),]
    start = start[!is.na(start)]

    if (nrow(e_i) < 2 )
    {
      print(paste0("Found a total of ", nrow(e_tmp), " species with annotation as ", sncrna_family, " and mean count >= ", minimum_mean_count, ". However, ", nrow(e_i), " species mapped to the parental RNA ", names(parental_rna_list[i]), " with the maximum allowed mismatch of ",  max.mismatch, ". Atleast two species are needed to call Mantel-Haenszel test"))
      next
    }
    tmp_result = data.frame(names(parental_rna_list[i]), calc.MH(e_i, cl), start, row.names = NULL )
    colnames(tmp_result) = c("parental_rna", "sequence", "or", "p", "start")

    result = rbind(result, tmp_result)
    ### Sort the dataframe first by `parental_rna`, then by `matched start position` and then by 'sequence'
    result <- result %>% arrange(.data$parental_rna, .data$start, .data$sequence)
  }

  if (length(result) == 0 ){
    print(paste0("Query sncRNA family (", sncrna_family, ") was found to have less than two sncRNA species mapped with the chosen thresholds (mismatch <=", max.mismatch, ", minimum_mean_count >=", minimum_mean_count, ".)"))
    return(NULL) }

### save to file
if (identical(output, "temp")) {
  output <- tempfile(fileext = ".txt")
  message("No output filename provided. Saving to temporary file: ", output)
}else if (is.null(output)) {
  output <- file.path(getwd(), paste0("qMAP_MH_output_", sncrna_family, ".txt"))
  message("No output filename provided. Saving to default: ", output)
} else {
  # Normalize path to absolute to avoid confusion
  output <- normalizePath(output, mustWork = FALSE)
}
# Try writing the file
tryCatch({
  write.table(result, file = output, sep = "\t", row.names = FALSE, quote = FALSE)
  message("Result saved to: ", output)
}, error = function(e) {
  warning("Failed to save output file. Error: ", conditionMessage(e))
})
cat("Files written to:", dirname(output), "\n\n")
return(result)
}

########## qMAP_MH function ends ##########
