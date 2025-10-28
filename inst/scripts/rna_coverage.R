library(Biostrings)

rna.coverage <- function(matrix_file, sncrna_family, parental_rna, max.mismatch = 1)
{
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

  sncrna_family = sub("mature-", "", sncrna_family)

  ### a: a user-provided read-count matrix with first and second column as "Sequence" and "Annotation", respectively. Rest of the columns having read counts from different samples under study. Sequence must be unique.###
  a = read.delim(matrix_file)
  anno = a[,1:2]
  e = a[,3:ncol(a)]
  rownames(e) = a[,1]

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

  find_position_with_mismatches <- function(long_seq, small_seq) {
    match_result <- matchPattern(small_seq, long_seq, max.mismatch = max.mismatch)
    if (length(match_result) > 0) {
      return(start(match_result)[1])  ### Return the first position
     } else {
      return(NA)  ### If no match is found
    }
  }

  start = mapply(find_position_with_mismatches, parental_rna, rownames(e_tmp))
  e_tmp = e_tmp[!is.na(start),]
  start = start[!is.na(start)]
  len = nchar(rownames(e_tmp))

  if (nrow(e_tmp) == 0)
  {
    print("No sncRNA species from the query matrix were mapped to the parental RNA")
    return(NULL)
  }

  coverage = c()
  for (i in 1:ncol(e_tmp))
  {
    coverage = cbind(coverage, calc.coverage(e_tmp[,i], start, len, nchar(parental_rna)))
  }
  colnames(coverage) = colnames(e_tmp)
  rownames(coverage) = 1:nchar(parental_rna)
  coverage[is.nan(coverage)] <- 0
  return(coverage)
}

matrix_file = "./extdata/sample_matrix1.txt"
sncrna_family = "5S-rRNA"
parental_rna = "GTCTACGGCCATACCACCCTGAACGCGCCCGATCTCGTCTGATCTCGGAAGCTAAGCAGGGTCGGGCCTGGTTAGTACTTGGATGGGAGACCGCCTGGGAATACCGGGTGCTGTAGGCTTT"
coverage = rna.coverage(matrix_file, sncrna_family, parental_rna)
print(coverage)
