
[User Manual](inst/docs/qMAP_v1.0.0_User_Manual.pdf)

# qMAP (quantitative mapping for RNA differential fragmentation) Version 1.0.0

## Install instructions:

#### To install qMAP directly from GitHub, run the following command in R:

      install.packages("remotes")

      remotes::install_github("cozyrna/qMAP")

#### Using source code

    1. Download the package source file qMAP_*.tar.gz
    2. Install it by opening R or RStudio and running:

    install.packages("~/qMAP_*.tar.gz", repos = NULL, type = "source")

(Replace ~/qMAP_*.tar.gz with the actual path and filename of the downloaded package.)

## Running qMAP modules:

#### There are two modules in this package : 

#### qMAP module - for the identification of differentially fragmented parental RNAs between groups. The qMAP module have two models or methods:
 	
 	    qMAP Model 1 identify differentially fragmented parental RNAs by comparing the t-statistic values from the real sample grouping against absolute values of the t-statistic based on the permutated grouping by shuffling the samples between groups.

	    qMAP Model 2 identify differentially fragmented parental RNAs by comparing the  overall difference in normalized coverage between the two groups against the overall difference in normalized coverage based on the permutated grouping by shuffling the samples between groups.

#### qMAP_mh module - for the identification of the sncRNA species that are involved in the observed differential fragmentation in a parental RNA by using Mantel-Haenszel procedure with continuity correction, and then calculating the Mantel-Haenszel statistic on a chi-square distribution.

#### After installing the qMAP package, load it in R with:

    library(qMAP)

To access the description and help documentation for the package and its functions, use:

    ?qMAP.single_family

    ?qMAP.MH.single_family


There are two variants of each of these two modules for a single sncRNA family:

1. Function call with read-count matrix file path as input:
   
      	qMAP.single_family.1(input_file, ...)
   
      	qMAP.MH.single_family.1(input_file, ...) 

2. Function call with read-count matrix as a data frame as input:      
 
      	qMAP.single_family.2(input_dataframe, ...)
   
      	qMAP.MH.single_family.2(input_dataframe, ...) 
                                                                      
The second variant is useful when one has to call qMAP or qMAP_mh function for multiple sncRNA families. 

#### Note: After installation, the sample matrix files used in the example runs can be found in the "../qMAP/extdata/" directory. 
#### To run the examples exactly as shown in the help documentation, you must first set the working directory to the base folder of the installed qMAP package. For example: setwd("./..../R/.../qMAP/"). Alternatively, you can directly provide the full path to the appropriate example matrix files without changing the working directory.

#### Here are the sample runs with example data sets:

    qMAP.single_family.1("./extdata/sample_matrix1.txt", cl=c(1,1,1,1,1,2,2,2,2), sncrna_family="5S-rRNA", parental_rna_filename="./extdata/mouse_rRNAs.fa")

    a <- read.delim("./extdata/sample_matrix1.txt", header = TRUE)
    qMAP.single_family.2(a, cl=c(1,1,1,1,1,2,2,2,2), sncrna_family="tRNA-His-GTG", parental_rna_filename="./extdta/mm10-tRNAs_CCA.fa", max.mismatch = 1, method = "2")


    qMAP.MH.single_family.1("./extdata/sample_matrix1.txt", cl=c(1,1,1,1,1,2,2,2,2), sncrna_family="5S-rRNA", parental_rna_filename="./extdata/mouse_rRNAs.fa")

    a <- read.delim("./extdata/sample_matrix1.txt", header = TRUE)
    qMAP.MH.single_family.2(a, cl=c(1,1,1,1,1,2,2,2,2), sncrna_family="tRNA-His-GTG", parental_rna_filename="./extdata/mm10-tRNAs_CCA.fa")

To run for multiple families, see example runs mentioned in the user manual.

