#' Fit the Zero-inflated Dirichlet-tree multinomial (ZIDTM) regression.
#'
#' @param otu.tab a data frame or matrix containing the count data. Rows of the matrix represent observations and columns are the taxa.
#' @param X.mean a data frame or matrix containing the covariates link to the mean parameter. Rows of the matrix represent observations.
#' @param X.disp a data frame or matrix containing the covariates link to the dispersion parameter. Rows of the matrix represent observations.
#' @param X.zero a data frame or matrix containing the covariates link to the zero prevalance. Rows of the matrix represent observations.
#' @param tree the phylogenetic tree
#' @param test.type f test.type = "Mean", the function will test for differential mean (Default). If test.type = "Disp", the function will test for differential dispersion. If test.type = "Freq", the function will test for differential presence-absence frequency.
#' @param index index position for the interest variable for testing and remain treat as confounders.
#' @return Returns a list, containing testing results of each subtree.
#' @references Zhou, C., Zhao, H., and Wang, T.: Transformation and differential abundance analysis of microbiome data incorporating phylogeny. Bioinformatics btab543, https://doi.org/10.1093/bioinformatics/btab543 (2021).
#' @examples
#' library(phyloseq)
#' data("combo.phyloseq.obj")
#' otu.tab <- t(otu_table(combo.phyloseq.obj))
#' tree <- phy_tree(combo.phyloseq.obj)
#' sodium <- sample_data(combo.phyloseq.obj)$sodium
#' fit <- ZIdtmReg(otu.tab, sodium, tree=tree)
#' @export
#' @import phyloseq miLineage



ZIdtmReg <- function(otu.tab, X.mean, X.disp=NULL, X.zero=NULL, tree, test.type='Mean', index=1) {
  ytree <- YtreeFun(otu.tab, tree)
  fit <- lapply(ytree, function(y)  {
    ZIGDM(y, X4freq=X.zero, X4mean=X.mean, X4disp=X.disp, test.type, X.index=index)
    })
  fit
}



