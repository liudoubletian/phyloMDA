#' @title Fit the Zero-inflated Dirichlet-tree multinomial (ZIDTM) distribution.
#'
#' @param otu.tab	 a data frame or matrix containing the count data. Rows of the matrix represent observations and columns are the taxa.
#' @param tree	 the phylogenetic tree (should be binary here).
#'
#' @return Returns a list, containing results of each subtree.
#' @references Zhou, C., Zhao, H., and Wang, T.: Transformation and differential abundance analysis of microbiome data incorporating phylogeny. Bioinformatics btab543, https://doi.org/10.1093/bioinformatics/btab543 (2021).
#' @examples
#' library(phyloseq)
#' data("combo.phyloseq.obj")
#' otu.tab <- t(otu_table(combo.phyloseq.obj))
#' tree <- phy_tree(combo.phyloseq.obj)
#' fit <- ZIdtmFit(otu.tab, tree)
#' Extract_logL(fit)
#' @export
#'
#' @import phyloseq adaANCOM methods
#'

ZIdtmFit <- function (otu.tab, tree) {
  ytree <- YtreeFun(otu.tab, tree)
  #setClass('res', representation(alpha="numeric",pi="numeric", logL='numeric'))
  fit <- lapply(ytree, function(x) {
     est.zidm.EM(x)
   # new('res', alpha=tmp$alpha, pi=tmp$pi, logL=tmp$loglik)
    })
  fit <- list("zidtm",fit)
  fit
}


