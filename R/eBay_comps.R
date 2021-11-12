#' An empirical Bayes approach for transforming raw counts into relative abundances
#'
#' @param otu_tab a data frame or matrix containing the count data. Rows of the matrix represent observations and columns are the taxa.
#' @param prior Dirichlet, Dirichlet-tree, or zero-inflated Dirichlet-tree prior for the proportions.
#' @param tree the phylogenetic tree.
#' @param model smooth methods for zero-inflated Dirichlet-tree mulinomial, including 0.5, DM, ZIDM, or MIX.
#' @return Returns the estimated microbial compositions by empirical Bayes.
#' @references Liu, T., Zhao, H., Wang, T.: An empirical Bayes approach to normalization and differential abundance testing for microbiome data. BMC Bioinformatics 21, 225 (2020).
#' @references Zhou, C., Zhao, H., and Wang, T.: Transformation and differential abundance analysis of microbiome data incorporating phylogeny. Bioinformatics btab543, https://doi.org/10.1093/bioinformatics/btab543 (2021).
#' @examples
#' library(phyloseq)
#' otu_tab <- t(combo.phyloseq.obj@otu_table@.Data)
#' tree <- phy_tree(combo.phyloseq.obj)
#' eBay.comps <- eBay_comps(otu_tab, prior = "Dirichlet", tree=NULL)
#' eBay.tree.comps <- eBay_comps(otu_tab, prior = "Dirichlet-tree", tree = tree)
#' eBay.zitree.comps <- eBay_comps(otu_tab, prior = "zero-inflated-Dirichlet-tree", tree = tree, model= "MIX")
#' @export
#' @importFrom MGLM  MGLMfit
#' @importFrom caper clade.matrix
#' @import phyloseq
#' @import magrittr
#' @import plyr
#' @import adaANCOM


eBay_comps <- function(otu_tab, prior, tree=NULL, model= "MIX"){
  if(prior == "zero-inflated-Dirichlet-tree"){
    eBay.comps <- PosteriorMeanTrans(otu_tab, tree, model = model)
    eBay.comps <- eBay.comps$probability$L
    colnames(eBay.comps) <- tree$tip.label
    rownames(eBay.comps) <- rownames(otu_tab)
  }
  if(prior == "Dirichlet"){
  fit_dm <- MGLMfit(data=otu_tab, dist='DM')
  alpha <- fit_dm@estimate
  eBay.comps <- t(apply(otu_tab,1,function(x){(x+alpha)/(sum(x+alpha))}))
  }
  if(prior == "Dirichlet-tree"){
    if(is.null(tree)) {
    stop("The tree structure must be provided")
    }
  inter_node <- unique(tree$edge[,1]) ##internal nodes of the tree
  leaf <- c()
  for (i in 1:length(inter_node)){
    leaf_node <- tree$edge[which(tree$edge[,1]==inter_node[i]),2]
    leaf <- append(leaf,leaf_node)
  }####the children node of each internal node

  ytree <- YtreeFun(otu_tab, tree)
  fit_dtm <- lapply(ytree, function(x) MGLMfit(x, dist = 'DM'))

  node.comps <- lapply(1:length(fit_dtm), function(x){
    alpha <- fit_dtm[[x]]@estimate
    data.frame((apply(ytree[[x]],1,function(j){(j+alpha)/(sum(j+alpha))})))
  }) %>%  do.call(rbind.fill,.) %>% t
  colnames(node.comps) <- leaf

  clade.mat <- clade.matrix(tree)$clade.matrix
  leaf.path <- apply(clade.mat,2,function(x) { which(x==1)} ) ###the path of each leaf node on the tree

  eBay.comps <- lapply(leaf.path,function(x){
    x <- x[-which(x==(length(tree$tip.label)+1))]
    sub.mat <- subset(node.comps, select=as.character(x))
    as.data.frame(t(apply(sub.mat,1,prod)))
  })

  eBay.comps <- do.call(rbind.fill,eBay.comps) %>% t
  colnames(eBay.comps) <- tree$tip.label
  rownames(eBay.comps) <- rownames(otu_tab)
  }

  return(eBay.comps)
}


