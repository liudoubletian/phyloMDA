#' @title Visualization on the tree
#'
#' @description The results of TASSO and TreeFusedlasso vasualized on the tree.
#'
#' @param results results from \link{TASSO} or \link{TreeFusedlasso}.
#' @param tree the phylogenetic tree.
#' @param model indicate the model used.
#' @param IC IC = BIC for Bayesian information criterion and IC = AIC for Akaike's information criterion.
#' @param tree.type a character string specifying the type of phylogeny to be drawn; it must be one of "phylogram" (the default), "cladogram", "fan", "unrooted", "radial" or any unambiguous abbreviation of these.
#' @param color the colours used for the tip labels or choosed pattern.
#'
#' @return A ggplot object.
#'
#' @examples
#' library(phyloseq)
#' otu.tab <- t(otu_table(combo.phyloseq.obj))
#' eBay.comps <- eBay_comps(otu.tab, prior = "Dirichlet")
#' tree <- phy_tree(combo.phyloseq.obj)
#' metadata <- sample_data(combo.phyloseq.obj)
#' y <- metadata$bmi
#' fit_tasso <- TASSO(y, eBay.comps, tree)
#' plot_onTree(fit_tasso, tree)
#'
#' fit_tflasso1 <- TreeFusedlasso(y, eBay.comps, tree)
#' plot_onTree(fit_tflasso1, tree, model='TreeFusedlasso')
#' @export
#' @import ggtree ggplot2
#'


plot_onTree <- function(results, tree, model='TASSO', IC='BIC', tree.type='fan', color='darkred') {
  tmp <- ICgenlasso(results, IC=IC)
  b1 <- results$beta[, which.min(tmp[,2])]
  b2 <- round(b1 - mean(b1) , 8)
  all_taxa <- tree$tip.label
  K <- length(all_taxa)
  if(model=='TASSO') {
    nonzreo <- all_taxa[b2!=0]
    cls <- list(nonzero=nonzreo, zero=setdiff(all_taxa,nonzreo))
    tree_cls <- groupOTU(tree, cls)
    g1 <- ggtree(tree_cls, aes(color=group), layout=tree.type, size=0.3, ladderize=F) +
      geom_tiplab(size=2, aes(color=group), hjust=-0.05, align = F) + theme_tree() +
      guides(color=guide_legend(title="Estimated coefficients"))
    res <- g1 + scale_color_manual(values=c(color,'black'))
  }
  GetFusedPattern <- function(beta, tree) {
    edge <- tree$edge
    all_taxa <- tree$tip.label
    K <- length(all_taxa)

    cls <- list()
    cls_names <- c()
    g <- i <- 1
    groups <- 1:tree$Nnode

    find_descent <- function(edge, v) {
      # find all the descents
      childs <- c(); par <- v
      while(length(par) >0){
        ch <- unlist(sapply(par, function(x) edge[edge[,1]==x,2]))
        par <- ch[ch>min(edge[,1])]
        childs <- c(childs, ch)
      }
      childs
    }
    while(length(groups) >0 ) {
      ds <- find_descent(edge, i+K)
      ds_l <- ds[ds<=K]
      ds_f <- ds[ds>K]
      bi <- beta[ds_l]
      groups <- setdiff(groups, i)
      if(length(unique(bi))==1){
        cls[[g]] <- all_taxa[ds_l]
        # cls_names <- c(cls_names, unique(bi))
        cls_names <- c(cls_names, i+K)
        groups <- setdiff(groups, ds_f-K)
        g <- g+1
        i <- i + length(ds_f)
      }
      i <- i+1
    }
    # cls[[length(cls)+1]] <- setdiff(all_taxa, unlist(cls))
    # names(cls) <- c(round(cls_names, 2), 'other') # sprintf("%0.3f", cls_names)
    names(cls) <- cls_names
    cls
  }

  if(model=='TreeFusedlasso') {
    cls <- GetFusedPattern(b2, tree)
    g1 <- ggtree(tree, ladderize=F, layout=tree.type, size=0.3) + theme_tree() +
      geom_tiplab(size=2, hjust=-0.05, align = F)
    cols <- rep('black', tree$Nnode)
    cols[as.numeric(names(cls))-K] <- color
    siz <- rep(0.1, 61)
    siz[as.numeric(names(cls))-K] <- sapply(cls, length)
    res <- g1 +  geom_point2(aes(subset=!isTip), color=cols, alpha=0.8, size=siz/2)
  }
  return(res)
}

