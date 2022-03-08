#' function to create ggtree with colors from MCC tree

#' @param mcc - mcc tree from BEAST
#' @param gps - mcc tree from BEAST


ggtree_mcc <- function(mcc,gps) {
  
  #mcc = treeio::read.beast("data/batRABV.mcc.tree")
  
  
  colnames(gps) = c("location", "lat", "lon")
  locations = gps$location
  n = length(locations)
  mrst = 2005.5
  
  
  # 5.Plot MCC tree using ggtree---------------------------------------
  # Location colors 
  color_states = c(brewer.pal(n = 8, name = 'Dark2'), brewer.pal(n = 12, name = 'Paired'))[1:n]
  names(color_states) = locations
  rootNode = unique(mcc@phylo$edge[!mcc@phylo$edge[,1] %in% mcc@phylo$edge[,2], 1])
  rootCol = color_states[mcc@data$state[mcc@data$node == rootNode]]
  
  # Piechart of root location
  mcc@data$state.prob = as.numeric(mcc@data$state.prob) 
  rootLocationProbs = 
    
    ancestorlocs = mapply(
      function(x,y) {
        nI = length(x)
        nT = length(color_states)
        x = as.numeric(x)
        names(x) = y
        if (nI < nT) {
          x = c(x, rep(0, nT-nI))
          names(x)[(nI+1):nT] = names(color_states)[!names(color_states) %in% y]
        }
        
        as.data.frame(t(x))
      },
      x=mcc@data$state.set.prob,
      y=mcc@data$state.set,
      SIMPLIFY=FALSE
    )
  ancestorlocs = do.call("rbind", ancestorlocs) 
  ancestorlocs$node = mcc@data$node
  rootloc = 
    
    ancestorlocs[ancestorlocs$node != mcc@phylo$Nnode+2, names(color_states)] = NA 
  ancestorlocs = ggtree::nodepie(ancestorlocs, cols = 1:6, color = color_states)
  
  mcc@data$posterior = as.numeric(mcc@data$posterior)
  mcc@data$state.prob = as.numeric(mcc@data$state.prob)
  
  final_tree <- ggtree::ggtree(mcc, aes(color = state, size = state.prob), mrsd=date_decimal(mrst), ladderize = T, right=T) +
    ggtree::geom_rootedge(rootedge = 3, colour = rootCol) +
    ggtree::geom_tippoint(size = 2) +
    ggtree::geom_nodepoint(size = 2) +
    ggtree::geom_tippoint(aes(color = state), size = 1) +
    # geom_inset(ancestorlocs, width=.1, height=.1, vjust=-5) +
    ggtree::scale_color_manual(values = color_states) +
    scale_size_continuous(range = c(0.2, 1)) +
    ggtree::theme_tree2() +
    labs(color = "State", size = "State probability")
  
 
 return(final_tree)
}