#
# Copyright (c) 2015 Ieuan Clay
# based on code from https://github.com/biolion/genbench
# Copyright (c) 2015-2016 BeDataDriven B.V.
# License: http://www.gnu.org/licenses/gpl.html GPL version 2 or higher
#
# ## a minimal script for extracting/cleaning/and visualising entrez RIF networks
# load iGraph network object from falt table
# do some basic analysis/coloring/plotting
# do some subsetting
#  - random 1%
#  - limit genes by GO terms (or anything other annotation)
#  - limit papers by MeSh terms

### set up session
#rm(list=ls())
# reproducibility
set.seed(8008)

## packages
library(igraph) # requires > R-15.2
library(XML)
library(R.utils)
library(plyr)
library(reshape)
library(utils)
library(sqldf)
## global vars
PATH <- file.path(normalizePath("./"), "generifs_basic.gz")

#### functions
do.plot <- function(network, title = "", layout = igraph::layout.kamada.kawai,
                    # node properties
                    shape = "circle", size = 1, colour = "black", label = NA,
                    # edge properties
                    weight = 1
                    ) 
{
  cat("> START: do.plot()\n")

  ### for simplicity, basic plotting of iGraph objects
  # provides default arguments for most parameters
  # otherwise can be passed from:
  # V(..igraph instance..)$..property name..
  # and
  # E(..igraph instance..)$..property name..
  # functions

  plot(
    network,                          # the graph to be plotted
    layout = layout,	                # kamada.kawai|fruchterman.reingold|lgl|auto -  the layout method. see the igraph documentation for details
    main = title,	                    # specifies the title
    vertex.shape = shape,             # nodeshape
    vertex.size = size,               # nodesize
    vertex.color = colour,            # node colour
    edge.width = weight,              # edge weight
    vertex.label.color = 'black',		  # the color of the name labels
    vertex.label.font = 2,			      # the font of the name labels
    vertex.label = label,		          # specifies the lables of the vertices. in this case the 'name' attribute is used
    vertex.label.cex = 1			        # specifies the size of the font of the labels. can also be made to vary

  )



  cat("> END: do.plot()\n")
}

do.load.edges <- function(PATH){

  # unpack data
  #tmpfile <- file.path(dirname(PATH), "rif.tmp")
  #if(file.exists(tmpfile)) file.remove(tmpfile)
  #gunzip(filename = PATH, remove = FALSE, destname = tmpfile)

  ## load and 'query' downloaded data
  # There is more information on sqldf on the sqldf home page:
  # http://sqldf.googlecode.com
  generifs_basic <- read.delim("generifs_basic", header = TRUE, stringsAsFactors = FALSE)
  names(generifs_basic) <- c("tax_id", "gene_id", "pubmed_ids", "timestamp", "annotation")
  # load statement and execute
  statement <- paste(readLines("get_all_RIF.sql"), collapse = "\n")
  #sqldf(drv = "SQLite","select * from generifs_basic limit 5") # 'head'
  edges <- sqldf(drv = "SQLite", statement) # 43491 as of 28/06/13

  ## remove tmp file
  #file.remove(tmpfile)

  return(edges) # first two columns are nodes, other columns considered edge annotations

}

do.load <- function(PATH, percentage = 5, plot_results = TRUE) 
{
  cat("> START: do.load()\n")
  
  edges <- do.load.edges(PATH)

  ### build data into a graph object
  ## make a bipartite graph of pubmed ids and genes
  # graph.data.frame() takes a data frame where the first two columns are nodes IDs (so each row is an edge)
  # and all the other columns are taken as edge attributes, node attributes are loaded separately - see below
  network <- graph.data.frame(
    d = subset(
      x = edges,
      subset = pubmed_ids %in% (sample(unique(edges$pubmed_ids), size = (abs(length(unique(edges$pubmed_ids)) / 100)) * percentage)),
      select = c("pubmed_ids", "gene_id", "annotation")), # load random subset of edges
    directed = FALSE)
  # add node attributes
  # V(network) # all node objects
  V(network)$shape <- c('circle', 'square')[1 + V(network)$name %in% edges$pubmed_ids]
  V(network)$color <- c('black', 'red')[1 + V(network)$name %in% edges$pubmed_ids] # black dots - genes, red boxes - papers
  V(network)$isa <- c('gene', 'paper')[1 + V(network)$name %in% edges$pubmed_ids]
  V(network)$degree <- scale(degree(network), center = FALSE, scale = TRUE )
  # degree(network)
  if(plot_results) {
    do.plot(
      network,
      title = "5%er",
      size = V(network)$degree,
      shape = V(network)$shape,
      colour = V(network)$color)
  }

  cat("> END: do.load()\n")
  str(network)
  return(network)
}

do.decompose <- function(network, plot_results = TRUE) 
{
  cat("> START: do.decompose()\n")
  ## pick out the biggest component
  network <- decompose.graph(network) # coverts to list of networks (each component as separate element)
  largest <- function(x) {
    # return largest sub-component
    sizes <- sapply(
      X = x,
      FUN = function(x) {return(length(V(x)))}
    )
    max.size <- sapply(
      X = sizes,
      FUN = function(i) {return(i == max(sizes))}
    ) # vector of T/F for (potentially more than one) biggest component(s)
    # just take largest sub-list (biggest component)
    if(sum(max.size == 1)) {
      return(x[max.size][[1]])
    }
    else {
      return(x[[max.size]])
    }
  }
  network <- largest(network)
  ## plot it out
  V(network)$degree <- scale(degree(network), center = FALSE, scale = TRUE)
  if(plot_results) {
    do.plot(
      network,
      title="Largest component",
      size = V(network)$degree,
      shape = V(network)$shape,
      colour = V(network)$color)
  }

  cat("> END: do.decompose()\n")
  return(network)
}


do.subset <- function(network, node_ids) 
{
  ## return igraph instance containing only the nodes listed
  # copying all properties
  network <- induced.subgraph(network, vids = node_ids, impl = "copy_and_delete")
  return(network)
}

do.cocitation <- function(network, plot_results = TRUE) 
{
  cat("> START: do.cocitation()\n")

  ### use cocitation to project graph into single type space
  # colour and shape according to node type
  genes <- get.data.frame(network, what = "vertices")
  genes <- genes[genes$isa == "gene", "name"]

  # do cocitation
  network <- graph.adjacency(
                cocitation(
                  graph = network,
                  v = V(network)),
                mode = "undirected",
                weighted = TRUE)
  # add annotations
  V(network)$degree <- scale(degree(network), center = FALSE, scale = TRUE)
  # black dots - genes, red boxes - papers
  V(network)$shape <- c('circle', 'square')[1 + !(V(network)$name %in% genes)]
  V(network)$color <- c('black', 'red')[1 + !(V(network)$name %in% genes)]
  V(network)$isa <- c('gene', 'paper')[1 + !(V(network)$name %in% genes)]

  # dump "paper" nodes to just examine "gene" nodes
  network <- do.subset(network, node_ids = genes)

  # plot it
  if(plot_results) {
    do.plot(
      network,
      title = "Genes only",
      weight = E(network)$weight,
      size = V(network)$degree,
      shape=V(network)$shape,
      colour = V(network)$color)
  }

  cat("> END: do.cocitation()\n")
  return(network)
}

do.phospho <- function(plot_results = TRUE) 
{
  cat("> START: do.phospho()\n")
  ### construct network, not with random 1%, but with phosphatase/kinase subset
  ## pull annotations for human PPases and kinases from GO, see later for use
  pk.ids <- read.delim(header = TRUE, file = "pk_pp_9606.txt")
  # how many of each gene type did we get?
  print(table(pk.ids[ ,c("go_term")]))
  # phosphoprotein phosphatase activity             protein kinase activity
  # 33                                 195
  ## collect a fresh edge list containing all the information
  edges <- do.load.edges(PATH)
  # subset to phospho interaction network
  edges <- subset(edges, gene_id %in% pk.ids$gene_id)
  # create graph
  network <- graph.data.frame(
    d = edges,
    directed = FALSE)
  V(network)$isa <- c("gene", "paper")[1 + V(network)$name %in% edges$pubmed_ids]

  ## decompose and run cocitation
  network <- do.decompose(network, plot_results = FALSE)
  network <- do.cocitation(network, plot_results = FALSE)

  ## plot
  # black dots = kinases, red boxes = PPases
  V(network)$shape <- c('circle', 'square')[1 + V(network)$name %in% pk.ids[pk.ids$go_term == "phosphoprotein phosphatase activity","gene_id"]]
  V(network)$color <- c('black', 'red')[1 + V(network)$name %in% pk.ids[pk.ids$go_term == "phosphoprotein phosphatase activity","gene_id"]]
  V(network)$degree <- scale(degree(network), center = FALSE, scale = TRUE)
  E(network)$weight.scaled <- scale(E(network)$weight, center = FALSE, scale = TRUE)
  V(network)$label <- V(network)$name # default
  V(network)[rank(degree(network), ties.method = "max") < length(V(network)) - 10]$label <- NA # replace everything but the top 10 with NA
  V(network)$shape.2 <- V(network)$shape
  V(network)[!is.na(V(network)$label)]$shape.2 <- 'none' # is it has a label, then don't plot a shape
  V(network)$label.pp <- V(network)$name # default
  V(network)[V(network)$shape == 'circle']$label.pp <- NA # replace kinases with NA
  V(network)$shape.pp <- V(network)$shape
  V(network)[!is.na(V(network)$label.pp)]$shape.pp <- 'none' # is it has a label, then don't plot a shape
  # plot
  if(plot_results) {
    plot(
      network,                                #the graph to be plotted
      layout = layout.kamada.kawai,           # kamada.kawai|fruchterman.reingold|lgl|auto -  the layout method. see the igraph documentation for details
      main = 'largest component, cocitation, kinases - black dots, PPases - red boxes',                    #specifies the title
      vertex.shape = V(network)$shape,        #nodeshape
      vertex.size = 3,                        #nodesize
      # vertex.label.dist=0.5,    	          #puts the name labels slightly off the dots
      vertex.color = V(network)$color,        #node colour
      # vertex.frame.color='blue', 		        #the color of the border of the dots
      vertex.label.color = V(network)$color,	#the color of the name labels
      vertex.label.font = 2,			            #the font of the name labels
      vertex.label = NA,		                  #specifies the lables of the vertices. in this case the 'name' attribute is used
      vertex.label.cex = 1,			              #specifies the size of the font of the labels. can also be made to vary
      edge.width = E(network)$weight.scaled   #specifies the thickness of the edges
    )
  }

  cat("> END: do.phospho()\n")
  return(network)
}

do.mesh <- function(term = "20161122132200_Wnt_Signaling_Pathway.xml", PATH, plot_results = TRUE) 
{
  cat("> START: do.mesh()\n")

  ### as for do.phospho,
  ### pulling out papers related to a given MeSh term

  # get all pmids linked to mesh term
  res <- xmlParse(term)
  res <- xmlToList(res)
  res <- as.integer(as.vector(c(res$IdList, recursive = TRUE)))
  cat(sprintf('%i PMIDS found for term \"%s\"\n',length(res), term))

  ## load edges, subset to retrieved terms and load graph
  edges <- do.load.edges(PATH)
  edges <- subset(edges, pubmed_ids %in% res)
  network<-graph.data.frame(
    d = edges,
    directed = FALSE)
  V(network)$isa <- c("gene", "paper")[1 + V(network)$name %in% edges$pubmed_ids]

  ## decompose and run cocitation
  network <- do.decompose(network, plot_results = FALSE)
  network <- do.cocitation(network, plot_results = FALSE)

  ## add some community detection calls!
  V(network)$pr <- scale(page.rank(network)$vector, center = FALSE, scale = TRUE) # page rank
  V(network)$walktrap <- walktrap.community(network)$membership # community memberships
  V(network)$eigen <- leading.eigenvector.community(network)$membership
  V(network)$spin <- spinglass.community(network)$membership

  # for plotting
  V(network)$shape <- 'circle'
  V(network)$color <- 'black'
  V(network)$degree <- scale(degree(network), center = FALSE, scale = TRUE)
  E(network)$weight.scaled <- scale(E(network)$weight, center = FALSE, scale = TRUE)
  V(network)$label <- V(network)$name # default
  V(network)[rank(degree(network), ties.method="max") < length(V(network)) - 10]$label <- NA # replace everything but the top 10 with NA
  V(network)$shape.2 <- V(network)$shape
  V(network)[!is.na(V(network)$label)]$shape.2 <- 'none' # is it has a label, then don't plot a shape

  # NB layout.lgl seems to work well with coloring by community
  if(plot_results) {
    plot(
      network,                                  #the graph to be plotted
      layout = layout.lgl,                      # kamada.kawai|fruchterman.reingold|lgl|auto -  the layout method. see the igraph documentation for details
      main = 'cocitation, top 10 by degree',    #specifies the title
      vertex.shape = V(network)$shape.2,        #nodeshape
      vertex.size = V(network)$degree*4,        #nodesize
      # vertex.label.dist = 0.5,  		          #puts the name labels slightly off the dots
      vertex.color = V(network)$spin,           #node colour
      # vertex.frame.color = 'blue', 		        #the color of the border of the dots
      vertex.label.color = V(network)$color,		#the color of the name labels
      vertex.label.font = 2,			              #the font of the name labels
      vertex.label = V(network)$label,		      #specifies the lables of the vertices. in this case the 'name' attribute is used
      vertex.label.cex = 1,			                #specifies the size of the font of the labels. can also be made to vary
      edge.width = E(network)$weight.scaled      #specifies the thickness of the edges
    )
  }
  cat("> END: do.mesh()\n")
  return(network)
}

### calls

# phospho network
network <- do.phospho(plot_results = FALSE)
print(get.data.frame(network, what = "vertices")[ ,c( "name", "degree")])


## as above using MEsh terms
# look up some interesting mesh terms
# http://www.ncbi.nlm.nih.gov/mesh?term=autism
lapply(c("20161122132200_Wnt_Signaling_Pathway.xml", "20161122132200_Autistic_Disorder.xml", "20161122132200_Melanoma.xml", "20161122132200_Hedgehog_Proteins.xml"),
       function(x) {
         do.mesh(term = x, PATH, plot_results = FALSE)
       }
)

# final clean up
#rm(list = ls())
#gc()
