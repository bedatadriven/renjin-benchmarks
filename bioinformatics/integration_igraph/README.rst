
Integration: iGraph
===================

Genes are sequences of DNA that encode functional molecules (such as RNAs or Proteins). 
These functional molecules most of time interact and affect eachothers function. 
Therefor, to understand the role of a specific gene in a disease, scientist need to 
know about all the published functions.

Given many published scientific manuscripts about any given gene, the `National 
Center for Biotechnology Information (NCBI)`_ has introduced `GeneRIF`_. GeneRIF 
provides scientist a simple mechanism to add functional annotation to known genes 
based on scientific publications. These annotations are than processed and approved 
by GeneRIF staff. All these functional annotations are freely accessible and can be 
used by anyone to extract functions of a gene or genes with specific function.

This workflow imports human Phophatases and Kinases and the publications. And uses 
'igraph' to store this information and to plot the relation between these publications 
as well as interaction network between molecules.

.. graphviz::
   :caption: Workflow for integration liver cohort data.

   digraph INTEGRATE_liver {
   		compound=true;

		print [label = "print()"];
		
		data_pk [shape = invhouse, label = "human PPases &  \nkinases from GO"];
		data_sql [shape = invhouse, label = "SQL query  \nget_all_RIF.sql"];
		data_rif [shape = invhouse, label = "RIFs \nfrom entrez"];
		read [shape = box, label = "Read\ngunzip()\nread.delim()"];
		edges [shape = box, label = "Extract edges\nsqldf()"];
		network [shape = box; label = "Make graph object\ngraph.data.frame()   "];
		subset [label = "subset"];
		network2 [shape = box; label = "Make graph object\ngraph.data.frame()   "];
		decomp [shape = box; label = "Decompose\ndecompose.graph()   "];
		largest [label = "select largest     \ncomponent"]
		vertices [shape = box; label = "extract vertices \nget.data.frame()   "];
		cocitation [shape = box; label = "peform cocitation\ngraph.adjacency(cocitation())        "];
		subgrph [shape = box; label = "subset gene nodes\ninduced.subgraph()"];
		nodeprop [label = "change node    \nproperties"];
		getdata [label = "get.data.frame()"];

		subgraph cluster_0 {
        	label = "Load data";
			style = filled;
			color = lightgrey;
			data_rif -> read;
			read -> query;
			query -> edges;
			edges -> network;
			data_sql -> edges;
		}
		subgraph cluster_1 {
        	label = "Decompose";
			style = filled;
			color = lightgrey;
			node [style = filled, color = white];
			decomp -> largest;
		}
		subgraph cluster_2 {
        	label = "Do cocitation";
			style = filled;
			color = lightgrey;
			node [style = filled, color = white];
			vertices -> cocitation;
			cocitation -> subgrph;
		}

		network -> subset [ltail = cluster_0];
		data_pk	-> subset -> network2;
		network2 -> decomp [lhead = cluster_1];
		largest -> vertices [ltail = cluster_1, lhead = cluster_2];
		subgrph-> nodeprop [ltail = cluster_2];
		nodeprop -> getdata -> print;
        
        { rank=same;
        network -> subgrph [style=invis];
        }

		# do.load: data_rif->read->query->edges->network;data_sql->edges;
		# do.mesh
		terms [shape = invhouse; label = "Ontology terms"];
		pubmed [shape = box; label = "Query PubMed"];
		listxml [shape = box; label = "xmlToList(xmlParse())"];
		network3 [shape = box; label = "Make graph object\ngraph.data.frame()   "];
		filter [label = "Filter"];
		terms -> pubmed -> listxml -> filter -> network3;
		edges -> filter [ltail = cluster_0];
		network3 -> decomp [lhead = cluster_1];
	}


Packages and Dependencies
----------------------------

There are 7 packages used in this workflow, which depend on 35 additional packages 
from CRAN (dependencies)

**Used packages:**

* *CRAN*: igraph, XML, R.utils, plyr, reshape, utils, sqldf

**Package dependencies:**

* *CRAN*: magrittr, irlba, Matrix, NMF, lattice, foreach, gridBase, pkgmaker, reshape2, stringr, colorspace, doParallel, digest, rngtools, ggplot2, RColorBrewer, cluster, registry, codetools, iterators, xtable, plyr, Rcpp, stringi, scales, gtable, MASS, proto, dichromat, labeling, munsell, RSQLite, gsubfn, chron, DBI

Data
-------

* Basic GeneRIF data:
  ftp://ftp.ncbi.nih.gov/gene/GeneRIF/generifs_basic.gz
* SQL script written by Ieuan Clay:
  get_all_RIF.sql
* Human kinases and phosphatases protein:
  GeneOntology.org list of human (Taxonomy id: 9606) genes using terms GO:0050222
  (protein kinase activity) and GO:0004721 (phosphoprotein phosphatase activity).

.. _National Center for Biotechnology Information (NCBI): http://www.ncbi.nlm.nih.gov
.. _GeneRIF: http://www.ncbi.nlm.nih.gov/gene/about-generif


License
-------

* Copyright (c) 2015 Ieuan Clay based on code from `genbench <https://github.com/biolion/genbench>`_
* Copyright (c) 2015-2016 BeDataDriven B.V.  License: `GPL version 2 or higher`_

.. _GPL version 2 or higher: http://www.gnu.org/licenses/gpl.html


.. raw:: latex

    \clearpage

