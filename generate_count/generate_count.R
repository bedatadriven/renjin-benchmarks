
##### set up session #####
#rm(list=ls())
set.seed(1000)

# Load used libraries
library(BiocGenerics)
library(GenomicRanges)
library(IRanges)
library(Rsamtools)
library(parallel)
library(rlogging)
SetLogFile("run.log")
options(warn=-1)
cat(">\tLoaded libraries...\n")
# Set options and variables
options( srapply_fapply = "parallel", mc.cores = detectCores() )
annotation_file <- 'mm9_genes.gtf'
bamfiles <- c( 'D2Q2.bam', 'D2Q3.bam','D3Q2.bam', 'D3Q3.bam' )

# Read in annotation
annotation <- read.table( annotation_file, sep = '\t', header = FALSE )

# use only exon coordinates
annotation <- annotation[ which( annotation[,3] == 'exon' ), ]

# IRanges object from exon coordinates
exons <- lapply(
	unique( annotation[,1] ),
	function( chr, ann ) {
		i <- which( ann[,1] == chr );
		IRanges( start = ann[i,4], end = ann[i,5], names = ann[i,9]) ;
	},
	ann = annotation
	)
names(exons) <- unique(annotation[,1])

# removed chrX_random chromosomes
exons$chr1_random <- NULL
exons$chr4_random <- NULL
exons$chr5_random <- NULL
exons$chr7_random <- NULL
exons$chr8_random <- NULL
exons$chr9_random <- NULL
exons$chr13_random <- NULL
exons$chr17_random <- NULL
exons$chrX_random <- NULL
exons$chrY_random <- NULL
exons$chrUn_random <- NULL

# Convert list of chromosome IRanges objects to a RangesList
exonsRangesList<- RangesList()
for( x in names(exons) ) {
	exonsRangesList[[x]] <- exons[[x]]
}

what <- c('pos')

# Read bam files and calculate counts
param <- ScanBamParam( what = "pos", which = exonsRangesList )
#bam <- NULL
#counts <- list()
bam <- scanBam( bamfiles[1], param = param )
counts <- lapply( bam, function(x){ length(x[[1]]) } )
grng <- as(exonsRangesList, 'GRanges')
elementMetadata(grng)[['name']] <- unlist( lapply( exonsRangesList,names))
elementMetadata(grng)[[sub(".bam","", bamfiles[1])]] <- unlist(counts)

for ( i in 1:length(bamfiles) ) {
	if ( i == 1 ) {
		bam <- scanBam( bamfiles[i], param = param )
		counts <- lapply( bam, function(x){ length(x[[1]]) } )
		grng <- as(exonsRangesList, 'GRanges')
		elementMetadata(grng)[['name']] <- unlist( lapply( exonsRangesList,names))
		elementMetadata(grng)[[sub(".bam","", bamfiles[i])]] <- unlist(counts)
	} else {
		bam <- scanBam( bamfiles[i], param=param )
		counts <- lapply( bam, function(x){ length(x[[1]]) } )
		elementMetadata(grng)[[sub(".bam","", bamfiles[i])]] <- unlist(counts)
	}
}


# Aggregate exon counts to gene level
dt <- as.data.frame( grng )
dt$gene_name <- unlist(lapply( strsplit(dt$name, ';'), function(x){ strsplit(x[[1]], " ")[[1]][2] } ))
data <- aggregate( cbind(dt$D2Q2, dt$D2Q3, dt$D3Q2, dt$D3Q3 ), by=list(dt$gene_name), FUN=sum, na.rm=TRUE )
colnames(data) <- c('gene', colnames(dt)[7:ncol(dt)])


print(data[nrow(data),ncol(data)])
