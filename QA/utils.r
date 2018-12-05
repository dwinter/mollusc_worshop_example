#LIbrary-like code for manipulating hapmap results as returned by uneak


#read a (possiblh gzipped) hapmap file into memory 
# returns a data.frame with one row per locus
read_hapmap <- function(fname, gzip=TRUE){
    if(gzip){
        fname = gzfile(fname)
    }
    read.table(fname, header=1, stringsAsFactors=FALSE)
}

# Get indices of inividual (not locus) data for a hapmap data.frame
data_cols <- function(hm){
    12:ncol(hm)
}

# extract individual names from genotype matrix
# NOTE: specific to te formatting used in the this project
get_inds <- function(x){
    sapply(strsplit(colnames(x), "_"), "[[", 1)
}

# extract population names from genotype matrix
# NOTE: specific to te formatting used in the this project
get_pop <- function(x){
    inds <- get_inds(x)
    substr(inds, 1,3)
}

# extract species codes from genotype matrix
# NOTE: specific to te formatting used in the this project
get_spp <- function(x){
    inds <- get_inds(x)
    substr(inds, 3,3)
}

#load IUPAC base names, used in het map and ampa below 
iupac_genotypes <- c("A", "C", "G", "K", "N", "R", "T", "W", "Y", "M", "S")

#map IUPAC base names to het/hom in hapmap
het_map <- structure(c("hom", "hom", "hom", "het", NA, "het", "hom", "het", "het", "het", "het"), 
                     .Names=iupac_genotypes)

# A matrix to map iupac base codes to the number of A,C,G and T bases in a genotype
amat <- cbind(
      #A #C #G #T       
    c( 2, 0, 0, 0), 
    c( 0, 2, 0, 0), 
    c( 0, 0, 2, 0), 
    c( 0, 0, 1, 1), 
    c( 0, 0, 0, 0), 
    c( 1, 0, 1, 0), 
    c( 0, 0, 0, 2), 
    c( 1, 0, 0, 1), 
    c( 0, 1, 0, 1), 
    c( 1, 1, 0, 0), 
    c( 0, 1, 1, 0)
)
colnames(amat) <- iupac_genotypes
rownames(amat) <- c("A", "C", "G", "T")

# map hapmap "Allele" column tp transition or transversion
ti_tv_map <- structure(c("tv",  "ti", "tv",   "tv",  "ti",  "ti",   "tv", "tv",  "tv", "ti"),
                      .Names=c("A/C", "A/G", "A/T", "C/G", "C/T", "G/A", "G/C", "G/T", "T/A", "T/C")
)

# map 'phased' substitution data from allele column to six 'reversable '
sub_spectrum <- structure(c("A/C", "A/G", "A/T", "C/G", "C/T", "A/G", "C/G", "G/T", "A/T", "C/T"),
                .Names=c("A/C", "A/G", "A/T", "C/G", "C/T", "G/A", "G/C", "G/T", "T/A", "T/C")
)




# for a row of genotype data, generate total counts of all 4 possible alleles
allele_counts <- function(locus_data){
    rowSums(amat[,locus_data])
}

# for a row of genotype data, generate allele frequency
# used by AF beloe
allele_freq <- function(locus_data){
    AC <- allele_counts(locus_data)
    AC/sum(AC)
}

# expected heterozygosity based on matrix of observed allele frequencies
He <- function(AF){
    apply(AF, 2, function(x) 1 - sum(x**2))
}

# observed heteozygosity based on called genotypes (takes genotype matrix as input)
Ho <- function(genotypes){
    apply(genotypes, 1, function(x) mean(het_map[x] == "het", na.rm=TRUE))
}

# Allele frequencies from a genotype matrix 
# (returns a 4xN.sample matrix, with rows for frequency A,C,G and T allele)
AF <- function(genotypes){
    apply(genotypes, 1, allele_freq)
}

#Convert 4xN AF matrix into a vector of minor allele frequencies
maf <- function(AF){
    1 - apply(AF, 2, function(x) max(x))
}
