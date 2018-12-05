## Code to reproduce data/figure presented an Genomics workshop

#load packages
library(dendextend, quietly=TRUE)
library(ggplot2, quietly=TRUE)
library(RColorBrewer)


##################################################
# read in SNV file and compare Ti/Tv             #
##################################################

var_spectrum <- table(sub_spectrum[ SNVs$alleles ])
barplot(var_spectrum)


#################################################
## identify missing data                        #
## calculate call-rate for individuals and loci #
#################################################

genotypes <- SNVs[ , data_cols(SNVs)]
missing <- genotypes == "N"
ind_call_rate <- median(1 - colMeans(missing))
loc_call_rate <- median(1 - rowMeans(missing))

###################################################
# Is missingness random with respect to           #
# population or species (i.e. check for allelic   #
# dropout)                                        #
###################################################

#easier to shade tips by population if we first extract
# the population name from the more verbose column names
colnames(missing) <- get_pop(missing)

#Now generate distance matrix and cluster with `hclust`
D <- dist(t(missing))
tr <- as.dendrogram(hclust(D))

#using dndexten to shade tip-lables
sample_pal <- structure(brewer.pal(n=8, name="Set1"), .Names = unique(get_pop(missing)))
labels_colors(tr) <- sample_pal[ labels(D)[order.dendrogram(tr)] ]
plot(tr)

###################################################
# Lots of allelic dropout in species, so focus on #
# one species for rest of this demo               #
###################################################
conv <- genotypes[,get_spp(genotypes) == "c"]
conv <- conv[!is.nan(maf(AF(conv))),]
dim(conv)


###################################################
# Observed v Expected heterozygosity              #
###################################################

conv_afreq  <- AF(conv)
conv_He  <- He(conv_afreq)
conv_Ho  <- Ho(conv)

hets <- data.frame(observed=conv_Ho, expected = conv_He)
ggplot(hets ,aes( observed, expected)) + 
    geom_jitter(width=0.003, height=0.003, alpha=0.2) + 
    geom_abline(colour="red", size=1) + 
    coord_equal()


###################################################
# Site frequency spectrum                         #
###################################################

hist(maf(conv_afreq), breaks=50)


###################################################
# Filter                                          #
###################################################

to_keep <-  (rowMeans(conv == "N") < 0.4) & (conv_Ho <= conv_He) & (maf(conv_afreq) > 0)
sum(to_keep)

###################################################
# Filtered Site frequency spectrum                #
###################################################

hist(maf(AF(conv[to_keep,])), breaks=50)

###################################################
# Horrible hacky PCA                              #
###################################################

filt <- conv[to_keep,]
res <- matrix(0, nrow=1604, ncol=ncol(conv))
for(i in 1:nrow(filt)){
   res[i,] <- colSums( amat[,unlist(filt[i,]) ][-major_allele[i],])
}
PCA <- prcomp(res)
PCA_df <- data.frame(PCA$rotation, pop = get_pop(filt), ind= get_inds(filt))
ggplot(PCA_df, aes(PC3, PC2, colour=pop,label=ind)) + geom_text() + coord_equal()




