## Import packages
library(phyloseq)
library(vegan)
library(dplyr)
setwd('/home/vdarbot/Bureau/FROGS/analyses_novogen/')
## Setting variables
## The Phyloseq object (format rdata)
# phyloseq <- ""

## The beta diversity distance matrix file
# distance <- ""

## The experiment variable that you want to analyse
# varExp <- ""

## Create input and parameters dataframe
# params <- data.frame( "phyloseq" = phylose, "distance" = distance, "varExp" = varExp)

## Load data
# the phyloseq object, nammed data in FROGSSTAT Phyloseq Import data
load("data/3j_Brune_v2.Rdata.Rdata")

# Convert sample_data to data.frame
metadata <- read.table("metadata/Galaxy107-[metadata_corrig_es.tsv].tabular", row.names=1, header = TRUE)
metadata <- metadata[metadata$Age == "3j" & metadata$Souche == "Brune",]
metadata <- as.data.frame(metadata)
# the distance matrix file
A        <- read.table("data/3j_Brune_v2.tsv", row.names=1)
dist     <- as.dist(A)

## Multivariate ANOVA performed with adonis
adonis_res <- adonis2(formula = dist ~ as.factor(Famille), data = metadata, permutations = 9999)
adonis_res


load("data/3j_Leghorn.Rdata")

# Convert sample_data to data.frame
metadata <- read.table("metadata/Galaxy107-[metadata_corrig_es.tsv].tabular", row.names=1, header = TRUE)
metadata <- metadata[metadata$Age == "3j" & metadata$Souche == "Leghorn",]
metadata <- as.data.frame(metadata)
# the distance matrix file
A        <- read.table("data/3j_Leghorn.tsv", row.names=1)
dist     <- as.dist(A)

## Multivariate ANOVA performed with adonis
adonis_res <- adonis2(formula = dist ~ as.factor(Famille), data = metadata, permutations = 9999)
adonis_res



load("data/16s_Leghorn.Rdata")

# Convert sample_data to data.frame
metadata <- read.table("metadata/Galaxy107-[metadata_corrig_es.tsv].tabular", row.names=1, header = TRUE)
metadata <- metadata[metadata$Age == "16s" & metadata$Souche == "Leghorn",]
metadata <- as.data.frame(metadata)
# the distance matrix file
A        <- read.table("data/16s_Leghorn.tsv", row.names=1)
dist     <- as.dist(A)

## Multivariate ANOVA performed with adonis
adonis_res <- adonis2(formula = dist ~ as.factor(Famille), data = metadata, permutations = 9999)



load("data/16s_Brune.Rdata")

# Convert sample_data to data.frame
metadata <- read.table("metadata/Galaxy107-[metadata_corrig_es.tsv].tabular", row.names=1, header = TRUE)
metadata <- metadata[metadata$Age == "16s" & metadata$Souche == "Brune",]
metadata <- as.data.frame(metadata)
# the distance matrix file
A        <- read.table("data/16s_Brune.tsv", row.names=1)
dist     <- as.dist(A)






########## Pairwise comparisons
library(vegan)
library(dplyr)
library(phyloseq)

con <- url("https://forgemia.inra.fr/lcauquil/16s/-/raw/main/data/16S_phyloseq.rdata")
load(file = con)
data
data_raref_f <- rarefy_even_depth(data, rngseed = 1234)

x <- data_raref_f
factors = "Famille"
sim.method = 'bray'
p.adjust.m = 'BH'

fact <- sample_data(x)[,colnames(sample_data(x)) == factors]
fact <- unlist(fact@.Data)
co <- combn(unique(as.character(fact)), 2)

pairs2 <- c()
total.DF2 <- c()
F.Model2 <- c()
R22 <- c()
p.value2 <- c()
x1 <- c()

for (i in c(1:ncol(co)))
{
  x1 = distance(subset_samples(x, get(factors) %in% c(co[1, i], co[2, i])), method = sim.method)
  ad2 <- adonis2(x1 ~ fact[fact %in% c(co[1, i], co[2, i])])
  pairs2 <- c(pairs2, paste(co[1, i], "vs", co[2, i]))
  total.DF2 <- c(total.DF2, ad2$Df[3])
  F.Model2 <- c(F.Model2, ad2$F[1])
  R22 <- c(R22, ad2$R2[1])
  p.value2 <- c(p.value2, ad2$`Pr(>F)`[1])
}

p.adjusted2 <- p.adjust(p.value2, method = "BH")
sig2 = c(rep("", length(p.adjusted2)))
sig2[p.adjusted2 <= 0.05] <- "."
sig2[p.adjusted2 <= 0.01] <- "*"
sig2[p.adjusted2 <= 0.001] <- "**"
sig2[p.adjusted2 <= 1e-04] <- "***"
pairw.res2 <- data.frame(pairs2, total.DF2, F.Model2, R22, p.value2, p.adjusted2, sig2)
pairw.res2


