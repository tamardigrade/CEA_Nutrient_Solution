#SAVING THE ENVIRONMENT IN R
#https://www.techcoil.com/blog/how-to-save-and-load-environment-objects-in-r/
save.image(file='/local1/workdir1/tw488/CEA_RData_RScript/CEA_manuscript_04232026.RData')
#VIEW THE SAVED FILE IN A THE DIRECTORY
#dir()
#TERMINATE R SESSION
#quit(save='no')
#Loading back the entire list environment of objects
load('/local1/workdir1/tw488/CEA_RData_RScript/CEA_manuscript_04232026.RData')

#https://benjjneb.github.io/dada2/tutorial.html
#https://benjjneb.github.io/dada2/ITS_workflow.html
#if (!requireNamespace("BiocManager", quietly = TRUE))y
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Maaslin2")
BiocManager::install("dada2")
BiocManager::install("decontam")
install.packages('theseus')
library("theseus")
library(dada2); packageVersion("dada2")
library(tidyverse)
library(ggplot2)
library(Maaslin2)
library(phyloseq)
library(vegan)
library(permute)
library(reshape2)
library(parallel)
library(pwr)


#workflow from  https://benjjneb.github.io/dada2/ITS_workflow.html
#map the path to the sequencing files on my drive
path <- "/local1/workdir1/tw488/CEA_16S_amplicon/For_Analysis_v2/"
list.files(path)

#generate matched lists of the forward and reverse read files
fnFs <- sort(list.files(path, pattern = "_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fastq.gz", full.names = TRUE))

list(fnFs)
list(fnRs)
#identify the primers used to make these amplicons

FWD <- "GTGCCAGCMGCCGCGGTAA"  ## 
REV <- "GGACTACHVGGGTWTCTAAT"  ##
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

#_________________
#remove poly Ns from reads, create output filtN
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))

#************this step takes a long time
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

#counting primers in F and R reads
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]),multithread = TRUE)

#the above step took awhile too
#what are the results? all F reads in F position, all R reads in R position?

#remove primers with cutadapt- 
cutadapt <-"/workdir1/tw488/CIDA_Sequencing_Runs/cutadapt-venv/bin/cutadapt"
system2(cutadapt, args = "--version")

#create output filenames for the cutadapt-ed files,
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads 
                             "-m", 1, "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

#double check for primers, output should all be 0,

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))


# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2.fastq.gz", full.names = TRUE))
list(cutFs)

q<-basename(file_path_sans_ext(cutFs))
#this only got rid of .gz, need to get rid of .fastq as well
sample.names <-(file_path_sans_ext(q))
#.fastq.gz should be removed
head(sample.names)
list(sample.names)
type(sample.names)



#inspect quality reads
plotQualityProfile(cutFs[1:2])
plotQualityProfile(cutRs[1:2])

#_____start following this workflow: https://benjjneb.github.io/dada2/tutorial.html
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs,  maxN = 0, maxEE = c(2, 2),truncLen=c(230,230), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = FALSE)

list(out)
plotQualityProfile(filtFs[1:2])
plotQualityProfile(filtRs[1:2])
class(out)

list(filtFs)
class(filtFs)
list(filtRs)
class(filtRs)
list(sample.names)

#Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)


#run DADA
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)


#Merge paired reads
#Most of your reads should successfully merge. If that is not the case upstream parameters may need to be revisited: Did you trim away the overlap between your reads?
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
#99% of reads non chimeras!

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
list(track)


#Assign taxonomy
# did you download the DADA2 SILVA dataset or use RDP database in workdir1?
#https://zenodo.org/record/4587955#.ZDVV7XYpDrc
#NOTE: These database files have a known problem in 3/895 families and 59/3936 genera. See https://github.com/mikemc/dada2-reference-databases/blob/main/silva-138.1/v1/bad-taxa.csv for a list of affected taxa and https://github.com/benjjneb/dada2/issues/1293 for more information.
taxa <- assignTaxonomy(seqtab.nochim, "/local1/workdir1/tw488/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


#___________________________________
## 1. ASV Table 
# Prep the asv table! 
samples_out <- rownames(seqtab.nochim)
head(samples_out)
# Pull out sample names from the fastq file name 
sample_names_reformatted <-gsub("_R1.fastq.gz","",samples_out)
head(sample_names_reformatted)


# Replace the names in our seqtable 
#check names and sort
seqtab.nochimdf<-as.data.frame(seqtab.nochim)
sample_names_reformatted.df<-as.data.frame(sample_names_reformatted)
#sort
#replace
rownames(seqtab.nochim) <- sample_names_reformatted
seqtab.nochimdfv2<-as.data.frame(seqtab.nochim)
head(seqtab.nochim)
#do the numbers match up with the new names? If so then continue if not then figure out where the sorting went wrong
### intuition check 
h<-stopifnot(rownames(seqtab.nochim) == sample_names_reformatted)
h
############## Modify the ASV names and then save a fasta file! 
# Give headers more manageable names
# First pull the ASV sequences
asv_seqs <- colnames(seqtab.nochim)
# make headers for our ASV seq fasta file, which will be our asv names
asv_headers <- vector(dim(seqtab.nochim)[2], mode = "character")
# loop through vector and fill it in with ASV names 
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep = "_")
}
# intuition check
asv_headers
##### Rename ASVs in table then write out our ASV fasta file! 
View(seqtab.nochim)
asv_tab <- t(seqtab.nochim)
View(asv_tab)
## Rename our asvs! 
row.names(asv_tab) <- sub(">", "", asv_headers)
View(asv_tab)

## 2. Taxonomy Table 
#prepare-tax-table
View(taxa)
# Add the ASV sequences from the rownames to a column 
new_tax_tab <- taxa %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASVseqs") 
head(new_tax_tab)
# intution check 
h<-stopifnot(new_tax_tab$ASVseqs == colnames(seqtab.nochim))
h
list(asv_tab)
type(asv_tab)


# Now let's add the ASV names 
rownames(new_tax_tab) <- rownames(asv_tab)
View(new_tax_tab)
### Final prep of tax table. Add new column with ASV names 
asv_tax <- 
  new_tax_tab %>%
  # add rownames from count table for phyloseq handoff
  mutate(ASV = rownames(asv_tab)) %>%
  # Resort the columns with select
  dplyr::select(Kingdom, Phylum, Class, Order, Family, Genus, ASV, ASVseqs)
View(asv_tax)
# Intution check
j<-stopifnot(asv_tax$ASV == rownames(asv_tax), rownames(asv_tax) == rownames(asv_tab))
j
#________________________________________________
#import sample df
samdf <-read.csv("/local1/workdir1/tw488/CEA_16S_amplicon/CEA_sampledata_forimport_07092025.csv", header=TRUE,row.names = 1)

#We now construct a phyloseq object directly from the dada2 outputs.
asv_tab<-as.matrix(asv_tab)
asv_tax<-as.matrix(asv_tax)
ps <- phyloseq(otu_table(asv_tab, taxa_are_rows=TRUE), 
               sample_data(samdf), 
               tax_table(as.matrix(asv_tax)))

sample_variables(ps)

ps


#code to remove ASVs==1, NC samples and chloroplast
ps.1 <- ps %>% subset_taxa( Family!= "mitochondria" | is.na(Family) & Class!="Chloroplast" | is.na(Class) )

ps.2 <- subset_samples(ps.1, SampleType!= "C")
ps.2samdata<-ps.2@sam_data

#make ps object of just NC samples
ps.NC <- subset_samples(ps.1, SampleType== "C")

# remove taxa with 1 read count, from https://deneflab.github.io/Diversity_Productivity/analysis/OTU_Removal_Analysis.html
ps.3<-prune_taxa(taxa_sums(ps.2) > 1, ps.2) 

# check to see if above worked by summing rows and columns, as well as 
smdta3<-as.data.frame(ps.3@sam_data)
otutbl3<-as.data.frame(ps.3@otu_table)
txtbl3<-as.data.frame(ps.3@tax_table)
ps.3otu_colsums<-as.data.frame(colSums(otutbl3))
ps.3otu_rowsums<-as.data.frame(rowSums(otutbl3))


#rarefy

ps.noncontam4 = rarefy_even_depth(ps.noncontam3, rngseed=1, sample.size=min(sample_sums(ps.noncontam3)), replace=TRUE, trimOTUs = TRUE, verbose = TRUE)
ps.noncontam4otu<-ps.noncontam4@otu_table
ps.noncontam4samdata<-as.data.frame(ps.noncontam4@sam_data)
ps.noncontam4otu_colsums<-as.data.frame(colSums(ps.noncontam4otu))
ps.noncontam4otu_rowsums<-as.data.frame(rowSums(ps.noncontam4otu))

class(ps.noncontam4otu) <- "matrix" # as.matrix() will do nothing
## you get a warning here, but this is what we need to have
ps.noncontam4otu<- t(ps.noncontam4otu)
ps.noncontam4otu.rarecurve = rarecurve(ps.noncontam4otu, step = 1, label = T)

#check distribution of 
otu_tbl_rowsums<-tibble::rownames_to_column(ps.noncontam4otu_rowsums, "ASV")
colnames(otu_tbl_rowsums)[2] <-"rowsums"
ggplot(otu_tbl_rowsums, aes(x=rowsums)) + geom_histogram()

#Good’s coverage is defined as 1 - (F1/N) where F1 is the number of singleton 
#OTUs and N is the sum of counts for all OTUs.
#Goods coverage for all samples
N=sum(rowSums(ps.noncontam4otu))
N
F1=88
#manually counted b/c its faster
Goodscoverage<-(1-(F1/N))
Goodscoverage
#If a sample has a Good’s coverage == .98, this means that 2% of reads in that 
#sample are from OTUs that appear only once in that sample.

#Goods coverage for Conventional_1_S
ps.noncontam4otu<-as.data.frame(ps.noncontam4@otu_table)
Conventional_1_S<-ps.noncontam4otu %>% select(3)
Conventional_1_S_colsums<-sum(Conventional_1_S)
Conventional_1_S_colsums
singletons=1
singletonsand0s=563
GoodscoverageConventional_1_S<-(1-(singletons/Conventional_1_S_colsums))
GoodscoverageConventional_1_S
GoodscoverageConventional_1_S_0sing<-(1-(singletonsand0s/Conventional_1_S_colsums))
GoodscoverageConventional_1_S_0sing

#Goods coverage for Organic_2_M
ps.noncontam4otu<-as.data.frame(ps.noncontam4@otu_table)
Organic_2_M<-ps.noncontam4otu %>% select(11)
Organic_2_M_colsums<-sum(Organic_2_M)
Organic_2_M_colsums
singletons=0
singletonsand0s=541
GoodscoverageOrganic_2_M<-(1-(singletons/Organic_2_M_colsums))
GoodscoverageOrganic_2_M
GoodscoverageOrganic_2_M_0sing<-(1-(singletonsand0s/Organic_2_M_colsums))
GoodscoverageOrganic_2_M_0sing

#### Calculate relative abundance ####

#count_seqs function for use in calculating relative abundance
count_seqs <- function(pool){
  count_dat <- sample_sums(pool)
  names <- names(count_dat)
  values <- unname(count_dat)
  counts_1 <- cbind(names, values) %>% as.data.frame()
  colnames(counts_1) = c("sample", "counts") 
  
  counts_1 <- counts_1 %>%
    arrange(as.numeric(counts)) %>%
    mutate(counts = as.numeric(counts))
  counts_1$n <- seq.int(nrow(counts_1))
  
  plot <- ggplot(counts_1, aes(n, counts)) +
    geom_line()+
    theme_classic()
  
  return(list("counts"=counts_1, "plot"=plot))
  
}
ps.noncontam4.counts <- count_seqs(ps.noncontam4@otu_table)

abund.ps.noncontam4 <- ps.noncontam4@otu_table %>% as.matrix() %>% as.data.frame() %>%
  rownames_to_column("ASV_") %>%
  tidyr::pivot_longer((-ASV_), names_to = "sample") %>%
  merge(ps.noncontam4.counts$counts, by="sample") %>%
  dplyr::select(-n) %>%
  mutate(newcol=sample)%>%
  tidyr::separate(newcol, c("Pond_Type","Pond_Number","Location"), "_",remove = TRUE)%>%
  group_by(sample) %>%
  mutate(rel_abund = value/counts) %>% #relative abundance by treatment/day/rep
  ungroup() 

#check that relative abundance is calculated correctly- everything should add up to one
abund.ps.noncontam4 %>% group_by(sample) %>% summarise(sum = sum(rel_abund))

#merge taxonomy with rarefied relative abundance
count_tax.ps.noncontam4 <- merge(abund.ps.noncontam4,ps.noncontam4tax, by="ASV_")
head(count_tax.ps.noncontam4)

# Get summary abundance of each taxa
taxa.meta.sumps.noncontam4 <- count_tax.ps.noncontam4%>%
  dplyr::select(-c(value, counts)) %>%
  tidyr::pivot_longer(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "ASV_"),
               names_to = "level", 
               values_to = "taxon")

#format phylum names for genus
genus_rel_abund_ps.noncontam4 <- taxa.meta.sumps.noncontam4 %>%
  filter(level == "Genus") %>%
  group_by(taxon,sample,Pond_Type,Pond_Number,Location) %>% 
  summarise(rel_abund = 100*sum(rel_abund), .groups="drop") %>%
  mutate(taxon = stringr::str_replace(taxon,
                             "^unclassified (.*)", "Unclassified *\\1*"),
         taxon = stringr::str_replace(taxon, "^(\\S*)$", "\\1")) %>%
  ungroup()
genus_rel_abund_ps.noncontam4$taxon <-genus_rel_abund_ps.noncontam4$taxon%>%tidyr::replace_na('Unclassified')


#pool taxon <5% relative abundance into "other" category
taxon_pool_ps.noncontam4 <- genus_rel_abund_ps.noncontam4 %>%
  group_by(sample, taxon, Pond_Type,Pond_Number,Location) %>%
  summarise(mean = mean(rel_abund), .groups="drop") %>%
  group_by(taxon) %>%
  summarise(pool = max(mean) <5,
            mean=mean(mean), 
            .groups="drop") 

#arrange data frame to group samples by relative abundance within each taxon; 
#this will aid visual comparison between the barplots
join_taxon_rel_abund_ps.noncontam4 <- inner_join(genus_rel_abund_ps.noncontam4, taxon_pool_ps.noncontam4, by="taxon") %>%
  mutate(taxon = ifelse(pool, "<5% Relative Abundance", taxon)) %>%
  group_by(sample, taxon,Pond_Type,Pond_Number,Location) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  mutate(taxon = factor(taxon),
         taxon = forcats::fct_reorder(taxon, rel_abund, .desc=TRUE),
         taxon = forcats::fct_shift(taxon, n=1)) %>%
  ungroup() %>%
  unique() %>%
  group_by(taxon) %>%
  arrange(desc(rel_abund), .by_group=TRUE) 



#table of genera for org and convent samples in genus_rel_abund_ps.noncontam4

genera_org<-genus_rel_abund_ps.noncontam4%>%filter(Pond_Type=="Organic")%>%
  group_by(taxon)%>%
  mutate(Avg = mean(rel_abund))%>%distinct(taxon,Avg)%>%
  arrange(desc(Avg))%>%filter(Avg!="0")
colnames(genera_org)[1] <- "Organic Sample Genera"
colnames(genera_org)[2] <- "Average Relative Abundance (%)"
#write.csv(genera_org,"/local1/workdir1/tw488/CEA_RData_RScript/manuscript_data_analysis/Organic_genera_relab182026", row.names = TRUE)


genera_convent<-genus_rel_abund_ps.noncontam4%>%
  filter(Pond_Type=="Conventional")%>%group_by(taxon)%>%
  mutate(Avg = mean(rel_abund))%>%distinct(taxon,Avg)%>%
  arrange(desc(Avg))%>%filter(Avg!="0")
colnames(genera_convent)[1] <- "Conventional Sample Genera"
colnames(genera_convent)[2] <- "Average Relative Abundance (%)"
write.csv(genera_convent,"/local1/workdir1/tw488/CEA_RData_RScript/manuscript_data_analysis/Conventional_genera_relab182026", row.names = TRUE)

#make taxa and sample data tables from phyloseq object
ps.noncontam4_sampledata<-data.frame(ps.noncontam4@sam_data)

#check to make sure all are actually data.frames
class(ps.noncontam4_sampledata)

#DCAST TO CONVERT THE DATA FRAME
#HERE
genus_rel_abund_ps.noncontam4_4maaslin<-dcast(genus_rel_abund_ps.noncontam4 ,sample~taxon,value.var="rel_abund")
genus_rel_abund_ps.noncontam4_4maaslin<-genus_rel_abund_ps.noncontam4_4maaslin%>% column_to_rownames(var="sample")

#run maaslin2 with relative abundance at the genus level
set.seed(12098)
ps.noncontam4_genus_diffab = Maaslin2(input_data = genus_rel_abund_ps.noncontam4_4maaslin, 
                            input_metadata = ps.noncontam4_sampledata,
                            analysis_method = "CPLM",
                            normalization  = "NONE", 
                            transform="NONE",
                            min_prevalence = 0,
                            min_abundance = 0,
                            save_models="TRUE",
                            save_scatter="TRUE",
                            random_effects = c("Pond_Type_Number"),
                            fixed_effects  = "Pond_Type",
                            reference      = c("Pond_Type,Organic"),
                            output         = "/local1/workdir1/tw488/CEA_RData_RScript/Maaslin2/ps.noncontam4_genus_diffab_10122026/",
                            cores=2) 

#Import results in and filter by a <0.05 FDR
ps.noncontam4_genus_diffab_results<-read.table("/local1/workdir1/tw488/CEA_RData_RScript/Maaslin2/ps.noncontam4_genus_diffab_10122026/all_results.tsv",header=TRUE)
ps.noncontam4_genus_diffab_results0.05<-ps.noncontam4_genus_diffab_results%>%filter(qval<0.05)%>%
  select("feature", "coef","qval")%>% mutate(log2_fold_change = log2(exp(coef)))%>%select("feature","qval","log2_fold_change")

genus_rel_abund_ps.noncontam4_TAXAmeanbyOvC<-genus_rel_abund_ps.noncontam4%>% group_by(taxon,Pond_Type) %>% summarise(AvgbyDay= mean(rel_abund))
genus_rel_abund_ps.noncontam4_TAXAmeanbyOvC1<-dcast(genus_rel_abund_ps.noncontam4_TAXAmeanbyOvC,taxon~Pond_Type,value.var="AvgbyDay")

colnames(ps.noncontam4_genus_diffab_results0.05)[1] <- "Genus"
colnames(genus_rel_abund_ps.noncontam4_TAXAmeanbyOvC1)[1] <- "Genus"
ps.noncontam4_genus_diffab_results0.05$Genus <- gsub(".", " ", ps.noncontam4_genus_diffab_results0.05$Genus, fixed = TRUE)
Genus_diffab_relab_table<-merge(genus_rel_abund_ps.noncontam4_TAXAmeanbyOvC1,ps.noncontam4_genus_diffab_results0.05,by="Genus")
colnames(Genus_diffab_relab_table)[2] <- "Conventional Sample Average Relative Abundance (%)"
colnames(Genus_diffab_relab_table)[3] <- "Organic Sample Average Relative Abundance (%)"
colnames(Genus_diffab_relab_table)[4] <- "Significance Level"
colnames(Genus_diffab_relab_table)[5] <- "Log2 Fold Change"
write.csv(Genus_diffab_relab_table,"/local1/workdir1/tw488/CEA_RData_RScript/manuscript_data_analysis/Genus_diffab_relab_table_10122026", row.names = TRUE)


#format phylum names for FAMILY________________________________________________________________

fam_rel_abund_ps.noncontam4 <- taxa.meta.sumps.noncontam4 %>%
  filter(level == "Family") %>%
  group_by(taxon,sample,Pond_Type,Pond_Number,Location) %>% 
  summarise(rel_abund = 100*sum(rel_abund), .groups="drop") %>%
  mutate(taxon = stringr::str_replace(taxon,
                                      "^unclassified (.*)", "Unclassified *\\1*"),
         taxon = stringr::str_replace(taxon, "^(\\S*)$", "\\1")) %>%
  ungroup()
fam_rel_abund_ps.noncontam4$taxon <-fam_rel_abund_ps.noncontam4$taxon%>%tidyr::replace_na('Unclassified')

#check how many taxa there are total
fam_rel_abund_ps.noncontam4_alltaxa <- fam_rel_abund_ps.noncontam4%>% filter(rel_abund!="0") %>%distinct(taxon)
#156 genera plus unclassified

#pool taxon <5% relative abundance into "other" category
taxon_pool_ps.noncontam4 <- fam_rel_abund_ps.noncontam4 %>%
  group_by(sample, taxon, Pond_Type,Pond_Number,Location) %>%
  summarise(mean = mean(rel_abund), .groups="drop") %>%
  group_by(taxon) %>%
  summarise(pool = max(mean) <5,
            mean=mean(mean), 
            .groups="drop") 

#arrange data frame to group samples by relative abundance within each taxon; 
#this will aid visual comparison between the barplots
join_taxon_rel_abund_ps.noncontam4 <- inner_join(fam_rel_abund_ps.noncontam4, taxon_pool_ps.noncontam4, by="taxon") %>%
  mutate(taxon = ifelse(pool, "<5% Relative Abundance", taxon)) %>%
  group_by(sample, taxon,Pond_Type,Pond_Number,Location) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  mutate(taxon = factor(taxon),
         taxon = forcats::fct_reorder(taxon, rel_abund, .desc=TRUE),
         taxon = forcats::fct_shift(taxon, n=1)) %>%
  ungroup() %>%
  unique() %>%
  group_by(taxon) %>%
  arrange(desc(rel_abund), .by_group=TRUE) 


color_map <- c( "<5% Relative Abundance"="black", "Chitinophagaceae" = "#665191",
                "Flavobacteriaceae"="#8DA0CB", "Micrococcaceae"= "#8F9E6F", "Devosiaceae"="#E9C46A", 
                "Pleomorphomonadaceae"= "#4D4D4D", "Lachnospiraceae"= "#E76F51","Mycobacteriaceae"="#8C564B",
                "Spirosomaceae"= "#3A86A8","env.OPS 17"="#AA4499","Sphingomonadaceae"= "#F95D6A",
                "Rhodobacteraceae"= "#009E73","Rhizobiaceae"="#E7298A","Unclassified"="#CDBE70",
                "Xanthobacteraceae" = "gray")

#display Unclassified rel abundances
join_taxon_rel_abund_ps.noncontam4 %>%filter(taxon=="Unclassified")%>%arrange(Pond_Type,rel_abund)
#display >5% rel ab rel abundances
join_taxon_rel_abund_ps.noncontam4 %>%filter(taxon=="<5% Relative Abundance")%>%arrange(Pond_Type,rel_abund)

#### Relative abundance plot ####

join_taxon_rel_abund_ps.noncontam4 %>%
  mutate(taxon = factor(taxon, levels = c(
    "<5% Relative Abundance","Chitinophagaceae","Devosiaceae","env.OPS 17",
    "Flavobacteriaceae","Lachnospiraceae",
    "Micrococcaceae","Mycobacteriaceae","Pleomorphomonadaceae",
    "Rhizobiaceae","Rhodobacteraceae",
    "Sphingomonadaceae","Spirosomaceae",
    "Xanthobacteraceae","Unclassified" ))) %>%
  ggplot() +
  geom_col(mapping = aes(x = sample, y = rel_abund, fill = taxon), 
           color = "black", position = "fill", show.legend = TRUE, width = .4) +
  facet_grid(cols = vars(Pond_Type), scales = "free_x", space = "free_x") +
  ylab("Relative Abundance") +
  xlab("CEA Pond Water Samples") +
  scale_x_discrete(labels = c(
    "Conventional_1_H"="Conventional Pond 1 H", 
    "Conventional_1_M"="Conventional Pond 1 M", 
    "Conventional_1_S"="Conventional Pond 1 S",
    "Conventional_2_H"="Conventional Pond 2 H",
    "Conventional_2_M"="Conventional Pond 2 M", 
    "Conventional_2_S"= "Conventional Pond 2 S",
    "Organic_1_H"="Organic Pond 1 H", 
    "Organic_1_M"="Organic Pond 1 M", 
    "Organic_1_S"="Organic Pond 1 S",
    "Organic_2_H"="Organic Pond 2 H", 
    "Organic_2_M"="Organic Pond 2 M",
    "Organic_2_S"="Organic Pond 2 S" )) +
  theme_linedraw() +
  scale_fill_manual(
    values = color_map,
    labels = c(
      "<5% Relative Abundance" = "<5% Relative Abundance",
      "Unclassified" = "Unclassified",
      "env.OPS 17" = expression(italic("env.OPS 17")),
      "Lachnospiraceae"=expression(italic("Lachnospiraceae")),
      "Chitinophagaceae"=expression(italic("Chitinophagaceae")),
      "Devosiaceae"=expression(italic("Devosiaceae")),
      "Flavobacteriaceae"=expression(italic("Flavobacteriaceae")),
      "Micrococcaceae"=expression(italic("Micrococcaceae")),
      "Mycobacteriaceae"=expression(italic("Mycobacteriaceae")),
      "Pleomorphomonadaceae"=expression(italic("Pleomorphomonadaceae")),
      "Rhizobiaceae"=expression(italic("Rhizobiaceae")),
      "Rhodobacteraceae"=expression(italic("Rhodobacteraceae")),
      "Sphingomonadaceae"=expression(italic("Sphingomonadaceae")),
      "Spirosomaceae"=expression(italic("Spirosomaceae")),
      "Xanthobacteraceae"=expression(italic("Xanthobacteraceae")) )) +
  theme(
    text = element_text(family = "Times New Roman"),
    axis.text.y = element_text(size = 50, color = "black"),
    strip.text = element_text(size = 60),
    axis.title.x = element_text(size = 60, color = "black", margin = margin(t = 30, r = 0, b = 30, l = 0)),
    axis.title.y = element_text(size = 60, color = "black", margin = margin(t = 0, r = 30, b = 0, l = 0)),
    axis.text.x = element_text(size = 50, angle = 90, color = "black"),   
    legend.text = element_text(size = 55), # Removed face = "italic"
    legend.position = "bottom",
    legend.spacing.x = unit(0.1, 'mm'),
    legend.spacing.y = unit(0.05, 'mm'),
    legend.location = "plot",
    legend.title = element_text(face = "bold", size = 55)
  ) +
  guides(fill = guide_legend(ncol = 5, byrow = TRUE)) + 
  labs(fill = 'Bacterial Groups')

ggsave(
  filename = "/local1/workdir1/tw488/CEA_RData_RScript/manuscript_data_analysis/samples_family_relabbarplot10122026.tiff", 
  plot = last_plot(),       
  device = "tiff",          
  units = "in", 
  width = 45, 
  height = 25, 
  dpi = 700, 
  compression = "lzw",
  type = "cairo")



#ALPHA DIVERSITY INDICES AND TFSOS
alphadiv<-estimate_richness(ps.noncontam4, split = TRUE, measures = NULL)
alphadiv$sample <- row.names(alphadiv)  
alphadiv <- alphadiv %>% separate(sample, c("Pond_Type", "Pond_Number", "Sample_Location"), "_")
alphadiv$sample <- row.names(alphadiv) 
alphadiv$Pond_Type_Number = paste(alphadiv$Pond_Type, alphadiv$Pond_Number, sep="_")

#calculate Pielou's Eveness
alphadiv <- alphadiv %>% mutate (Pielou=(alphadiv$Shannon)/log(alphadiv$Observed))

#re normality- there are under 30 data pts. 

#ALPHA DIV
r<-plot_richness(ps.noncontam4, color="Pond_Type", measures=c("Observed","Eveness"))+ 
  geom_boxplot()+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))+
  theme(axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))
r
#ggsave("alphadiv.tiff", units="in", width=12, height=6, dpi=300, compression = 'lzw')




#____________________________________________________________________________________
#Pond Type: Conventional vs Organic
#Sample Location: Harvest, Middle, Seeding
#Need to account for organic/ conventional and which pond it came from, or just organic/ conventional


#PERMANOVA for pond type

ps.noncontotu<-pstoveg_otu(ps.noncontam4)
ps.noncontsample<-pstoveg_sample(ps.noncontam4)
ps.noncontsd<-pstoveg_sd(ps.noncontam4)

#check samples are in the same order in both
all(rownames(ps.noncontotu) == rownames(ps.noncontsd))   #must be TRUE

#make into distance object
ps.noncontotu_vegdist<-vegdist(ps.noncontotu, method="bray")

#combine otu and sd
ps.noncontsd_otu<-cbind(ps.noncontsd,ps.noncontotu)

set.seed(2026)

#Pond type - permute whole ponds
how_pondtype <- with(ps.noncontsd, how(plots  = Plots(strata = Pond_Type_Number, type = "free"),
                                       within = Within(type = "none")))
how_location <- with(ps.noncontsd, how(blocks = Pond_Type_Number, within = Within(type = "free")))

numPerms(nrow(ps.noncontsd), control = how_pondtype)   #should be 24

ps.noncontotu_perm_pondtype <- adonis2(ps.noncontotu_vegdist ~ Pond_Type, data = ps.noncontsd,
                                       permutations = how_pondtype)
ps.noncontotu_perm_pondtype


#ensure replicates are permuted together
testsd<-ps.noncontsd_otu
testsd <- as_tibble(testsd, rownames = "SampleID")
#the number below should be the number of samples you have (the number of rows in testsd)
testsd$permute1<-shuffle(12,control=how_pondtype)
testsd$origorder<-1:12
testsd<-as.data.frame(testsd)

#visualize permutations, change 'group' and 'color' to whatever variable that you want to group samples by
testsd%>%pivot_longer(permute1:origorder, names_to="column", values_to="position", cols_vary = "slowest") %>%
  ggplot(aes(x=column, y=position, group=SampleID, color=Pond_Type_Number))+geom_line()


#check to see even # of samples
table(ps.noncontsd$Pond_Type_Number)
table(ps.noncontsd$Sample_Location)

#Sample location - permute within each pond only 
how_location <- with(ps.noncontsd, how(blocks = Pond_Type_Number, within = Within(type = "free")))
numPerms(nrow(ps.noncontsd), control = how_location)   #should be 1296

ps.noncontotu_perm_location <- adonis2(ps.noncontotu_vegdist ~ Pond_Type_Number + Sample_Location,
                                       data = ps.noncontsd, permutations = how_location,
                                       by = "terms")
ps.noncontotu_perm_location     

#ensure replicates are permuted together
testsd<-ps.noncontsd_otu
testsd <- as_tibble(testsd, rownames = "SampleID")
#the number below should be the number of samples you have (the number of rows in testsd)
testsd$permute1<-shuffle(12,control=how_location)
testsd$origorder<-1:12
testsd<-as.data.frame(testsd)

#visualize permutations, change 'group' and 'color' to whatever variable that you want to group samples by
testsd%>%pivot_longer(permute1:origorder, names_to="column", values_to="position", cols_vary = "slowest") %>%
  ggplot(aes(x=column, y=position, group=SampleID, color=Sample_Location))+geom_line()

#________________________________________________________-
# Simulation-based power analysis for PERMANOVA: conventional vs organic ponds
#
# Purpose: estimate how many ponds per pond type are needed for PERMANOVA to
# detect a difference in community composition with adequate power (e.g., 80%),
# using the present study (2 ponds/type, 3 samples/pond) as pilot data.
#
# Design matched to the study:
#   - Pond is the experimental unit for pond type; samples are nested in ponds.
#   - Each simulated pond gets its own composition (between-pond variation),
#     and each sample within a pond varies around it (within-pond variation).
#   - Each simulated dataset is analyzed exactly as planned: Bray-Curtis +
#     adonis2 with whole ponds permuted (permute::how with Plots).
#
# Model: hierarchical Dirichlet-multinomial
#   pond composition  ~ Dirichlet( pi_type * (1 - theta_b) / theta_b )
#   sample counts     ~ Dirichlet-multinomial( depth, pond comp * (1 - theta_w) / theta_w )

library(vegan)
library(permute)
library(parallel)

set.seed(2026)

##Inputs

# otu  : ps.noncontotu
# data frame with columns SampleID, Pond, PondType, Location: ps.noncontsd

#combined otu and sd
ps.noncontsd_otu


n_grid      <- c(2, 3, 4, 5, 6, 8, 10)  # ponds per pond type to evaluate
w_grid      <- c(0.5, 0.75, 1)          # effect size: 1 = pilot difference, 0.5 = half of it
samples_per_pond <- 3                    # as in the present study
nsim        <- 200                       # simulations per scenario (500+ for final numbers)
nperm       <- 999                       # permutations per PERMANOVA
alpha       <- 0.05
top_taxa    <- 100                       # keep most abundant taxa; rest pooled as "Other"
n_cores     <- max(1, detectCores() - 1) # set to 1 on Windows


## PREPARE PILOT DATA

meta <- ps.noncontsd[order(ps.noncontsd$Pond_Type, ps.noncontsd$Sample_Location), ]
meta <-meta%>%rownames_to_column(var = "SampleID")
otu  <- as.matrix(ps.noncontotu[as.character(meta$SampleID), ])

# Keep the most abundant taxa, pool the rest (reduces sparsity, keeps compositions whole)
ord  <- order(colSums(otu), decreasing = TRUE)
keep <- ord[seq_len(min(top_taxa, ncol(otu)))]
otu  <- cbind(otu[, keep, drop = FALSE],
              Other = rowSums(otu[, -keep, drop = FALSE]))

depth <- round(median(rowSums(otu)))   # sequencing depth used in simulations

meta$Pond     <- meta$Pond_Type_Number          # pond ID: must have 4 levels
meta$PondType <- as.character(meta$Pond_Type)   # Conventional / Organic
meta$Location <- meta$Sample_Location

table(meta$Pond, meta$PondType)   # should show 4 ponds x 3 samples, each under one type

#3. ESTIMATE PARAMETERS

# Weir & Hill method-of-moments estimate of Dirichlet-multinomial overdispersion
weir_theta <- function(X) {
  X  <- X[rowSums(X) > 0, , drop = FALSE]
  J  <- nrow(X); N <- rowSums(X)
  P  <- X / N
  pi <- colSums(X) / sum(X)
  MSP <- sum(N * sweep(P, 2, pi)^2) / (J - 1)
  MSG <- sum(N * P * (1 - P)) / (sum(N) - J)
  Nc  <- (sum(N) - sum(N^2) / sum(N)) / (J - 1)
  th  <- (MSP - MSG) / (MSP + (Nc - 1) * MSG)
  max(th, 1e-6)                         # theta cannot be negative
}

pond_counts <- rowsum(otu, group = meta$Pond)                 # ponds x taxa
pond_type   <- tapply(meta$PondType, meta$Pond, function(x) x[1])
types       <- unique(as.character(pond_type))                # e.g. Conventional, Organic

# Mean composition per pond type (pseudocount avoids zero probabilities)
type_pi <- lapply(types, function(t) {
  x <- colSums(pond_counts[pond_type == t, , drop = FALSE]) + 0.5
  x / sum(x)
})
names(type_pi) <- types

# Between-pond overdispersion within each type (averaged across types)
theta_b <- mean(sapply(types, function(t)
  weir_theta(pond_counts[pond_type == t, , drop = FALSE])))

# Within-pond overdispersion among samples (averaged across ponds)
theta_w <- mean(sapply(split(seq_len(nrow(otu)), meta$Pond), function(i)
  weir_theta(otu[i, , drop = FALSE])))

cat(sprintf("Depth = %d | theta_between = %.4f | theta_within = %.4f\n",
            depth, theta_b, theta_w))
if (length(types) != 2) stop("Expected exactly 2 pond types; found: ", paste(types, collapse = ", "))
if (!all(is.finite(c(theta_b, theta_w))) || theta_b >= 1 || theta_w >= 1)
  stop("Overdispersion estimates are not valid numbers between 0 and 1; check otu and meta.")


# Pilot effect size, analyzed the planned way (for reporting alongside power)
bc_pilot <- vegdist(otu / rowSums(otu), method = "bray")
h_pilot  <- how(plots = Plots(strata = meta$Pond, type = "free"),
                within = Within(type = "none"))
pilot    <- adonis2(bc_pilot ~ PondType, data = meta, permutations = h_pilot)
print(pilot)

## 4. SIMULATION FUNCTIONS
rdirichlet1 <- function(alpha) { g <- rgamma(length(alpha), shape = alpha); g / sum(g) }

sim_dataset <- function(n_ponds, w) {
  pi_A <- type_pi[[1]]
  pi_B <- (1 - w) * type_pi[[1]] + w * type_pi[[2]]   # w scales the pilot difference
  rows <- list(); md <- list(); k <- 0
  for (t in 1:2) {
    pi_t <- if (t == 1) pi_A else pi_B
    for (p in seq_len(n_ponds)) {
      k <- k + 1
      pond_comp <- rdirichlet1(pi_t * (1 - theta_b) / theta_b)
      for (s in seq_len(samples_per_pond)) {
        samp_comp <- rdirichlet1(pond_comp * (1 - theta_w) / theta_w)
        rows[[length(rows) + 1]] <- as.vector(rmultinom(1, depth, samp_comp))
        md[[length(md) + 1]] <- data.frame(Pond = k, Type = types[t])
      }
    }
  }
  list(x = do.call(rbind, rows), md = do.call(rbind, md))
}

one_test <- function(n_ponds, w) {
  d  <- sim_dataset(n_ponds, w)
  bc <- vegdist(d$x / rowSums(d$x), method = "bray")
  h  <- how(plots = Plots(strata = factor(d$md$Pond), type = "free"),
            within = Within(type = "none"), nperm = nperm)
  p  <- suppressMessages(adonis2(bc ~ Type, data = d$md, permutations = h))$`Pr(>F)`[1]
  p < alpha
}

# Single test run outside mclapply: shows any error message in full
print(one_test(4, 1))

#check
rownames(pond_counts)    # 4 pond IDs
pond_type                # 4 entries: 2 Conventional, 2 Organic
c(theta_b, theta_w)      # two numbers between 0 and 1, no NaN
print(one_test(4, 1))    # TRUE or FALSE

## 5. RUN

grid <- expand.grid(n_ponds = n_grid, w = w_grid)
grid$power <- NA
grid$mc_se <- NA

for (i in seq_len(nrow(grid))) {
  hits <- unlist(mclapply(seq_len(nsim), function(j) one_test(grid$n_ponds[i], grid$w[i]),
                          mc.cores = n_cores, mc.set.seed = TRUE))
  grid$power[i] <- mean(hits)
  grid$mc_se[i] <- sqrt(grid$power[i] * (1 - grid$power[i]) / nsim)  # simulation error
  cat(sprintf("n = %2d ponds/type, w = %.2f: power = %.2f\n",
              grid$n_ponds[i], grid$w[i], grid$power[i]))
}

## 6. RESULTS
## ---------------------------------------------------------------------------
print(grid)

# Smallest number of ponds per type reaching 80% power, for each effect size
needed <- sapply(split(grid, grid$w), function(g) {
  ok <- g$n_ponds[g$power >= 0.80]
  if (length(ok)) min(ok) else NA
})
print(needed)   # NA means 80% was not reached within n_grid; extend n_grid

plot(NA, xlim = range(n_grid), ylim = c(0, 1),
     xlab = "Ponds per pond type", ylab = "Power (PERMANOVA, alpha = 0.05)")
abline(h = 0.80, lty = 2)
for (j in seq_along(w_grid)) {
  g <- grid[grid$w == w_grid[j], ]
  lines(g$n_ponds, g$power, type = "b", pch = j, lty = j)
}
legend("bottomright", legend = paste0("Effect = ", w_grid * 100, "% of pilot"),
       pch = seq_along(w_grid), lty = seq_along(w_grid), bty = "n")

write.csv(grid, "permanova_power_results.csv", row.names = FALSE)
sessionInfo()

pwr.anova.test()



## ---- inputs pulled directly from your nested adonis2 object ----
Rsquared_only <- adonis2(ps.noncontotu_vegdist ~ Pond_Type + factor(Pond_Type_Number) + Sample_Location,
                          data = ps.noncontsd, permutations = 0, by = "terms")
tab <- as.data.frame(Rsquared_only)

SS_type <- tab["Pond_Type", "SumOfSqs"];                df_type <- tab["Pond_Type", "Df"]
SS_pond <- tab["factor(Pond_Type_Number)", "SumOfSqs"]; df_pond <- tab["factor(Pond_Type_Number)", "Df"]
stopifnot(length(SS_type) == 1, length(SS_pond) == 1)   # row names must match your model terms

F_type  <- (SS_type / df_type) / (SS_pond / df_pond)                 # ~22.9
omega2p <- (df_type * (F_type - 1)) / (df_type * (F_type - 1) + (df_type + df_pond + 1))
f_obs   <- sqrt(omega2p / (1 - omega2p))                              # Cohen's f

## ---- effect sizes to evaluate (pilot estimate is very uncertain) ----
f_grid <- c(pilot = f_obs, half_pilot = f_obs / 2, large_0.40 = 0.40, medium_0.25 = 0.25)

perm_floor <- 4   # minimum ponds per type for a whole-pond permutation test to reach p < 0.05

res <- do.call(rbind, lapply(names(f_grid), function(nm) {
  f <- f_grid[[nm]]
  n <- tryCatch(pwr.anova.test(k = 2, f = f, sig.level = 0.05, power = 0.80)$n,
                error = function(e) NA)
  data.frame(effect = nm, cohens_f = round(f, 2),
             ponds_per_type_pwr = ceiling(n),
             ponds_per_type_final = max(ceiling(n), perm_floor))
}))
print(c(F = F_type, partial_omega2 = omega2p, f = f_obs))
print(res)

## ---- power curve: power vs ponds per type ----
n_grid <- 2:20
curve <- sapply(f_grid, function(f) sapply(n_grid, function(n)
  pwr.anova.test(k = 2, n = n, f = f, sig.level = 0.05)$power))
curve[n_grid < perm_floor, ] <- 0     # permutation test cannot reach p < 0.05 here
matplot(n_grid, curve, type = "b", pch = 1:4, lty = 1:4, ylim = c(0, 1),
        xlab = "Ponds per pond type", ylab = "Power (alpha = 0.05)")
abline(h = 0.80, lty = 2)
legend("bottomright", legend = paste0(names(f_grid), " (f = ", round(f_grid, 2), ")"),
       pch = 1:4, lty = 1:4, bty = "n")

#k=# of grps: 2, organic and conventional

# n = # observations per group: 3 samples/pond, but 2 ponds per group- should this be 3 or 6 or 2?

pwr.anova.test(k=2,n=3,power=0.80,sig.level=0.05)

#___________________________________________________
#FIND OUT NUM OR PERCENTAGE OF READS THAT HAVE THIS CLASSIFICATION

#make data frame of otus and tax 
ps.noncontam4TAX_withNAs<-as.data.frame(ps.noncontam4@tax_table)

asvscounts<- merge(ps.noncontam4otu,ps.noncontam4TAX_withNAs, by='row.names', all = TRUE)
asvscounts<-column_to_rownames(asvscounts, var = "Row.names")

#filter out asvs with genus only classification
genusonly<-asvscounts%>%filter(Genus!="NA")
#find total # reads w/ family only classification
genusonly_nums <- genusonly[,-(13:20),drop=FALSE]
genusonly_colsums<-as.data.frame(colSums(genusonly_nums))
f<-sum(genusonly_colsums)
f
#178995 ASVs matched to genus classification
#find the number of ASVs matched to genus only
#total ASVs in genusonly
dim(genusonly)
#282 ASvs in genus only but 2 are archaea

#filter out asvs with family only classification
familyonly<-asvscounts%>%filter(Family!='NA')
familyonly<-familyonly[ !(familyonly$Genus %in% genusonly$Genus), ]
#find total # reads w/ family only classification
familyonly_nums <- familyonly[,-(13:20),drop=FALSE]
familyonly_colsums<-as.data.frame(colSums(familyonly_nums))
e<-sum(familyonly_colsums)
e
#7307 ASVs matched to family classification (do not have genus classification)

#find the number of ASVs matched to family only
#total ASVs in familyonly
dim(familyonly)
#176 ASvs in genus only but 2 are archaea


#filter out asvs with order only classification
orderonly<-asvscounts%>%filter(Order!='NA')
orderonly<-orderonly[ !(orderonly$Genus %in% genusonly$Genus), ]
orderonly<-orderonly[ !(orderonly$Family %in% familyonly$Family), ]


#find total # reads w/ order only classification
orderonly_nums <- orderonly[,-(13:20),drop=FALSE]
orderonly_colsums<-as.data.frame(colSums(orderonly_nums))
d<-sum(orderonly_colsums)
d
#2761 ASVs matched to order classification (do not have genus/fam classification)
#find the number of ASVs matched to order only
#total ASVs in orderonly
dim(orderonly)
#88 ASvs in order only no archaea


#find total # reads in this whole data set
asvscounts_nums <- asvscounts[,-(13:20),drop=FALSE]
asvscounts_colsums<-as.data.frame(colSums(asvscounts_nums))
g<-sum(asvscounts_colsums)

#percentage of order only reads
orderpercent=(d/g)*100
#percentage of fam only reads
fampercent=(e/g)*100
#percentage of genus only reads
genuspercent=(f/g)*100

taxlevels<-data.frame("Taxonomic Level"=c('Genus','Family','Order'),
               "Number of ASVs that match each taxonomic level exclusively"=c(282,176,88),
               "Number of ASVs that are archaea"=c(2,2,0),
               "Reads per taxonomic level"=c(f,e,d),
               "Percentage per taxonomic level"=c(genuspercent,fampercent, orderpercent),
               "Total ASvs in rarified data frame"=c(g))
taxlevels
write.csv(taxlevels,"/local1/workdir1/tw488/CEA_RData_RScript/manuscript_data_analysis/taxlevelsdf", row.names = TRUE)

#NC sample analysis, for unrarified NC samples:ps.noncontamNC 
ps.NCsampledata<-as.data.frame(ps.NC@sam_data)
ps.NCtbl<-as.data.frame(ps.NC@otu_table)
ps.NCtxtbl<-as.data.frame(ps.NC@tax_table)
ps.NCotu_colsums<-as.data.frame(colSums(ps.NCtbl))
ps.NCotu_rowsums<-as.data.frame(rowSums(ps.NCtbl))

ps.NC.counts <- count_seqs(ps.NC@otu_table)

abund.ps.noncontamNC <- ps.NC@otu_table %>% as.matrix() %>% as.data.frame() %>%
  rownames_to_column("ASV_") %>%
  tidyr::pivot_longer((-ASV_), names_to = "sample") %>%
  merge(ps.NC.counts$counts, by="sample") %>%
  dplyr::select(-n) %>%
  mutate(newcol=sample)%>%
  tidyr::separate(newcol, c("Pond_Type","Pond_Number","Location"), "_",remove = TRUE)%>%
  group_by(sample) %>%
  mutate(rel_abund = value/counts) %>% #relative abundance by treatment/day/rep
  ungroup() 

#check that relative abundance is calculated correctly- everything should add up to one
abund.ps.noncontamNC %>% group_by(sample) %>% summarise(sum = sum(rel_abund))

#get taxonomy
ps.noncontamNCtax <- ps.NC@tax_table %>%
  as.data.frame() %>%
  rownames_to_column("ASV_") 
head(ps.noncontamNCtax )

#merge taxonomy with rarefied relative abundance
count_tax.ps.noncontamNC <- merge(abund.ps.noncontamNC,ps.noncontamNCtax, by="ASV_")
head(count_tax.ps.noncontamNC)

# Get summary abundance of each taxa
taxa.meta.sumps.noncontamNC <- count_tax.ps.noncontamNC%>%
  dplyr::select(-c(value, counts)) %>%
  tidyr::pivot_longer(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "ASV_"),
                      names_to = "level", 
                      values_to = "taxon")

#format phylum names for GENUS
genus_rel_abund_ps.noncontamNC <- taxa.meta.sumps.noncontamNC %>%
  filter(level == "Genus") %>%
  group_by(taxon,sample,Pond_Type) %>% 
  summarise(rel_abund = 100*sum(rel_abund), .groups="drop") %>%
  mutate(taxon = stringr::str_replace(taxon,
                                      "^unclassified (.*)", "Unclassified *\\1*"),
         taxon = stringr::str_replace(taxon, "^(\\S*)$", "\\1")) %>%
  ungroup()
genus_rel_abund_ps.noncontamNC$taxon <-genus_rel_abund_ps.noncontamNC$taxon%>%tidyr::replace_na('Unclassified')

#check how many taxa there are total
genus_rel_abund_ps.noncontamNC_alltaxa <- genus_rel_abund_ps.noncontamNC%>% filter(rel_abund!="0") %>%distinct(taxon)
#10 "genera"

#filter out genera that have 0% rel ab
#genus_rel_abund_ps.noncontamNC_no0relab <- genus_rel_abund_ps.noncontamNC%>% filter(rel_abund!="0") 



#rel ab of genera for org and convent samples in join_taxon_rel_abund_ps.noncontam4

genera_org<-genus_rel_abund_ps.noncontamNC_no0relab

genus_rel_abund_ps.noncontamNC_no0relab_4export<-dcast(genus_rel_abund_ps.noncontamNC_no0relab,sample~taxon,value.var="rel_abund")
genus_rel_abund_ps.noncontamNC_no0relab_4export<-genus_rel_abund_ps.noncontamNC_no0relab_4export%>%replace(is.na(.), 0)
write.csv(genus_rel_abund_ps.noncontamNC_no0relab_4export,"/local1/workdir1/tw488/CEA_RData_RScript/manuscript_data_analysis/NC_samples_relab", row.names = TRUE)

#family analysis for NC samples
#NC sample analysis, for unrarified NC samples:ps.noncontamNC 

#format phylum names for family
fam_rel_abund_ps.noncontamNC <- taxa.meta.sumps.noncontamNC %>%
  filter(level == "Family") %>%
  group_by(taxon,sample,Pond_Type) %>% 
  summarise(rel_abund = 100*sum(rel_abund), .groups="drop") %>%
  mutate(taxon = stringr::str_replace(taxon,
                                      "^unclassified (.*)", "Unclassified *\\1*"),
         taxon = stringr::str_replace(taxon, "^(\\S*)$", "\\1")) %>%
  ungroup()
fam_rel_abund_ps.noncontamNC$taxon <-fam_rel_abund_ps.noncontamNC$taxon%>%tidyr::replace_na('Unclassified')


#pool taxon <5% relative abundance into "other" category
fam_rel_abund_ps.noncontamNC_pool <- fam_rel_abund_ps.noncontamNC %>%
  group_by(sample, taxon, Pond_Type) %>%
  summarise(mean = mean(rel_abund), .groups="drop") %>%
  group_by(taxon) %>%
  summarise(pool = max(mean) <5,
            mean=mean(mean), 
            .groups="drop")

join_fam_rel_abund_ps.noncontamNC <- inner_join(fam_rel_abund_ps.noncontamNC,fam_rel_abund_ps.noncontamNC_pool, by="taxon") %>%
  mutate(taxon = ifelse(pool, "<5% Relative Abundance", taxon)) %>%
  group_by(sample, taxon,Pond_Type) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  mutate(taxon = factor(taxon),
         taxon = forcats::fct_reorder(taxon, rel_abund, .desc=TRUE),
         taxon = forcats::fct_shift(taxon, n=1)) %>%
  ungroup() %>%
  unique() %>%
  group_by(taxon) %>%
  arrange(desc(rel_abund), .by_group=TRUE) 


#rel ab genera for org and convent samples in join_taxon_rel_abund_ps.noncontam4
#filter out families that have 0% rel ab
fam_rel_abund_ps.noncontamNC_no0relab <- fam_rel_abund_ps.noncontamNC%>% filter(rel_abund!="0") 

fam_rel_abund_ps.noncontamNC_no0relab_4export<-dcast(fam_rel_abund_ps.noncontamNC_no0relab,sample~taxon,value.var="rel_abund")
fam_rel_abund_ps.noncontamNC_no0relab_4export<-fam_rel_abund_ps.noncontamNC_no0relab_4export%>%replace(is.na(.), 0)
write.csv(fam_rel_abund_ps.noncontamNC_no0relab_4export,"/local1/workdir1/tw488/CEA_RData_RScript/manuscript_data_analysis/NC_samples_relab_fam", row.names = TRUE)
