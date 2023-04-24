#Clear out R before starting anything
rm(list=ls())
#Set the working directory to my file
setwd("/Volumes/Projects/8/8289c66d-2d56-4706-a307-5a9a3eb3747e/data/Capture_experiment_data/Capture_probes_list")
getwd()
#Import libraries 
library(dplyr)
library(readr)
library(data.table)
library(tximport)
library(tidyr)
library(ggplot2)
library(viridis)
library(ggpmisc)
library(ggExtra)

lncRNA <- read.delim("lincRNA_probes.txt", header=FALSE)
probe_num <- read.delim("probe_gene_overlap.txt", header=FALSE)
colnames(lncRNA) <- "Gene"
colnames(probe_num) <- c("Gene", "Probes")
tx2gene <- import_gtf("./tama2geneid_transcriptome.sequin.gtf")
tx2gene <- move_columns_to_front(df = tx2gene, columns = c("transcript_id", "gene_id"))
probe_num <- dplyr::mutate(probe_num, Category = ifelse((probe_num$Gene%in%lncRNA$Gene), "lncRNA", "Coding"))
txi <- tximport(files, type="salmon", tx2gene=tx2gene[,c("transcript_id", "gene_id")])
# gene <- summarizeToGene(txi,tx2gene=tx2gene[,c("transcript_id", "gene_id")], countsFromAbundance = "scaledTPM")
gene_counts<- as.dataframe(txi$counts)
gene_counts<-  tibble::rownames_to_column(gene_counts, "Gene")
probe_num <- left_join(probe_num,gene_counts,by = "Gene")
probe_num[is.na(probe_num)]<-0
probe_num_long <- gather(probe_num, sample, expression, X1:X49, factor_key=TRUE)

# Scatterplot
theme_set(theme_bw())  # pre-set the bw theme.
probe_num_long %>%
  ggplot(aes(x=Probes, y=expression, group=Category, color=Category)) +
  scale_color_viridis(discrete = TRUE) + geom_point(size = 0.2, alpha = 0.5 ) + 
  scale_fill_viridis(discrete = TRUE ) +
  ylim(-5, 15) +
  stat_poly_line() +
  stat_poly_eq() +
  xlab("Number of probes") + theme_bw() +
  ylab(expression(Measured~TPM~(log[2])))  
probe_num %>% 
  ggplot(aes(x=Probes, color=Category)) +
  scale_color_viridis(discrete = TRUE) +
  geom_histogram(fill="white", alpha=0.5, position="identity") + theme_bw() 
probe_num_log <-  probe_num_long%>% mutate_at(vars(matches("expr")), log2)
is.na(probe_num_log)<-sapply(probe_num_log, is.infinite)
probe_num_log[is.na(probe_num_log)]<-0

group.colors <- c(Coding = "#460B6AFF", lncRNA = "#FCA007FF")
g <- ggplot(probe_num_log, aes(x=Probes, y=expression,group=Category, color=Category)) + 
  geom_count() + 
  # geom_smooth(method="lm", se=F) +
  scale_fill_manual(values=group.colors) + 
  scale_color_manual(values=group.colors) + geom_point(alpha = 1/20) +
  xlab("Number of probes") + theme_bw() +
  stat_poly_eq(label.x = "right",
               label.y = "top") +
  ylab(expression(Measured~TPM~(log[2])))  

ggMarginal(g, type = "histogram", groupColour = TRUE, fill="white")

dev.print(pdf,"probe_per_gene.pdf") 

capture_counts<- read.delim("../capture_counts_WH.txt", header=TRUE)
capture_counts <- left_join(probe_num,capture_counts,by = "Gene")
capture_counts[is.na(capture_counts)]<-0
capture_counts <- capture_counts %>% select(!starts_with("X"))
capture_counts <- gather(capture_counts, sample, expression, 4:52, factor_key=TRUE)
capture_counts <- capture_counts%>% mutate_at(vars(matches("expr")), log2)
is.na(capture_counts)<-sapply(capture_counts, is.infinite)
capture_counts[is.na(capture_counts)]<-0
g <- ggplot(capture_counts, aes(x=Probes, y=expression,group=Category, color=Category)) + 
  geom_count() + 
  scale_color_viridis(discrete = TRUE) + 
  scale_fill_viridis(discrete = TRUE ) + geom_point(alpha = 1/10) +
  xlab("Number of probes") + theme_bw() +
  ylab("expression")  

ggMarginal(g, type = "histogram", groupColour = TRUE, fill="white")

g <- ggplot(capture_counts, aes(x=Probes, y=expression,group=Category, color=Category)) + 
  geom_count() + 
  scale_fill_manual(values=group.colors) + 
  scale_color_manual(values=group.colors) + geom_point(alpha = 1/10) +
  xlab("Number of probes") + theme_bw() +
  ylab("expression")  

ggMarginal(g, type = "histogram", groupColour = TRUE, fill="white")




