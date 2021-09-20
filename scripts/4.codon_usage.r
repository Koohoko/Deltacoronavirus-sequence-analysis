options(scipen=100)
library(factoextra)
library(tidyverse)
library(ade4)
library(ggplot2)
# library(rayshader)
library(Biostrings)
library(SynMut)
library(seqinr)
# library(future.apply)
library(ggrepel)
library(ggalt)
source("../../2021-08-26_save_pptx/scripts/save_pptx.r")

# input data --------------------------------------------------------------
meta_data_raw <- read_csv("../data/df_metadata.csv")
cds_omsn <- readDNAStringSet("../data/seq_spike_aln.fasta")

colors <- c("#8f45b5",
"#aada53",
"#9887bc",
"#5fa15d",
"#c54b54",
"#91c2bb",
"#4c393f",
"#bc8e4e")
names(colors) <- unique(meta_data_raw$Host_sim)
meta_data_raw$color <- colors[match(meta_data_raw$Host_sim, names(colors))]

df_cu_ori <- as_tibble(get_cu(cds_omsn))

## remove stp codon
df_cu <- df_cu_ori %>% select(-TAA, -TGA, -TAG)

# overall analysis of the genome ------------------------------------------
ggplot(meta_data_raw)+
    geom_histogram(aes(x =`Host_sim`, fill = `Host_sim`), stat="count")+
    # scale_fill_viridis_d(name="Host")+
    scale_fill_manual(name="Host", values=meta_data_raw$color)+
	xlab("Host")+
	ylab("Count")+
    coord_flip()
ggsave("../results/composition_spike_analysis.pdf", width=6, height=6)

## GC content
tmp <- alphabetFrequency(cds_omsn)
gc_content <- apply(tmp, 1, function(x) {
    # (G + C)/(A + T + G + C)
    sum(x[2:3]) / sum(x)
})
ggplot() +
    geom_density(aes(x = gc_content)) +
    geom_vline(xintercept = 0.5, linetype = 2) +
    xlab("GC content")
meta_data <- left_join(meta_data_raw, tibble(Strain_Name_sim=names(cds_omsn), gc_content=gc_content))
writexl::write_xlsx(meta_data, "../results/metadata.xlsx")

df_plot <- meta_data %>% group_by(Host_sim) %>% filter(n()>1)
ggplot(df_plot) +
    geom_density(aes(x = gc_content, fill = Host_sim), n = 60, alpha=0.8) +
    geom_vline(xintercept = 0.5, linetype = 2) +
    facet_wrap(vars(Host_sim), scales = "free_y") +
    # scale_fill_viridis_d() +
	scale_fill_manual(name="Host", values=df_plot$color)+
    xlab("GC content")
ggsave("../results/gc_content_by_host_spike.pdf",
       width = 8, height = 6)
save_pptx("../results/gc_content_by_host_spike.pptx",
       width = 8, height = 6)

# Global CA model --------------------------------------------------------

## CA on CU
tuco.coa <- dudi.coa(df_cu, scannf = FALSE, nf = 5)
(pti <-100*tuco.coa$eig[1:5]/sum(tuco.coa$eig))
cumsum(pti)

get_eigenvalue(tuco.coa)
fviz_screeplot(tuco.coa,addlabels = TRUE, ylim = c(0, 30))
ggsave("../results/codon_usage_spike_scree_plot.pdf", width=8, height=6)

### We may say here,
# for instance, that F1 takes into account 39.2% of the
# structured variability. This approach is not standard but
# could be useful as a safeguard: if the cumulated structured
# variability of the remaining factors is over 100%, then it
# means that we are digging too far and that less factors should
# be considered. 
((nrow(df_cu) - 1)*(ncol(df_cu) - 1)/sum(df_cu) -> exptoti)
(pti2 <-100*tuco.coa$eig[1:5]/(sum(tuco.coa$eig) - exptoti))
cumsum(pti2)

### Coding sequence point of view
row <- get_ca_row(tuco.coa)
head(row$coord)
head(tuco.coa$li[,1])

F1 <- tuco.coa$li[,1]
qual_1 <- row$cos2[,1]
contrib_1 <- row$contrib[,1]

#### Coding sequence point of view
F1 <- tuco.coa$li[,1]
F2 <- tuco.coa$li[,2]
F3 <- tuco.coa$li[,3]
tmp_df <- left_join(meta_data_raw, tibble(Strain_Name_sim=names(cds_omsn), F1=F1, F2=F2))
# tmp_df$Type <- ifelse(tmp_df$`Virus Species`=="2019-nCoV", "2019-nCoV", "Others")
# tmp_df$Type <- factor(tmp_df$Type, levels = c("Others", "2019-nCoV"))
# idx_s <- which(tmp_df$id=="MN908947" & tmp_df$Gene=="spike")
# tmp_df <- tmp_df %>% mutate(dist_s = sqrt((F1-F1[idx_s])^2+(F2-F2[idx_s])^2))

tmp_df <- tmp_df %>% filter(!is.na(F1))
ggplot(tmp_df)+ 
    geom_point(aes(x = F1, y = F2, color = Host_sim), alpha = 0.8) +
    geom_text_repel(aes(x = F1, y = F2, label = Strain_Name_sim), data = filter(tmp_df, Strain_Name_sim %in% c("ISU73347", "F230-2006", "CHN-AH-2004", "HNZK-02", "0256-1-2015")), min.segment.length = 0.1)+
    scale_color_manual(name="Host", values=tmp_df$color)+
	NULL

ggsave("../results/codon_usage_spike.pdf", width=8, height=8)
save_pptx("../results/codon_usage_spike.pptx", width=8, height=8)


# within and between CA ---------------------------------------------------
fac_codon <- factor(Biostrings::GENETIC_CODE[colnames(df_cu)])
tuco <- df_cu

## spike
ttuco.coa_s <- dudi.coa(t(df_cu), scannf = FALSE, nf = 5)
ttuco.wca_s <- wca(ttuco.coa_s, fac_codon, scan = FALSE, nf = 5)
ttuco.bca_s <- bca(ttuco.coa_s, fac_codon, scan = FALSE, nf = 5)
##  variability at the synonymous level & at the amino acid level
100 * sum(ttuco.wca_s$eig) / sum(ttuco.coa_s$eig)
100 * sum(ttuco.bca_s$eig) / sum(ttuco.coa_s$eig)
# Synonymous codon usage (WCA) --------------------------------------------
## F1
### Coding sequences point of view

F1 <- rep(NA, nrow(meta_data))
F1 <- ttuco.wca_s$co[, 1]
F2 <- rep(NA, nrow(meta_data))
F2 <- ttuco.wca_s$co[, 2]
F3 <- rep(NA, nrow(meta_data))
F3 <- ttuco.wca_s$co[, 3]
tmp_df <- left_join(meta_data_raw, tibble(Strain_Name_sim=names(cds_omsn), F1=F1, F2=F2, F3=F3))

kmodel <- kmeans(tmp_df %>% select(F1, F2) %>% filter(!is.na(F1)), centers = 4, nstart = 2000, iter.max = 1000)
tmp_df$cluster <- NA
tmp_df$cluster[!is.na(tmp_df$F1)] <- kmodel$cluster
tmp_df$cluster <- factor(tmp_df$cluster)

ggplot(tmp_df) +
    geom_encircle(aes(x = F1, y = F2, group = cluster), linetype = 2, alpha = 0.6, expand = 0.01, spread = 0.001) +
    geom_point(aes(x = F1, y = F2, color = Host_sim, shape = cluster), alpha = 0.8) +
    geom_text_repel(aes(x = F1, y = F2, label = Strain_Name_sim), data = filter(tmp_df, Strain_Name_sim %in% c("ISU73347", "F230-2006", "CHN-AH-2004", "HNZK-02", "0256-1-2015")), min.segment.length = 0.1)+
    scale_color_manual(name="Host", values=tmp_df$color)+
    scale_shape_discrete(name="Cluster", na.translate = F)+
    scale_x_continuous(limits = c(-max(abs(tmp_df$F1), na.rm=T), max(abs(tmp_df$F1), na.rm=T)))+
    scale_y_continuous(limits = c(-max(abs(tmp_df$F2), na.rm=T), max(abs(tmp_df$F2), na.rm=T)))+
    ggtitle("A. WCA (synonymous codon usage)")

ggsave("../results/codon_usage_spike_WCA.pdf", width=8, height=8)
save_pptx("../results/codon_usage_spike_WCA.pptx", width=8, height=8)

# Amino acid usage (BCA) --------------------------------------------------

### Coding sequences point of view
#### F1 F2
F1 <- rep(NA, nrow(meta_data))
F1 <- ttuco.bca_s$co[, 1]
F2 <- rep(NA, nrow(meta_data))
F2 <- ttuco.bca_s$co[, 2]
F3 <- rep(NA, nrow(meta_data))
F3 <- ttuco.bca_s$co[, 3]
tmp_df <- left_join(meta_data_raw, tibble(Strain_Name_sim=names(cds_omsn), F1=F1, F2=F2, F3=F3))

kmodel <- kmeans(tmp_df %>% select(F1, F2) %>% filter(!is.na(F1)), centers = 4, nstart = 2000, iter.max = 1000)
tmp_df$cluster <- NA
tmp_df$cluster[!is.na(tmp_df$F1)] <- kmodel$cluster
tmp_df$cluster <- factor(tmp_df$cluster)

ggplot(tmp_df) +
    geom_encircle(aes(x = F1, y = F2, group = cluster), linetype = 2, alpha = 0.6, expand = 0.01, spread = 0.001) +
    geom_point(aes(x = F1, y = F2, color = Host_sim, shape = cluster), alpha = 0.8) +
    geom_text_repel(aes(x = F1, y = F2, label = Strain_Name_sim), data = filter(tmp_df, Strain_Name_sim %in% c("ISU73347", "F230-2006", "CHN-AH-2004", "HNZK-02", "0256-1-2015")), min.segment.length = 0.1)+
    scale_color_manual(name="Host", values=tmp_df$color)+
    scale_shape_discrete(name="Cluster", na.translate = F)+
    scale_x_continuous(limits = c(-max(abs(tmp_df$F1), na.rm=T), max(abs(tmp_df$F1), na.rm=T)))+
    scale_y_continuous(limits = c(-max(abs(tmp_df$F2), na.rm=T), max(abs(tmp_df$F2), na.rm=T)))+
    ggtitle("B. BCA (amino acid usage)")

ggsave("../results/codon_usage_spike_BCA.pdf", width=8, height=8)
save_pptx("../results/codon_usage_spike_BCA.pptx", width=8, height=8)
