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

# Global CA model --------------------------------------------------------

## CA on CU
tuco.coa <- dudi.coa(df_cu, scannf = FALSE, nf = 5)
(pti <-100*tuco.coa$eig[1:5]/sum(tuco.coa$eig))
cumsum(pti)

get_eigenvalue(tuco.coa)
fviz_screeplot(tuco.coa,addlabels = TRUE, ylim = c(0, 30))

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

# within and between CA ---------------------------------------------------

ttuco.coa <- dudi.coa(t(df_cu), scannf = FALSE, nf = 5)
fac_codon <- factor(Biostrings::GENETIC_CODE[colnames(df_cu)])
ttuco.wca <- wca(ttuco.coa, fac_codon, scan = FALSE, nf = 5)
ttuco.bca <- bca(ttuco.coa, fac_codon, scan = FALSE, nf = 5)


##  variability at the synonymous level & at the amino acid level
100*sum(ttuco.wca$eig)/sum(ttuco.coa$eig)
100*sum(ttuco.bca$eig)/sum(ttuco.coa$eig)

##  express the contributions to the structured variability
(nrow(df_cu) - 1)*(ncol(df_cu) - 1)/sum(df_cu) -> exptoti
(nrow(df_cu) - 1)*(ncol(df_cu) - length(levels(fac_codon)))/sum(df_cu) ->
    exptotiw
(nrow(df_cu) - 1)*(length(levels(fac_codon)) - 1)/sum(df_cu) -> exptotib
### all.equal(exptoti, exptotiw + exptotib) is TRUE
###  there is more variability taken into account at the synonymous level than at the amino acid level
100*(sum(ttuco.wca$eig) - exptotiw)/(sum(ttuco.coa$eig) - exptoti)
100*(sum(ttuco.bca$eig) - exptotib)/(sum(ttuco.coa$eig) - exptoti)


# Synonymous codon usage (WCA) --------------------------------------------
## F1
### Coding sequences point of view
screeplot(ttuco.wca)
(pti <-100*ttuco.wca$eig[1:5]/sum(ttuco.wca$eig))
cumsum(pti)

F1 <- ttuco.wca$co[,1]
hist(F1)
ggplot()+
    geom_density(aes(x = F1, fill = meta_data$Gene), alpha = 0.3) 
ggplot()+
    geom_density(aes(x = F1, fill = meta_data$Host), alpha = 0.3) 
ggplot(meta_data)+
    geom_density(aes(x = F1, fill = meta_data$Host), alpha = 0.8) +
    facet_wrap(vars(Gene))
ggplot(meta_data)+
    geom_density(aes(x = F1, fill = meta_data$Host), alpha = 0.8) +
    facet_grid(vars(Host), vars(Gene), scales = "free_y")

ggplot(meta_data)+ # gc_content
    geom_point(aes(x = F1, y = gc_content, color = meta_data$Gene), 
               alpha = 0.3)
tmp <- lm(F1~gc_content)
summary(tmp)

F2 <- ttuco.wca$co[,2]
F3 <- ttuco.wca$co[,3]

ggplot(meta_data)+
    geom_density(aes(x = F2, fill = meta_data$Gene), alpha = 0.8) +
    facet_wrap(vars(Gene))

ggplot(meta_data)+
    geom_density(aes(x = F2, fill = meta_data$Host), alpha = 0.8) + 
    facet_wrap(vars(Host))

ggplot(meta_data)+
    geom_density(aes(x = F2, fill = meta_data$Host), alpha = 0.8) + 
    facet_wrap(vars(Gene))

tmp_df <- bind_cols(meta_data, F1=F1, F2=F2)
tmp_df$Type <- ifelse(tmp_df$`Virus Species`=="2019-nCoV", "2019-nCoV", "Others")
tmp_df$Type <- factor(tmp_df$Type, levels = c("Others", "2019-nCoV"))
idx_s <- which(tmp_df$id=="MN908947" & tmp_df$Gene=="spike")
tmp_df <- tmp_df %>% mutate(dist_s = sqrt((F1-F1[idx_s])^2+(F2-F2[idx_s])^2))

ggplot(tmp_df)+ 
    geom_point(aes(x = F1, y = F2, color = Host, shape = Type), 
               alpha = 0.3) +
    geom_text_repel(aes(x = F1, y = F2, label = "2019-nCoV"), 
                    data = filter(tmp_df, id=="MN908947"),
                    nudge_x = 0.5, nudge_y = 0.5)+
    geom_text_repel(aes(x = F1, y = F2, label = Host), 
                    data = tmp_df %>% filter(Gene == tmp_df$Gene[idx_s],
                                             `Virus Species`!="2019-nCoV") %>% 
                        arrange(dist_s) %>% .[1:2,],
                    nudge_x = 0.5, nudge_y = -0.1)+
    scale_color_viridis_d()

ggplot(tmp_df)+ 
    geom_point(aes(x = F1, y = F2, color = `Virus Species`, shape = Type), 
               alpha = 0.3) +
    geom_text_repel(aes(x = F1, y = F2, label = "2019-nCoV"), 
                    data = filter(tmp_df, id=="MN908947"),
                    nudge_x = 0.5, nudge_y = 0.5)+
    geom_text_repel(aes(x = F1, y = F2, label = `Virus Species`), 
                    data = tmp_df %>% filter(Gene == tmp_df$Gene[idx_s],
                                             `Virus Species`!="2019-nCoV") %>% 
                        arrange(dist_s) %>% .[1:2,],
                    nudge_x = 0.5, nudge_y = -0.1)+
    scale_color_viridis_d()

ggplot(meta_data)+ # synonymous codon usage difference between avian and swine
    geom_point(aes(x = F1, y = F2, color = meta_data$species_major), 
               alpha = 0.3) +
    facet_wrap(vars(Gene))

### Codon point of view
### F1
x <- ttuco.wca$li[,1]
total_num <- apply(df_cu,2,sum)
x_levels <- rownames(tuco.coa$co)[order(x)]
total_num <- total_num[order(x)]
y <- factor(x_levels, levels = x_levels)
F1 <- x[order(x)]
group <- sapply(x_levels, function(x){
    seqinr::s2c(x)[3]
})
# group <- ifelse(group %in% c("A", "T"), "A or T", "C or G")

ggplot() + # codon end with C 
    geom_point(aes(x = F1, y = y, color = group,
                   size = total_num), alpha = 0.5) +
    ylab("")
ggsave("../results/CLEVELANDs_dot_plot_WCA_F1.tiff",
       width = 5,
       height = 9,
       dpi = 300,
       compress = "lzw")   
### F2
x <- ttuco.wca$li[,2]
total_num <- apply(df_cu,2,sum)
x_levels <- rownames(tuco.coa$co)[order(x)]
total_num <- total_num[order(x)]
y <- factor(x_levels, levels = x_levels)
F2 <- x[order(x)]
group <- sapply(x_levels, function(x){
    seqinr::s2c(x)[3]
})
# group <- ifelse(group %in% c("C", "T"), "C or T", "A or G")

ggplot() + # codon end with A
    geom_point(aes(x = F2, y = y, color = group,
                   size = total_num), alpha = 0.5) +
    ylab("")
ggsave("../results/CLEVELANDs_dot_plot_WCA_F2.tiff",
       width = 5,
       height = 9,
       dpi = 300,
       compress = "lzw")  

# Amino acid usage (BCA) --------------------------------------------------
screeplot(ttuco.bca)
(pti <-100*ttuco.bca$eig[1:5]/sum(ttuco.bca$eig))
cumsum(pti)

### Coding sequences point of view
#### F1 F2
F1 <- ttuco.bca$co[,1]
F2 <- ttuco.bca$co[,2]

ggplot()+ # nucleocapsid
    geom_density(aes(x = F1, fill = meta_data$Gene), alpha = 0.3) 
ggplot()+ # membrane and spike
    geom_density(aes(x = F2, fill = meta_data$Gene), alpha = 0.3) 

ggplot() +
    geom_point(aes(x = F1, y = kd, color = meta_data$Gene), alpha = 0.3) + 
    geom_smooth(aes(x = F1, y = kd), method = "lm")
tmp <- lm(kd~F1)
summary(tmp)

# ggplot(meta_data)+
#     geom_density(aes(x = F1, fill = meta_data$Gene), alpha = 0.8) +
#     facet_wrap(vars(Gene))

tmp_df <- bind_cols(meta_data, F1=F1, F2=F2)
tmp_df$Type <- ifelse(tmp_df$`Virus Species`=="2019-nCoV", "2019-nCoV", "Others")
tmp_df$Type <- factor(tmp_df$Type, levels = c("Others", "2019-nCoV"))
idx_s <- which(tmp_df$id=="MN908947" & tmp_df$Gene=="spike")
tmp_df <- tmp_df %>% mutate(dist_s = sqrt((F1-F1[idx_s])^2+(F2-F2[idx_s])^2))

ggplot(tmp_df)+ 
    geom_point(aes(x = F1, y = F2, color = Host, shape = Type), 
               alpha = 0.3) +
    geom_text_repel(aes(x = F1, y = F2, label = "2019-nCoV"), 
                    data = filter(tmp_df, id=="MN908947"),
                    nudge_x = 0.5, nudge_y = 0.5)+
    geom_text_repel(aes(x = F1, y = F2, label = Host), 
                    data = tmp_df %>% filter(Gene == tmp_df$Gene[idx_s],
                                             `Virus Species`!="2019-nCoV") %>% 
                        arrange(dist_s) %>% .[1:2,],
                    nudge_x = 0.5, nudge_y = -0.1)+
    scale_color_viridis_d()

ggplot(tmp_df)+ 
    geom_point(aes(x = F1, y = F2, color = `Virus Species`, shape = Type), 
               alpha = 0.3) +
    geom_text_repel(aes(x = F1, y = F2, label = "2019-nCoV"), 
                    data = filter(tmp_df, id=="MN908947"),
                    nudge_x = 0.5, nudge_y = 0.5)+
    geom_text_repel(aes(x = F1, y = F2, label = `Virus Species`), 
                    data = tmp_df %>% filter(Gene == tmp_df$Gene[idx_s],
                                             `Virus Species`!="2019-nCoV") %>% 
                        arrange(dist_s) %>% .[1:2,],
                    nudge_x = 0.5, nudge_y = -0.1)+
    scale_color_viridis_d()

### Amino acid point of view
#### F1 F2
F1_aa <- ttuco.bca$li[,1]
F2_aa <- ttuco.bca$li[,2]
levels(fac_codon)
ttuco.bca$lw
EXP$KD ##TODO

ggplot()+ 
    geom_point(aes(x = F1, y = F2, color = meta_data$Gene), 
               alpha = 0.3, data = meta_data)+
    geom_point(aes(x = F1_aa, y = F2_aa, size = 10*ttuco.bca$lw),alpha = 0.3)
