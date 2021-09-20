library(tidyverse)
library(ggtree)
library(ggrepel)
source("../../2021-08-26_save_pptx/scripts/save_pptx.r")
# library("colorspace")

df_meta <- read_csv("../data/df_metadata.csv")
names(df_meta)[names(df_meta)=="Strain_Name_sim"] <- "label"
colors <- c("#8f45b5",
"#aada53",
"#9887bc",
"#5fa15d",
"#c54b54",
"#91c2bb",
"#4c393f",
"#bc8e4e")
names(colors) <- unique(df_meta$Host_sim)
write_csv(tibble(host = names(colors), colors=colors), "../results/colors.csv")

# spike tree
tree <- read.tree("../data/seq_spike_aln.fasta.treefile")
tree <- phytools::reroot(tree, 1)
p <- ggtree(tree, size=0.1)

p$data <- left_join(p$data, df_meta)
x_max <- max(p$data$x)
colors_t <- colors[names(colors) %in% p$data$Host_sim]

p$data$label

p + geom_tiplab(aes(color = Host_sim), show.legend=FALSE)+
 	geom_tippoint(aes(color = Host_sim))+
	# geom_nodelab(aes(x=branch, label=label), vjust=-.5, size=3)+
	geom_text_repel(aes(x=x, label=label), size=2, data=. %>% filter(!isTip), min.segment.length = 0.01, nudge_x = -0.1, segment.size=0.3)+
	# scale_color_discrete_qualitative(palette = "Harmonic", na.translate = F)+
	scale_colour_manual(name="Host", values=colors_t, na.translate=F)+
	xlim(x_max*0.5, x_max*1.1)+
	geom_treescale(x=x_max*0.6, y = 20)+
	guides(color = guide_legend(override.aes = list(size = 5)))+
	NULL

ggsave("../results/spike_mltree.pdf", height = 6*sqrt(2), width = 6, scale=1.2)
save_pptx("../results/spike_mltree.pptx", height = 6*sqrt(2), width = 6)

# full genome tree
tree <- read.tree("../data/seq_full_genome_aln.fasta.treefile")
tree <- phytools::reroot(tree, 1)
p <- ggtree(tree, size = 0.1)

p$data <- left_join(p$data, df_meta)
x_max <- max(p$data$x)
colors_t <- colors[names(colors) %in% p$data$Host_sim]
p + geom_tiplab(aes(color = Host_sim), show.legend=FALSE)+
 	geom_tippoint(aes(color = Host_sim))+
	geom_text_repel(aes(x=x, label=label), size=2, data=. %>% filter(!isTip), min.segment.length = 0.01, nudge_x = -0.1, segment.size=0.3)+
	# scale_color_discrete_qualitative(palette = "Harmonic", na.translate = F)+
	scale_colour_manual(name="Host", values=colors_t, na.translate=F, drop=T)+
	xlim(x_max*0.4, x_max*1.1)+
	geom_treescale(x=x_max*0.6, y = 20)+
	guides(color = guide_legend(override.aes = list(size = 5)))+
	NULL

ggsave("../results/full_genome_mltree.pdf", height = 6*sqrt(2), width = 6, scale=1.2)
save_pptx("../results/full_genome_mltree.pptx", height = 6*sqrt(2)*1.2, width = 6*1.2)
