library(tidyverse)
library(ggtree)
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

# spike tree
tree <- read.tree("../data/seq_spike_aln.fasta.treefile")
tree <- phytools::reroot(tree, 1)
p <- ggtree(tree)

p$data <- left_join(p$data, df_meta)
x_max <- max(p$data$x)
colors_t <- colors[names(colors) %in% p$data$Host_sim]
p + geom_tiplab(aes(color = Host_sim), show.legend=FALSE)+
 	geom_tippoint(aes(color = Host_sim))+
	# scale_color_discrete_qualitative(palette = "Harmonic", na.translate = F)+
	scale_colour_manual(name="Host", values=colors_t, na.translate=F)+
	xlim(x_max*0.5, x_max*1.1)+
	geom_treescale(x=x_max*0.6, y = 20)+
	guides(color = guide_legend(override.aes = list(size = 5)))+
	NULL

ggsave("../results/spike_mltree.pdf", height = 8*sqrt(2), width = 8, scale=1.2)

# full genome tree
tree <- read.tree("../data/seq_full_genome_aln.fasta.treefile")
tree <- phytools::reroot(tree, 1)
p <- ggtree(tree)

p$data <- left_join(p$data, df_meta)
x_max <- max(p$data$x)
colors_t <- colors[names(colors) %in% p$data$Host_sim]
p + geom_tiplab(aes(color = Host_sim), show.legend=FALSE)+
 	geom_tippoint(aes(color = Host_sim))+
	# scale_color_discrete_qualitative(palette = "Harmonic", na.translate = F)+
	scale_colour_manual(name="Host", values=colors_t, na.translate=F, drop=T)+
	xlim(x_max*0.4, x_max*1.1)+
	geom_treescale(x=x_max*0.6, y = 20)+
	guides(color = guide_legend(override.aes = list(size = 5)))+
	NULL

ggsave("../results/full_genome_mltree.pdf", height = 8*sqrt(2), width = 8, scale=1.2)

