library(Biostrings)
library(tidyverse)
library(seqcombo)
source("../../2021-08-26_save_pptx/scripts/save_pptx.r")
df_meta <- read_csv("../data/df_metadata.csv")
# seq spike selected
# seq_name_few <- df_meta %>% filter(Host_sim %in% "Other wild birds") %>% .$Strain_Name_sim
seq_name_few <- c("ISU73347", "ISU690-4", "F230-2006", "CHN-AH-2004", "HNZK-02", "0256-1-2015")
seq_spike <- readDNAStringSet("../data/seq_spike.fasta")
seq_spike_few <- seq_spike[names(seq_spike) %in% seq_name_few]

writeXStringSet(seq_spike_few, "../data/seq_spike_few.fasta")
system("mafft ../data/seq_spike_few.fasta > ../data/seq_spike_few_aln.fasta")

simplot("../data/seq_spike_few_aln.fasta", 'ISU73347')
ggsave("../results/sim_plot.pdf", width=8, height=6)
save_pptx("../results/sim_plot.pptx", , width=8, height=6)

# seq full genome selected
seq_name_few <- c("ISU73347", "ISU690-4", "F230-2006", "CHN-AH-2004", "HNZK-02", "0256-1-2015")
seq_fg <- readDNAStringSet("../data/seq_full_genome.fasta")
seq_fg_few <- seq_fg[names(seq_fg) %in% seq_name_few]

writeXStringSet(seq_fg_few, "../data/seq_fg_few.fasta")
system("mafft --auto --thread -1 ../data/seq_fg_few.fasta > ../data/seq_fg_few_aln.fasta")

simplot("../data/seq_fg_few_aln.fasta", 'ISU73347')
ggsave("../results/sim_plot_fullgenome.pdf", width=8, height=6)
save_pptx("../results/sim_plot_fullgenome.pptx", width=8, height=6)
