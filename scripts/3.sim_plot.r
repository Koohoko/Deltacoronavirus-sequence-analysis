library(Biostrings)
library(tidyverse)
library(seqcombo)

df_meta <- read_csv("../data/df_metadata.csv")
# seq spike selected
# seq_name_few <- df_meta %>% filter(Host_sim %in% "Other wild birds") %>% .$Strain_Name_sim
seq_name_few <- c("ISU73347", "F230-2006", "CHN-AH-2004", "HNZK-02", "0256-1-2015")
seq_spike <- readDNAStringSet("../data/seq_spike.fasta")
seq_spike_few <- seq_spike[names(seq_spike) %in% seq_name_few]

writeXStringSet(seq_spike_few, "../data/seq_spike_few.fasta")
system("mafft ../data/seq_spike_few.fasta > ../data/seq_spike_few_aln.fasta")

simplot("../data/seq_spike_few_aln.fasta", 'ISU73347')
ggsave("../results/sim_plot.pdf", width=8, height=6)
