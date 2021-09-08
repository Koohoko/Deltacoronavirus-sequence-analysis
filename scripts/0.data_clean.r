library(tidyverse)
library(Biostrings)

# merge data
files_spike <- list.files("../data/0903-δ属毒株序列信息/δ属毒株序列信息-Fasta格式/S 基因全序/", full.names = T)
files_fullg <- list.files("../data/0903-δ属毒株序列信息/δ属毒株序列信息-Fasta格式/全序/", full.names = T)

seq_spike <- lapply(files_spike, function(x){
	readDNAStringSet(x)
})
seq_spike <- do.call("c", seq_spike)
seq_fullg <- lapply(files_fullg, function(x){
	readDNAStringSet(x)
})
seq_fullg <- do.call("c", seq_fullg)

names(seq_spike) <- gsub("\"", "", names(seq_spike))
names(seq_spike) <- gsub("S gene-", "", names(seq_spike))
names(seq_spike) <- gsub(" *\\(.*\\)", "", names(seq_spike))
names(seq_spike) <- gsub("_", "-", names(seq_spike))

names(seq_fullg) <- gsub("\"", "", names(seq_fullg))
names(seq_fullg) <- gsub(" *\\(.*\\)", "", names(seq_fullg))
names(seq_fullg) <- gsub("_", "-", names(seq_fullg))

# metadata 
df_meta <- readxl::read_excel("../data/0903-δ属毒株序列信息/0903-δ属毒株序列信息.xlsx")
df_meta <- df_meta %>% filter(!is.na(`Country（国家）`))
names(df_meta) <- gsub(" *（.*\\）", "", names(df_meta))
names(df_meta) <- gsub(" ", "_", names(df_meta))
df_meta$Strain_Name <- gsub("/", "-", df_meta$Strain_Name)
df_meta$Strain_Name <- gsub("_", "-", df_meta$Strain_Name)
df_meta$Strain_Name <- gsub("_", "", df_meta$Strain_Name)

df_meta$Strain_Name_sim <- df_meta$Strain_Name
df_meta$Strain_Name_sim <- gsub("Haiti-Human-", "", df_meta$Strain_Name_sim)
df_meta$Strain_Name_sim <- gsub("PDCoV-", "", df_meta$Strain_Name_sim)
df_meta$Strain_Name_sim <- gsub("USA-", "", df_meta$Strain_Name_sim)
df_meta$Strain_Name_sim <- gsub("Guangxi-", "", df_meta$Strain_Name_sim)
df_meta$Strain_Name_sim <- gsub("Swine-Vietnam-", "", df_meta$Strain_Name_sim)
df_meta$Strain_Name_sim <- gsub("Swine-Thailand-", "", df_meta$Strain_Name_sim)

df_meta$Country[df_meta$Country == "Viet-Nam"] <- "Vietnam"

df_meta$Host_sim <- df_meta$Host
df_meta$Host_sim <- gsub(" *（.*）", "", df_meta$Host_sim)
df_meta$Host_sim <- gsub("（夜鹭)", "", df_meta$Host_sim)
df_meta$Host_sim <- gsub("\\(鸡\\)", "", df_meta$Host_sim)
df_meta$Host_sim[1:14] <- "Other wild birds"

write_csv(df_meta, "../data/df_metadata.csv")

# reconcile strain name
check <- sapply(names(seq_fullg), function(x){
	tmp <- sapply(df_meta$Strain_Name_sim, function(y){
		grepl(y, x)
	})
	if(sum(tmp)==0){return(NA)}
	which(tmp)
})
name_seq_fullg <- df_meta$Strain_Name_sim[check]
names(seq_fullg)[!is.na(name_seq_fullg)] <- name_seq_fullg[!is.na(name_seq_fullg)]

names(seq_fullg)[!names(seq_fullg) %in% df_meta$Strain_Name_sim]
names(seq_fullg)[names(seq_fullg) == "Chinese ferret badger coronavirus"] <- df_meta$Strain_Name_sim[df_meta$Host=="Chinese Ferret Badger（中国雪貂）"]
names(seq_fullg)[names(seq_fullg) == "PDCoV-2016-Lao"] <- "P1-16-BTL-0115-2016-Lao"
names(seq_fullg)[names(seq_fullg) == "PDCoV-CH01"] <- "CH-01"
names(seq_fullg)[names(seq_fullg) == "PDCoV-Swine-Thailand-S5011"] <- df_meta$Strain_Name_sim[df_meta$Country=="Thailand"]
names(seq_fullg)[names(seq_fullg) == "PDCoV-Swine-Vietnam-Binh21"] <- df_meta$Strain_Name_sim[df_meta$Country=="Vietnam"]
names(seq_fullg)[names(seq_fullg) == "PDCoV-USA-Iowa136"] <- "Iowa136-2015"
names(seq_fullg)[names(seq_fullg) == "PDCoV-USA-Ohio137"] <- "Ohio137-2014"


check <- sapply(names(seq_spike), function(x){
	tmp <- sapply(df_meta$Strain_Name_sim, function(y){
		grepl(y, x)
	})
	if(sum(tmp)==0){return(NA)}
	if(sum(tmp)>1){return(which(df_meta$Strain_Name_sim==x))}
	which(tmp)
})
name_seq_spike <- df_meta$Strain_Name_sim[check]
names(seq_spike)[!is.na(name_seq_spike)] <- name_seq_spike[!is.na(name_seq_spike)]

names(seq_spike)[!names(seq_spike) %in% df_meta$Strain_Name_sim]
names(seq_spike)[names(seq_spike) == "Asian leopard cat coronavirus"] <- df_meta$Strain_Name_sim[df_meta$Host=="Asian leopard cat（亚洲豹猫）"]
names(seq_spike)[names(seq_spike) == "PDCoV-2016-Lao"] <- "P1-16-BTL-0115-2016-Lao"
names(seq_spike)[names(seq_spike) == "PDCoV-CH01"] <- "CH-01"
names(seq_spike)[names(seq_spike) == "PDCoV-Swine-Thailand-S5011"] <- df_meta$Strain_Name_sim[df_meta$Country=="Thailand"]
names(seq_spike)[names(seq_spike) == "PDCoV-Swine-Vietnam-Binh21"] <- df_meta$Strain_Name_sim[df_meta$Country=="Vietnam"]
names(seq_spike)[names(seq_spike) == "PDCoV-USA-Iowa136"] <- "Iowa136-2015"
names(seq_spike)[names(seq_spike) == "PDCoV-USA-Ohio137"] <- "Ohio137-2014"

# add reference genome, IBV
seq_ref <- readDNAStringSet("../data/0903-δ属毒株序列信息/M95169.fasta")
seq_spike_ref <- readDNAStringSet("../data/0903-δ属毒株序列信息/M95169_cds.fasta")
seq_spike_ref <- seq_spike_ref[grepl("spike", tolower(names(seq_spike_ref)))]

names(seq_spike_ref) <- "M95169-Spike"
names(seq_ref) <- "M95169-IBV"

writeXStringSet(c(seq_ref, seq_fullg), "../data/seq_full_genome.fasta")
writeXStringSet(c(seq_spike_ref, seq_spike), "../data/seq_spike.fasta")
