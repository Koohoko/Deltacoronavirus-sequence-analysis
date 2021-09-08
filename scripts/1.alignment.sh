cd "../data/"

mafft --auto --thread -1 ./seq_spike.fasta > ./seq_spike_aln.fasta
mafft --auto --thread -1 ./seq_full_genome.fasta > ./seq_full_genome_aln.fasta

~/softwares/iqtree-2.1.3-MacOSX/bin/iqtree2 -s ./seq_spike_aln.fasta --alrt 1000 -B 1000 -T AUTO -o M95169-Spike
~/softwares/iqtree-2.1.3-MacOSX/bin/iqtree2 -s ./seq_full_genome_aln.fasta --alrt 1000 -B 1000 -T AUTO -o M95169-IBV