mafft all_seqs.fasta > all_seqs_aligned.fasta
iqtree -s all_seqs_aligned.fasta -o NC_009527.1_eblv1 -asr -bb 1000
