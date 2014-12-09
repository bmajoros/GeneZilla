perl -e 'foreach $i (1..50) {print "\n\n\n===============================================\nCHUNK $i\n\n";system("get-transcripts.pl $i.gff ../chunks/$i.fasta $i.transcripts")}' > tmp.1
