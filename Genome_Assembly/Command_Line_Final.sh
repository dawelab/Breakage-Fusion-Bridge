###These commands were run on in an interactive job directly on the command line. 

#Combining the new c2 assembly with the prior ABS assembly to obtain the diploid assembly needed for the alignments
cat AbsGenomePBHIFI_version_1.fa ragtag.scaffold.fasta > ABSxc2v2.fasta

#Remove unplaced contigs from the c2 assembly
samtools faidx ABSxc2v2.fasta
cut -f 1 ABSxc2v2.fasta.fai | grep -v "h2tg" > headers_to_keep.txt
samtools faidx ABSxc2v2.fasta.fasta -r headers_to_keep.txt > ABSxc2v2.1.fasta
