bowtie2 -f --end-to-end --sensitive -a --time --threads 1 -x AB469193.fasta -U AB469193_windows.fasta | samtools view -b - > mapping.bam
~/micropopgen/src/SeqTools/gbk2fasta.pl -i AB469193.gbk -o AB469193.fasta
~/micropopgen/src/grins/detection/produce_windows_from_bam.py --input mapping.bam --output mapping_grins.gff3
