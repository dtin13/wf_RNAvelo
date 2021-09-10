## PSEUDOALIGNMENT

pip install kb-python
kb count -i mm_cDNA_introns_97.idx -g tr2g.tsv -x 10xv3 -o kb -c1 cDNA_tx_to_capture.txt
-c2 introns_tx_to_capture.txt --lamanno sample.R1.fastq.gz sample.R2.fastq.gz
