## PREPROCESSING SCRNA-SEQ DATA USING STARSOLO

STAR-master/source/STAR --genomeDir
"../res/index/" --readFilesIn
"../data/SRR10511597_3.fastq.gz"
"../data/SRR10511597_2.fastq.gz"
--soloType CB_UMI_Simple --soloCBwhitelist
"../res/3M-february-2018.txt"
--soloFeatures Gene Velocyto --readFilesCommand zcat --soloUMIlen 12
