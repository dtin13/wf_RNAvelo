## OBTAINING OF FASTQ FILE

# 1 . Download SRA file
prefetch SRR10511597 --max-size u

# 2. Convert to FASTQ
fastq-dump --origfmt --split-files --gzip SRR10511597
-O /local/ajgsanchez/skin/data_validation/SRR10511597

