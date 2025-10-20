nextflow run nf-core/hic \
   -profile singularity \
   --input samplesheet_Sakamoto_2.csv \
   --fasta /rds/project/rds-CXkCVqbkvMY/Col-0.ragtag_chrs.mito.chloro.fa \
   --outdir DRR327470_Sakamoto -resume --digestion "dpnii"
