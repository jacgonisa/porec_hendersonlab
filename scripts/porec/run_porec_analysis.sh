
enzy=$1
barcode=$2
#Need to import the singularity module 
#Run this in the nextflow environment 
#echo $enzy'/all.reads.fastq.gz'
echo 'results_$enzy'
nextflow run epi2me-labs/wf-pore-c -profile singularity \
        --fastq /rds/project/rds-CXkCVqbkvMY/demux_output/barcode"$barcode".fastq \
        --chunk_size 20000 \
        --cutter "$enzy" \
        --hi_c \
	--mcool \
	--pairs \
        --coverage \
        --ref 'Col-0.ragtag_chrs.mito.chloro.fa' \
        --threads 19 \ #-resume \
	--out_dir 'results_all_formats_april_'$enzy 
        #--minimap2_settings "-x map-ont -w 13"\
              #-profile standard

