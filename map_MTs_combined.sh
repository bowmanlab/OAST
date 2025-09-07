#!/bin/bash

## note that bwa will support interleaved file format

CPUS=64
REF=combined_bins_MT.fasta
PATHIN=/data_store/seq_data/2024_jgi_oast/qc_filtered_data/

python3 build_MT_database.py

bwa index ${REF}

cat orca_salcedo_bins_MT.csv temp.csv > combined_bins_MT.jgicounts.csv

for f in `ls $PATHIN*fastq.gz`;do
	fmap=`basename $f .filtered.fastq.gz`
	bwa mem -t ${CPUS} -p ${REF} ${f} | samtools view -F 260,279,269 -q 20 | gzip > ${fmap}_combined_map.sam.gz
	python3 count_mapped_reads.py ${fmap} combined_bins_MT.jgicounts.csv	
done
