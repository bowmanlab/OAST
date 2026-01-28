#!/bin/bash

## conda activate anvio-8

ncbi-genome-download --assembly-level all \
					archaea \
					-A anme2_genomes_outgroup.txt \
                     --metadata metadata.txt \
					 --flat-output \
					 -s genbank
					 
anvi-script-process-genbank-metadata -m metadata.txt \
                                     -o anme2a \
                                     --output-fasta-txt anme2a-fasta.txt \
                                     --exclude-gene-calls-from-fasta-txt
									 
cp anme2a/*fa ./

rm -rf *db
rm -rf *clean.fasta

## be sure to rename genomes so that they don't start with numbers or contain any characters other than _

for f in `ls *fa`;do
fname=`basename $f .fa` 
mv $f $fname.fasta
done

for f in `ls *contigs*fasta` X231106-OB-CTD-1-B3-DNA_filter_8.fasta;do

fname=`basename $f .fasta`
fname=$(printf '%s' "$fname" | sed 's/[^[:alnum:]_]/_/g')
anvi-script-reformat-fasta $f -o ${fname}_clean.fasta --simplify-names --seq-type NT
anvi-gen-contigs-database -T 12 --force-overwrite --contigs-fasta ${fname}_clean.fasta --project-name $fname -o $fname.db
anvi-run-ncbi-cogs -T 12 -c $fname.db
anvi-run-hmms -T 12 -c $fname.db

done

## pause here and create the external-genomes.txt file.  This could be automated but it's a nice opportunity for a sanity check.
## make sure that the names here don't contain any punctuation beside underscore, and don't start with a number.

anvi-gen-genomes-storage -e external-genomes.txt -o anme2a-GENOMES.db

## just in case you're redoing the analysis

rm -rf anme2a_genomes_compare

anvi-pan-genome -T 12 --force-overwrite --genomes-storage anme2a-GENOMES.db --project-name anme2a_genomes_compare

rm -rf ANI

anvi-compute-genome-similarity -e external-genomes.txt -o ANI -p anme2a_genomes_compare/anme2a_genomes_compare-PAN.db -T 12

### steps from https://merenlab.org/data/spiroplasma-pangenome/

## extract single copy marker genes

anvi-get-sequences-for-gene-clusters -p anme2a_genomes_compare/anme2a_genomes_compare-PAN.db \
                                     -g anme2a-GENOMES.db \
                                     --min-num-genomes-gene-cluster-occurs 27 \
                                     --max-num-genes-from-each-genome 1 \
                                     --concatenate-gene-clusters \
									 --force-overwrite \
                                     --output-file anme2a_genomes_compare/anme2a_genomes_compare-SCGs.fa
									 
## this command requires anvio-9, which is incompatible with ncbi-genome-download, not using output of this anyway

anvi-export-gene-clusters \
  -p anme2a_genomes_compare/anme2a_genomes_compare-PAN.db \
  -o anme2a_genomes_compare/gene_clusters.tsv

## rest is anvio-8
									 
trimal -in anme2a_genomes_compare/anme2a_genomes_compare-SCGs.fa \
       -out anme2a_genomes_compare/anme2a_genomes_compare-SCGs-trim.fa \
       -gt 0.50
	   
iqtree -s anme2a_genomes_compare/anme2a_genomes_compare-SCGs-trim.fa \
	-nt 8 \
	-m WAG \
	-bb 1000 \
	-redo
	
## create layers file and add to db
	
echo -e "item_name\tdata_type\tdata_value" \
         > anme2a_genomes_compare/anme2a_genomes_compare_compare-phylogenomic-layer-order.txt
		 
echo -e "SCGs_Bayesian_Tree\tnewick\t`cat anme2a_genomes_compare/anme2a_genomes_compare-SCGs-trim.fa.contree`" \
        >> anme2a_genomes_compare/anme2a_genomes_compare_compare-phylogenomic-layer-order.txt
		
anvi-import-misc-data -p anme2a_genomes_compare/anme2a_genomes_compare-PAN.db \
                      -t layer_orders \
					  anme2a_genomes_compare/anme2a_genomes_compare_compare-phylogenomic-layer-order.txt
					  
anvi-display-pan -p anme2a_genomes_compare/anme2a_genomes_compare-PAN.db -g anme2a-GENOMES.db -I localhost

anvi-interactive -p anme2a_genomes_compare-PAN.db \
                 -t anme2a_genomes_compare/anme2a_genomes_compare-SCGs-trim.fa.contree \
                 --title "Phylogenomics Tutorial Example #4" \
                 --manual
