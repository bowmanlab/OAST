#!/bin/bash

## conda activate anvio-8

## be sure to rename genomes so that they don't start with numbers or contain any characters other than _

rm *db
rm *clean.fasta

for f in `ls *Desulfobia.fasta`;do

fname=`basename $f .fasta`
anvi-script-reformat-fasta $f -o $fname.clean.fasta --simplify-names
anvi-gen-contigs-database -T 12 --force-overwrite --contigs-fasta $fname.clean.fasta --project-name $fname -o $fname.db
anvi-run-ncbi-cogs -T 12 -c $fname.db
anvi-run-hmms -T 12 -c $fname.db

done

## create external-genomes

anvi-gen-genomes-storage -e external-genomes.txt -o orca_select-GENOMES.db

## just in case you're redoing the analysis

rm -rf orca_select_genomes_compare

anvi-pan-genome -T 12 --force-overwrite --genomes-storage orca_select-GENOMES.db --project-name orca_select_genomes_compare

anvi-compute-genome-similarity -e external-genomes.txt -o ANI -p orca_select_genomes_compare/orca_select_genomes_compare-PAN.db -T 12

### steps from https://merenlab.org/data/spiroplasma-pangenome/

## extract single copy marker genes

anvi-get-sequences-for-gene-clusters -p orca_select_genomes_compare/orca_select_genomes_compare-PAN.db \
                                     -g orca_select-GENOMES.db \
                                     --min-num-genomes-gene-cluster-occurs 6 \
                                     --max-num-genes-from-each-genome 1 \
                                     --concatenate-gene-clusters \
                                     --output-file orca_select_genomes_compare/orca_select_genomes_compare-SCGs.fa
									 
trimal -in orca_select_genomes_compare/orca_select_genomes_compare-SCGs.fa \
       -out orca_select_genomes_compare/orca_select_genomes_compare-SCGs-trim.fa \
       -gt 0.50
	   
iqtree -s orca_select_genomes_compare/orca_select_genomes_compare-SCGs-trim.fa \
	-nt 8 \
	-m WAG \
	-bb 1000
	
## create layers file and add to db
	
echo -e "item_name\tdata_type\tdata_value" \
         > orca_select_genomes_compare/orca_select_genomes_compare-phylogenomic-layer-order.txt
		 
echo -e "SCGs_Bayesian_Tree\tnewick\t`cat orca_select_genomes_compare/orca_select_genomes_compare-SCGs-trim.fa.contree`" \
        >> orca_select_genomes_compare/orca_select_genomes_compare-phylogenomic-layer-order.txt
		
anvi-import-misc-data -p orca_select_genomes_compare/orca_select_genomes_compare-PAN.db \
                      -t layer_orders orca_select_genomes_compare/orca_select_genomes_compare-phylogenomic-layer-order.txt
					  
anvi-display-pan -p orca_select_genomes_compare/orca_select_genomes_compare-PAN.db -g orca_select-GENOMES.db -I localhost