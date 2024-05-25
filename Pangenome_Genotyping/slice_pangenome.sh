#!/bin/bash

# Define the input VCF file and the maximum number of variants per file
pangenome_vcf=$1
number_of_slices=$2
output_folder=$3

# Define chromosome names
CHROMOSOMES=('chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chr20' 'chr21' 'chr22' 'chrX' 'chrY')

# Ensure output folder exists
mkdir -p $output_folder

bcftools view  --header-only $pangenome_vcf  > $output_folder/header


 
# Loop over each chromosome
for chrom in ${CHROMOSOMES[@]}
do
    # Extract variants for the current chromosome
    chrom_vcf="$output_folder/${chrom}"
    bcftools view -r $chrom $pangenome_vcf |grep -vP "^#" > $chrom_vcf
    split -d  -n l/$number_of_slices  $chrom_vcf $chrom_vcf.splitted.
done



ls $output_folder/chr1.splitted.* | grep -oP "splitted.*" |cut -f2 -d'.' >  $output_folder/slices

parallel --gnu cp $output_folder/header $output_folder/slice_{}.vcf :::: $output_folder/slices

parallel --gnu cat "$output_folder/chr1.splitted.{} $output_folder/chr2.splitted.{} $output_folder/chr3.splitted.{} $output_folder/chr4.splitted.{} $output_folder/chr5.splitted.{} $output_folder/chr6.splitted.{} $output_folder/chr7.splitted.{} $output_folder/chr8.splitted.{} $output_folder/chr9.splitted.{}   $output_folder/chr10.splitted.{} $output_folder/chr11.splitted.{} $output_folder/chr12.splitted.{}                      $output_folder/chr13.splitted.{} $output_folder/chr14.splitted.{} $output_folder/chr15.splitted.{} $output_folder/chr16.splitted.{} $output_folder/chr17.splitted.{} $output_folder/chr18.splitted.{} $output_folder/chr19.splitted.{} $output_folder/chr20.splitted.{} $output_folder/chr21.splitted.{} $output_folder/chr22.splitted.{} $output_folder/chrX.splitted.{} $output_folder/chrY.splitted.{} >> $output_folder/slice_{}.vcf" :::: $output_folder/slices
