# Fine Mapping of GWAS SNPS using SVs from the 4k reference panel

## 1. Data Preparation
I am getting the gwas snps in vcf format Download both the GWAS catalog and dpsnp. Then 

```
wget https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL.tar
tar -xvf GTEx_Analysis_v8_eQTL.tar


wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz.tbi

cat GWAS_catalog | \
cut -f12,13 | \
sort -k1,1 -k2,2n -t $'\t' | \
uniq | \
grep -vP "^\t$" | \
grep -v "x"   | \
grep -v ";" | \
grep -v "CHR_ID" > variant_list.txt

bcftools view -R variant_list.txt ~/TheGreatGenotyper/FineTunning/00-All.vcf.gz | \
bcftools norm -m +any > GWAS.vcf

bcftools annotate --force --rename-chrs chr_name_conv.txt GWAS.vcf | \
grep  -v "chrMT"> GWAS.chr.vcf

./addDummySample.sh  GWAS.chr.vcf GWAS.chr.dummy.vcf
```

## 2. Running Workflow
1. Edit config.yaml to configure input and output, as well as programs.

| Field             | Description                                                                                            |
|-------------------|--------------------------------------------------------------------------------------------------------|
| INPUT_DIR         | Input folder containing the sliced pangenome                                                           |
| PANGENOME         | Input Folder containing pangenome produced by pangenome genotyping [use cases](https://github.com/dib-lab/TheGreatGenotyper_usecases/tree/main/Pangenome_Genotyping/)
| TEMP_FOLDER       | Folder to store temporary files.                                                                       |
| INPUT_REFERENCE   | Input genome reference                                                                                 |
| INPUT_INDEX       | Txt file containing a list of CCDG indexes.                                                            |
| BEAGLE            | Binary Path for Beagle                                                                                 |
| BEAGLE_MAP        | Beagle MAP file. Download here([maps](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/)) |
| TheGreatGenotyper | Binary Path for The Great Genotyper                                                                    |
| bcftools          | Binary path for Bcftools                                                                               |
| OUTPUT_dir        | Output Directory                                                                                       |


2. Run the workflow
```
snakemake --configfile config.yaml -np
```
