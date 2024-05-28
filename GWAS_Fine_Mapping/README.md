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
