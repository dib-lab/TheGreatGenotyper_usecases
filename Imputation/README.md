# Imputation Using 4K reference panel

## 1. Creating pseudo-microarray variant calls for HG002

I am adopting an experiment done in this [paper](https://academic.oup.com/bioinformatics/article/36/24/5582/6064144)

1. Download Infinium OmniExpress-24 Kit product files from [here](https://support.illumina.com/array/array_kits/humanomniexpress-24-beadchip-kit/downloads.html)

2. Convert The manifest file to vcf file using [GTCtoVCF](https://github.com/Illumina/GTCtoVCF)
```
git clone https://github.com/Illumina/GTCtoVCF.git
conda create -c miniconda -c  bioconda  -y -n GTC_TO_VCF python=2.7 numpy=1.11.2 pyvcf=0.6.8 pysam=0.9.0

python gtc_to_vcf.py --manifest-file InfiniumOmniExpress-24v1-4_A1.csv  \
                     --genome-fasta-file hs37d5.fa \
                     --output-vcf-path InfiniumOmniExpress-24v1-4_A1.vcf \
                     --log-file InfiniumOmniExpress-24v1-4_A1.log
```

3. Liftover the vcf file from GRCh37 to GRCh38 using [picard](https://broadinstitute.github.io/picard/)
```
wget https://hgdownload.soe.ucsc.edu/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz

for i in {1..22} X Y MT; do echo "$i chr$i" >> chr_name_conv.txt ; done
bcftools annotate --force \
                  --rename-chrs chr_name_conv.txt \
                  InfiniumOmniExpress-24v1-4_A1.vcf > InfiniumOmniExpress-24v1-4_A1.renamed_chr.vcf

picard -Xmx20g LiftoverVcf I=InfiniumOmniExpress-24v1-4_A1.renamed_chr.vcf \
                           O=InfiniumOmniExpress-24v1-4_A1.GRCh38.vcf \
                           CHAIN=hg19ToHg38.over.chain  \
                           REJECT=rejected_variants.vcf \
                           R=GRCh38_full_analysis_set_plus_decoy_hla.fa

```

4. Create the pseudo microarray variant call by intersecting InfiniumOmniExpress-24v1-4_A1.GRCh38.vcf with the truth small variant set for HG002 from [GIAB](https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/).

```
wget 
https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz

bcftools merge -0 InfiniumOmniExpress-24v1-4_A1.GRCh38.vcf.gz HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz | \
bcftools view -e 'ID="."'  | \
sed -e 's/\.\t\.\t\./.\tGT\t0\/0/' | \
grep -P "\t\.$"| \
bgzip > HG002_GRCh38.InfiniumOmniExpress-24v1-4_A1.vcf.gz

tabix -p vcf HG002_GRCh38.InfiniumOmniExpress-24v1-4_A1.vcf.gz
```


## 2. Imputation


1. Download 4k reference panel from [here](https://farm.cse.ucdavis.edu/~mshokrof/4k_reference_panel/)

```
wget https://farm.cse.ucdavis.edu/~mshokrof/4k_reference_panel/chr20.tagged.vcf.gz
wget https://farm.cse.ucdavis.edu/~mshokrof/4k_reference_panel/chr20.tagged.vcf.gz.tbi
```

2. Impute Chr20 only
```
bcftools view HG002_GRCh38.InfiniumOmniExpress-24v1-4_A1.vcf.gz chr20 >  HG002_GRCh38.InfiniumOmniExpress-24v1-4_A1.chr20.vcf

java -Xmx20G -jar beagle.22Jul22.46e.jar gp=true \
                                         gt=HG002_GRCh38.InfiniumOmniExpress-24v1-4_A1.chr20.vcf \
                                         ref=$chr20.tagged.vcf.gz \
                                         out=tmp \
                                         nthreads=8  \
                                         map=plink.autsomal.map

bcftools reheader -h  header  tmp.vcf.gz | \
    bcftools norm -m -any -Ou | \
    bcftools +setGT  -Ou -- -t . -n . | \
    bcftools +fill-tags -Ou -- -t AC | \
    bcftools view -e 'AC=0 ' -Oz -o  HG002_imputed.vcf.gz

tabix -p vcf HG002_imputed.vcf.gz
rm tmp.vcf.gz
```

## 3. Benchmarking

### Benchmark small variants using truth set from [GIAB NISTv4.2.1](https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/).

```
wget  https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed

zcat HG002_imputed.vcf.gz | \
python3 assign-variant-len.py | \
bcftools view -e 'VARLEN>50' -  | \
bgzip > out3.vcf.gz

singularity exec   -B "${PWD}":"/data"  \
                   -B "${PWD}/reference":"/reference" \
                   docker://jmcdani20/hap.py:v0.3.12 \
                   /opt/hap.py/bin/hap.py   \
                   /data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz   \
                   /data/out3.vcf.gz   \
                   -f /data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
                   -r /reference/GRCh38_full_analysis_set_plus_decoy_hla.fa   \
                   -o /data/happy_pg_chr20  --engine=vcfeval    \
                   -l chr20 --threads 16 
```

### Benchmark SV using truth set from [NIST0.6](https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6/README_SV_v0.6.txt)

1. Download Liftover benchmark data to GRCH38
```
wget https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz.tbi
wget https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6/HG002_SVs_Tier1_v0.6.bed

sed -e 's/^/chr/' HG002_SVs_Tier1_v0.6.bed  > HG002_SVs_Tier1_v0.6.chr.bed
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
./liftOver HG002_SVs_Tier1_v0.6.chr.bed hg19ToHg38.over.chain HG002_SVs_Tier1_v0.6.GRCh38.bed unmapped

grep -P "^chr20" HG002_SVs_Tier1_v0.6.GRCh38.bed > HG002_SVs_Tier1_v0.6.GRCh38.chr20.bed

bcftools annotate --force \
                  --rename-chrs chr_name_conv.txt \
                  HG002_SVs_Tier1_v0.6.vcf.gz  | \
bgzip  >  HG002_SVs_Tier1_v0.6.chrrenamed.vcf.gz

zcat HG002_SVs_Tier1_v0.6.chrrenamed.vcf.gz | \
sed -e 's/END=[^;]*;//'   | \
bgzip > HG002_SVs_Tier1_v0.6.chrrenamed.noEND.vcf.gz

picard -Xmx20g LiftoverVcf I=HG002_SVs_Tier1_v0.6.chrrenamed.noEND.vcf.gz
                           O=HG002_SVs_Tier1_GRCH38.vcf.gz
                           CHAIN=hg19ToHg38.over.chain
                           REJECT=rejected_variants.vcf
                           R=GRCh38_full_analysis_set_plus_decoy_hla.fa
```

2. Benchmark
```
bcftools view  HG002_imputed.vcf.gz | grep -v "0|0"  |bgzip > out2.vcf.gz
bcftools reheader -h  header out2.vcf.gz > out3.vcf.gz
bcftools norm -m -any out3.vcf.gz |grep -v "0|0"  > out.4.vcf
python  unphase.py out.4.vcf  out.5.vcf.gz
tabix -p vcf out.5.vcf.gz

truvari bench -b HG002_SVs_Tier1_GRCH38.vcf.gz  \
              -c out.5.vcf.gz  \
              --includebed HG002_SVs_Tier1_v0.6.GRCh38.chr20.bed
              -o truvari_output \
              --passonly   \
              -r 2000 -C 3000  \
              --reference GRCh38_full_analysis_set_plus_decoy_hla.fa
```
