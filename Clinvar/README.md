# Population Genotyping of the Clinvar database
1. Download and preprocess the Clinvar database:
```
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz

```
2. run the snakemake workflow

```
snakemake --configfile config.GRCh38.yaml
```

The final output will be at "{outputfolder}/genotyped/clinvar.tagged.vcf.gz". The final output can be downloaded from [vcf](https://farm.cse.ucdavis.edu/~mshokrof/The_great_genotyper_clinvar/)

