# Population Genotyping of the Clinvar database
1. Download and preprocess the Clinvar database:
```
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz

```
2. run the snakemake workflow

```
snakemake --configfile config.GRCh38.yaml
```

Final output will be at "{outputfolder}/genotyped/clinvar.tagged.vcf.gz"
