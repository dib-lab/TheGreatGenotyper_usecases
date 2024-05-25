# Pangenome Genotyping

1. Downlaod Decomposed HPRC from https://zenodo.org/records/6797328

```
wget https://zenodo.org/records/6797328/files/cactus_filtered_ids.vcf.gz
```

2. Download and Build [TheGreatGenotyper](https://github.com/dib-lab/TheGreatGenotyper), [beagle](https://faculty.washington.edu/browning/beagle/beagle.html), [bcftools](https://samtools.github.io/bcftools/howtos/install.html) and [AnnotSV](https://lbgi.fr/AnnotSV/).


3. Split the pangenome into 10 slices for faster computation. The following script splits variants from each chromosome into 10 chunks and forms a slice using a chunk from each chromosome.

```
./slice_pangenome.sh  cactus_filtered_ids.vcf.gz 10 sliced_pangenome_10/
```

4. Edit config.yaml to configure input and output, as well as programs.

| Field             | Description                                                                                            |
|-------------------|--------------------------------------------------------------------------------------------------------|
| INPUT_DIR         | Input folder containing the sliced pangenome                                                           |
| TEMP_FOLDER       | Folder to store temporary files.                                                                       |
| INPUT_REFERENCE   | Input genome reference                                                                                 |
| INPUT_INDEX       | Txt file containing a list of CCDG indexes.                                                            |
| BEAGLE            | Binary Path for Beagle                                                                                 |
| BEAGLE_MAP        | Beagle MAP file. Download here([maps](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/)) |
| TheGreatGenotyper | Binary Path for The Great Genotyper                                                                    |
| AnnotSV           | Binary Path for AnnotSV                                                                                |
| bcftools          | Binary path for Bcftools                                                                               |
| OUTPUT_dir        | Output Directory                                                                                       |


5. Run the workflow
```
snakemake --configfile config.yaml -np
```

