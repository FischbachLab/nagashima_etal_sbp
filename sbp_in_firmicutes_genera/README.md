# SBP in Firmicutes genera

Here are the steps for this analysis.

## UHGG Metadata

We obtained the metadata for UHGG database using the following link on March 25, 2022

```bash
https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0/genomes-all_metadata.tsv
```

## UHGG Isolates

For the purposes of this analysis we focused on isolate genomes

```bash
head -n 1 genomes-all_metadata.tsv > UHGG.isolates_metadata.tsv
awk '$2 == "Isolate" {print $0}' genomes-all_metadata.tsv >>UHGG.isolates_metadata.tsv
```

### Genera in firmicutes with >= 10 genomes

```bash
cut -f 15 UHGG.isolates_metadata.tsv | grep "p__Firmicutes" | cut -f 6 -d ';' | sort | uniq -c | awk '{printf "%s\t%s\n",$2,$1}' | sort -k2,2nr | awk '$2>=10' | cut -f 1 > UHGG.isolates.selected_firmicutes_genera.list
wc -l UHGG.isolates.selected_firmicutes_genera.list 
# found 60 genera
```

### How many species in the 60 genera?

```bash
cut -f 15 UHGG.isolates_metadata.tsv | grep "p__Firmicutes" | cut -f 6 -d ';' | sort | uniq -c | awk '{printf "%s\t%s\n",$2,$1}' | sort -k2,2nr | awk '$2>=10' | awk '{sum+=$2;} END{print sum;}'
# 2174 genomes
```

## Extract metadata for the shortlisted `2174` genomes

```bash
grep -Fwf UHGG.isolates.selected_firmicutes_genera.list UHGG.isolates_metadata.tsv > UHGG.isolates.selected_firmicutes_genera.metadata.tsv
```

## Download the genomes

```bash
python download_from_UHGG.py --outdir firmicutes_genera.gte10 --metadata UHGG.isolates.selected_firmicutes_genera.metadata.tsv
```

This script will download the `gff.gz` file from the `FTP_download` column in the metadata and split the file to separate the contig sequences from the gff.

## Extract genes from the contigs using the corresponding GFF files

Use bedtools

something like (not sure I need the `-s`, since blast should work regardless.)

```bash
bedtools getfasta -fi contigs.fna -bed contigs.gff -s -fullHeader
```

## Create a blastDB with the 3 genes in the operon

## BlastX with operon DB against all the genes files

## For each genera calculate how many genomes had hits to the entire operon

## Build a tree with a bar chart depicting how many genomes in the genera had a hit
