# SBP in Firmicutes genera

Here are the steps for this analysis.

## Data Location

```bash
s3://maf-versioned/nagashima_etal_SBP/genus_tree/
```

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

### Generate a Genome-Genus map file for annotating the tree

```bash
cut -f 1,15 UHGG.isolates.selected_firmicutes_genera.metadata.tsv | sed "s/;/\t/g" | cut -f 1,7 | sed "s/g__//" > firmicutes_genome_genera_map.tsv
```

## Download the genomes

```bash
python download_from_UHGG.py \
    --outdir firmicutes_genera.gte10 \
    --metadata UHGG.isolates.selected_firmicutes_genera.metadata.tsv
```

This script will download the `gff.gz` file from the `FTP_download` column in the metadata and split the file to separate the contig sequences from the gff.

## Extract genes from the contigs using the corresponding GFF files

Use the `run_bedtools_getfasta.sh` script to extract the genes from the genomes.

```bash
ls fna/*.gz | \
    parallel \
        -j 4 \
        --joblog run_bedtools_getfasta.joblog \
        --plus \
        "bash run_bedtools_getfasta.sh {} gff/{/..}.gff.gz" &> run_bedtools_getfasta.log
```

### Build a genus tree from 1 species rep in each genus

Get one genome for each genera:

```bash
cat UHGG.isolates.selected_firmicutes_genera.list \
    | parallel -j 4 \
    "grep -m 1 -Fw {} UHGG.isolates.selected_firmicutes_genera.metadata.tsv | cut -f 14" \
    > UHGG.isolates.selected_firmicutes_genera.sampled.list
```

Make sure we have the same number of genomes as genus:

```bash
wc -l UHGG.isolates.selected_firmicutes_genera.sampled.list
# 60 UHGG.isolates.selected_firmicutes_genera.sampled.list
```

Run BinQC to place these genomes on the GTDBtk tree

```bash
aws batch submit-job \
    --job-name nf-binqc-kazuki-firmicutes \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command=FischbachLab/nf-binqc,\
"--fastas","s3://maf-versioned/nagashima_etal_SBP/genus_tree/firmicutes_genera.gte10/genera_species_reps_fna",\
"--project","Kazuki_firmicutes_genera_tree"
```

Create a mapping file such that each row is tab seperated and contains two columns the species rep and it's genus.

```bash
echo -e "user_genome\tgenus" > selected_firmicutes_genera.formatted_names.tsv
grep -Fwf UHGG.isolates.selected_firmicutes_genera.sampled.list firmicutes_genome_genera_map.tsv >> selected_firmicutes_genera.formatted_names.tsv
head selected_firmicutes_genera.formatted_names.tsv
# user_genome     genus
# MGYG000000002   Blautia_A
# MGYG000000004   Anaerotruncus
# MGYG000000006   Staphylococcus
# MGYG000000007   Lactobacillus
# MGYG000000009   Ligilactobacillus
# MGYG000000014   Clostridium
# MGYG000000018   Coprococcus
# MGYG000000021   Limosilactobacillus
# MGYG000000022   Faecalibacterium
```

```bash
python ../treeViz.py gtdb.Kazuki_firmicutes_genera_tree.bac120.classify.tree gtdb.Kazuki_firmicutes_genera_tree.bac120.summary.tsv firmicutes_genera selected_firmicutes_genera.formatted_names.tsv
```

## Create a nucleotide blastDB with the SBP gene

### Extract the genes of interest

```bash
bash -x run_bedtools_getfasta.sh Clostridium-nexile-DSM-1787-MAF-2.fna.gz genes_of_interest.gff
```

### Make blastDB

```bash
bash -x makeblastdb.sh ../genes_of_interest/genes/Clostridium-nexile-DSM-1787-MAF-2.genes.fna sbp_operon_db
```

## BlastN with SBP DB against all the genes files

```bash
mv genes gene_queries
ls gene_queries/*fna \
    | parallel -j 5 \
    --joblog run_blastn.sbp_operon_db.joblog \
    "bash run_blast.sh {} $(pwd)/blastdb_custom/sbp_operon_db" &> run_blastn.sbp_operon_db.log
```

## For each genera calculate how many genomes had hits to the entire operon

### Do any genomes contain the entire operon?

```bash
cat gene_results/*tsv \
    | awk '$3 >= 70 {printf "%s\t%s\n", $1, $2}' \
    | sort -u \
    | sed "s/:/\t/" \
    | cut -f 1 \
    | sort \
    | uniq -c \
    | awk '$1 >= 3 {printf "%s\t%s\n",$2, $1}'
# 16
```

### Which Genera do they belong to?

## Build a tree with a bar chart depicting how many genomes in the genera had a hit

## Protein database

 bash -x makeblastdb.sh ../genes_of_interest/proteins/Clostridium-nexile-DSM-1787-MAF-2_protein.faa sbp_operon_db prot
