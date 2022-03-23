# Operon Viz

## Steps

BaseDir: `s3://maf-users/Sunit_Jain/Kazuki/Gene_Neighborhood_Viz/v202111.1`

### Aggregate all the genomes

- Obtained most genomes from `s3://czbiohub-microbiome/ReferenceDBs/NinjaMap/Index/SCv2_4_20210212/fasta/`.
- Obtained the following 3 from `s3://czbiohub-microbiome/ReferenceDBs/NinjaMap/Hybrid_Closed_Genomes/`.
  - Holdemanella-biformis-DSM-3989-MAF-2
  - Solobacterium-moorei-DSM-22971-MAF-2
  - Anaerobutyricum-hallii-DSM-3353-MAF-2
- Used `Clostridium-hathewayi-DSM-13479-MAF-NJ35` from `s3://czbiohub-microbiome/ReferenceDBs/NinjaMap/Index/SCv2_4_20210212/fasta/` instead of the requested `Clostridium-hathewayi-DSM-13479-MAF-2`
- Upload everything to the BaseDir.

### Run BinQC on the genomes

```bash
aws batch submit-job \
    --profile maf \
    --job-name nf-binqc-SCv2_4_20210212 \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command=s3://nextflow-pipelines/nf-binqc,\
        "--fastas","s3://maf-users/Sunit_Jain/Kazuki/Gene_Neighborhood_Viz/v202111.1/genomes",\
        "--project","Kazuki_Gene_Neighborhood_Viz_v202111.1"
```

Output: `s3://genomics-workflow-core/Results/BinQC/Kazuki_Gene_Neighborhood_Viz_v202111.1/`

### Extract gff and gene files from the BinQC - GTDBtk outputs

```bash
aws s3 cp --recursive \
    s3://genomics-workflow-core/Results/BinQC/Kazuki_Gene_Neighborhood_Viz_v202111.1/04_GTDBtk/gtdbtk-results/identify/intermediate_results/marker_genes/ . \
    --exclude '*' --include '*_protein.fna' --include '*_protein.gff' --include '*_protein.faa'
```

### Blast the target sequence against proteins from each genome

BlastX with query as the antigen gene and database as the concatenated proteins file (`all_proteins.faa`) from all genomes.

(_If I were to do this again, I would verify that contigs in all genomes follow the standard naming convention. In this case the genomes obtained from the `Hybrid_Closed_Genomes` location did not follow convention._)

```bash
cat *_protein.faa > all_proteins.faa
bash -x makeblastdb.sh genomes/all_proteins.faa proteins prot
export TASK="blastx"; export QUERY_EXT=".fna";export NUM_ALN=100; bash -x run_blast.sh Blast/genome_db/queries/gene_of_interest.fna /mnt/efs/scratch/sunit/Kazuki/blastdb_custom/proteins
```

A cutoff of 70% identity was chosen to allow for a permissive search (Kazuki).

```bash
awk '$3>70 {print $0}' Blast/genome_db/results/gene_of_interest.blastx.outFmt_6.tsv | sort -k2,2 -k3,3nr | cut -f 2 >  antigen_hits.blastx_70.list
```

Based on a blastN performed earlier (not described here) where the antigen sequence was the database. We know that the following non-standard names map to these genomes:

```bash
Holdemanella-biformis-DSM-3989-MAF-2    1_2200
Solobacterium-moorei-DSM-22971-MAF-2    1_560
Anaerobutyricum-hallii-DSM-3353-MAF-2   1_519
```

Repeat this step for gene 2 and gene 3 in the cluster.

### Use the `gene_neighborhood_viz.py` script

Get the GFF for all genomes

```bash
aws s3 cp --recursive --profile maf \
    s3://genomics-workflow-core/Results/BinQC/Kazuki_Gene_Neighborhood_Viz_v202111.1/04_GTDBtk/gtdbtk-results/identify/intermediate_results/marker_genes/ . \
    --exclude '*' --include '*_protein.gff'
```

Create a tabular file with name of the genome and corresponding locus of interest. `pub/genomes_list.tsv`

Use this tabular file to draw the figure.

```bash
python gene_neighborhood_viz.py --seed pub/genomes_list.tsv --prefix pub/genomes_list
```
