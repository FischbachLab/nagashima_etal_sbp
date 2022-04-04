# Steps to run this analysis

## Run GTDBtk on all the genomes in hCom1 and hCom2

```bash
aws batch submit-job \
    --job-name nf-binqc-hCom1d_hCom2d \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command=FischbachLab/nf-binqc,\
"--fastas","s3://maf-versioned/ninjamap/Index/Kazuki_hCom1d_hCom2d/fasta",\
"--project","Kazuki_hCom1d_hCom2d"
```

## Download the .tree and .summary.tsv outputs

## Create the `formatted_names.tsv`

This was a manual process.

## Prune the gtdbtk tree to the shortlist in the formatted names

```bash
python treeViz.py \
    data/gtdb.Kazuki_hCom1d_hCom2d.bac120.classify.tree \
    data/gtdb.Kazuki_hCom1d_hCom2d.bac120.summary.tsv \
    pub/hCom1d_hCom2d \
    pub/formatted_names.tsv
```
