aws batch submit-job \
    --job-name nf-binqc-hCom1d_hCom2d \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command=FischbachLab/nf-binqc,\
"--fastas","s3://maf-versioned/ninjamap/Index/Kazuki_hCom1d_hCom2d/fasta",\
"--project","Kazuki_hCom1d_hCom2d"