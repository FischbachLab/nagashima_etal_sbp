#!/usr/bin/bash -x

FNA="${1}"
GFF="${2}"
UNZIPPED_FNA=$(echo "${FNA}" | sed "s/.gz//")
OUT=$(basename "${UNZIPPED_FNA}" | sed "s/.fna//").genes.fna

OUTDIR="genes"
CONTAINER="quay.io/biocontainers/bedtools"
CONTAINER_VERSION="2.23.0--h5b5514e_6"

mkdir -p ${OUTDIR}

gunzip ${FNA}

docker container run --rm \
    -v $(pwd):$(pwd) \
    -w $(pwd) \
    ${CONTAINER}:${CONTAINER_VERSION} \
    bedtools \
        getfasta \
            -fullHeader \
            -fi ${UNZIPPED_FNA} \
            -bed ${GFF} \
            -fo ${OUTDIR}/${OUT}

gzip "${UNZIPPED_FNA}"