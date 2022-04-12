#!/bin/bash -x
# shellcheck disable=SC2086
set -e
set -u
set -o pipefail

TASK=${TASK:-"blastn"}
QUERY_EXT=${QUERY_EXT:-".fna"}
THREADS=${THREADS:-6}
MAX_ALN=${MAX_ALN:-5}

QUERY_PATH=${1}
LOCAL_DB=${2:-"/mnt/efs/databases/Blast/nt/db/nt"}
RESULTS_PATH=$(echo "${QUERY_PATH}" | sed -e "s/queries/results/" -e "s/${QUERY_EXT}/.${TASK}/")

BASE_PATH=$(pwd)
LOCAL_QUERY_FASTA="${BASE_PATH}/${QUERY_PATH}"
LOCAL_RESULT_FILE="${BASE_PATH}/${RESULTS_PATH}"

LOCAL_DB_DIR="$(dirname ${LOCAL_DB})"
DBNAME="$(basename ${LOCAL_DB})"

LOCAL_QUERY_DIR="$(dirname ${LOCAL_QUERY_FASTA})"
QUERY_FILE_NAME="$(basename ${LOCAL_QUERY_FASTA})"

LOCAL_RESULTS_DIR="$(dirname ${LOCAL_RESULT_FILE})"
RESULTS_FILE_NAME="$(basename ${LOCAL_RESULT_FILE})"

mkdir -p ${LOCAL_RESULTS_DIR}
DOCKER_DB_DIR="/blast/blastdb_custom"
DOCKER_QUERIES_DIR="/blast/queries"
DOCKER_RESULTS_DIR="/blast/results"

# docker run --rm \
#     -v "${BASE_PATH}:${BASE_PATH}":rw \
#     -w /blast/blastdb_custom \
#     ncbi/blast ${TASK} -version > ${TASK}.version.txt

docker container run --rm \
    -v ${LOCAL_DB_DIR}:${DOCKER_DB_DIR}:ro \
    -v ${LOCAL_QUERY_DIR}:${DOCKER_QUERIES_DIR}:ro \
    -v ${LOCAL_RESULTS_DIR}:${DOCKER_RESULTS_DIR}:rw \
    ncbi/blast:latest \
        ${TASK} \
            -num_threads ${THREADS} \
            -query ${DOCKER_QUERIES_DIR}/${QUERY_FILE_NAME} \
            -db ${DBNAME} \
            -dbsize 1000000 \
            -num_alignments ${MAX_ALN} \
            -outfmt '6 std qlen slen qcovs ppos' \
            -out ${DOCKER_RESULTS_DIR}/${RESULTS_FILE_NAME}.outFmt_6.tsv