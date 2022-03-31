#!/usr/bin/env python3
# Given the metadata file
# Map the genome GFF file FTP URL to the acutal faa and fna files.
## Example:
## URL in file: ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0/all_genomes/MGYG0000015/MGYG000001531/genomes1/MGYG000001531.gff.gz
## Convert to:  ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0/species_catalogue/MGYG0000015/MGYG000001531/genome/MGYG000001531.faa
## All available file types:
# MGYG000001531.faa
# MGYG000001531.fna
# MGYG000001531.fna.fai
# MGYG000001531.gff
# MGYG000001531_InterProScan.tsv
# MGYG000001531_annotation_coverage.tsv
# MGYG000001531_cazy_summary.tsv
# MGYG000001531_cog_summary.tsv
# MGYG000001531_eggNOG.tsv
# MGYG000001531_kegg_classes.tsv
# MGYG000001531_kegg_modules.tsv
# MGYG000001531_rRNAs.fasta
import os
import logging
import argparse
import shutil
import urllib.request as request
from contextlib import closing
from numpy import random
from time import sleep

URI_BASE = (
    "ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0"
)


def usage():
    parser = argparse.ArgumentParser(description="Download from UHGG")

    parser.add_argument("-m", "--metadata", help="UHGG metadata file")
    parser.add_argument(
        "-e",
        "--extension",
        help="portion of the filename that is NOT the identifier.",
        default=".faa",
    )
    parser.add_argument("-o", "--outdir", help="output dir")

    return parser.parse_args().__dict__


def get_url_from_metadata(filename, ext, uri_base=URI_BASE):
    with open(filename, "r") as file:
        for line in file:
            # skip header
            if line[0] == "#":
                continue

            columns = line.strip().split("\t")
            identifier = columns[0]
            identifier_base = identifier[:-2]
            new_uri = f"{uri_base}/species_catalogue/{identifier_base}/{identifier}/genome/{identifier}{ext}"

            yield new_uri


def download_from_uhgg(mfile, ext, outdir=None):
    if outdir is None:
        outdir = "."
    else:
        os.makedirs(outdir, exist_ok=True)

    urls = get_url_from_metadata(mfile, ext)
    err_counter = 0
    for url in urls:
        sleeptime = random.uniform(2, 4)
        filename = os.path.basename(url)
        filepath = f"{outdir}/{filename}"
        if not os.path.exists(filepath):
            logging.info(f"Downloading {filename} in {sleeptime}sec from {url}")
            sleep(sleeptime)
            try:
                with closing(request.urlopen(url)) as r:
                    with open(f"{filepath}", "wb") as f:
                        shutil.copyfileobj(r, f)
            except:
                err_counter += 1
                logging.error(f"Error #{err_counter}: Failed to download {url}")
                if err_counter == 10:
                    break
    return


def main():
    args = usage()
    metadata_file = args["metadata"]
    extension = args["extension"]
    outdir = args["outdir"].rstrip("/")
    download_from_uhgg(metadata_file, extension, outdir)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s\t[%(levelname)s]:\t%(message)s",
    )
    main()
