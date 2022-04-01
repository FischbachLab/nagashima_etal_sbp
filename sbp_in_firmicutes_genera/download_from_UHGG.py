#!/usr/bin/env python3
# Given the metadata file
# Map the genome GFF file FTP URL to the acutal faa and fna files.
## Example:
## URL Example 1: ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0/all_genomes/MGYG0000015/MGYG000001531/genomes1/MGYG000001531.gff.gz
## URL Example 2:ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0/all_genomes/MGYG0000001/MGYG000000199/genomes1/MGYG000004930.gff.gz
## URL syntax: <URI_BASE>/all_genomes/<SPECIES_REP_PREFIX>/<SPECIES_REP>/genomes1/<IDENTIFIER>.gff.gz


import argparse
import gzip
import logging
import os
import shutil
import urllib.request as request
from contextlib import closing
from time import sleep

from numpy import random


def usage():
    parser = argparse.ArgumentParser(description="Download from UHGG")

    parser.add_argument("-m", "--metadata", help="UHGG metadata file")
    parser.add_argument("-o", "--outdir", help="output dir")

    return parser.parse_args().__dict__


def get_url_from_metadata(filename):
    with open(filename, "r") as file:
        for line in file:
            # skip header
            if line[0] == "#":
                continue

            columns = line.strip().split("\t")
            yield columns[-1]


def download_ftp_link(url, filepath, max_retries=0):
    success = False
    retries = 0
    while (not success) and (retries < max_retries):
        try:
            with closing(request.urlopen(url)) as r:
                with open(f"{filepath}", "wb") as f:
                    shutil.copyfileobj(r, f)
        except:
            retries += 1
            logging.warning(
                f"Retry [{retries}/{max_retries}]: Failed to download {url}, backing off ..."
            )
            sleeptimer = retries ** 2
            logging.warning(f"Will try again in {sleeptimer} seconds")
            sleep(sleeptimer)
        else:
            success = True
            break

    return success


def split_gff_fasta_file(filename, gff_dir, fna_dir):
    assert filename.endswith("gz") or filename.endswith(
        "gzip"
    ), f"File name must end with 'gz or gzip', provided name: {filename}"

    file_prefix = os.path.basename(filename).replace(".gff.gz", "")

    gff_file = os.path.join(gff_dir, f"{file_prefix}.gff.gz")
    fasta_file = os.path.join(fna_dir, f"{file_prefix}.fna.gz")

    with gzip.open(filename, "rb") as gzip_file:
        fh = gzip.open(gff_file, "wb")
        for line in gzip_file:
            if line.startswith(b"##FASTA"):
                fh.close()
                # switch to writing to a fasta file
                fh = gzip.open(fasta_file, "wb")
                continue

            fh.write(line)
        fh.close()

    return (gff_file, fasta_file)


def download_from_uhgg(mfile, outdir=None):
    subdirs = ["downloads", "gff", "fna"]
    if outdir is None:
        outdir = "."

    download_dir, gff_dir, fna_dir = [
        os.path.join(outdir, subdir) for subdir in subdirs
    ]

    _ = [
        os.makedirs(subdir, exist_ok=True)
        for subdir in [download_dir, gff_dir, fna_dir]
    ]

    err_counter = 0
    urls = get_url_from_metadata(mfile)
    for idx, url in enumerate(urls, 1):
        sleeptime = random.uniform(1, 3)
        filename = os.path.basename(url)
        filepath = os.path.join(download_dir, filename)
        if not os.path.exists(filepath):
            logging.info(f"#{idx} Downloading {filename} in {sleeptime}sec from {url}")
            sleep(sleeptime)
            if download_ftp_link(url, filepath, max_retries=5):
                _, _ = split_gff_fasta_file(filepath, gff_dir, fna_dir)
            else:
                err_counter += 1
                logging.error(f"#{idx} Error #{err_counter}: Failed to download {url}")
                if err_counter == 10:
                    break

    return


def main():
    args = usage()
    metadata_file = args["metadata"]
    outdir = args["outdir"].rstrip("/")
    download_from_uhgg(metadata_file, outdir)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s\t[%(levelname)s]:\t%(message)s",
    )
    main()
