#!/usr/bin/env python3
# (conda: ninjamap)

# Given a list of GFF files and their corresponding anchor gene
# visualize the neighborhood of the anchor gene upto X kb (default 25kb)
# for each gff file.
# Just one final figure is needed such that:
#  the anchor gene is aligned in the same direction 5’ to 3’ in all tracks.
#  the anchor gene is aligned at the start codon in all tracks.

import argparse
import colorsys
import logging
import os
import re
import urllib.parse
import random

from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import FeatureLocation, SeqFeature
from reportlab.lib.units import inch


def usage():
    parser = argparse.ArgumentParser(description="Gene Neighborhood Visualization")

    parser.add_argument("--seed", help="list of gff files with locus of interest")
    parser.add_argument("--prefix", help="specify the output prefix")

    return parser.parse_args().__dict__


#### from FlaGs ####
# Color generator
def random_color(c=None):
    """Generates a random color in RGB format."""
    if not c:
        c = random.random()
    d = 0.5
    e = 0.5
    return _hls2hex(c, d, e)


def _hls2hex(c, d, e):
    return "#%02x%02x%02x" % tuple(
        map(lambda f: int(f * 255), colorsys.hls_to_rgb(c, d, e))
    )


#### from FlaGs ####


def define_genome_track(d, track_name, track_num):
    feature_track = d.new_track(track_num, name=track_name, greytrack=False)
    return feature_track.new_set()


def add_feature_to_track(track, features_dict, track_metadata):
    f_start = features_dict["start"] - track_metadata["start"]
    f_end = features_dict["end"] - track_metadata["start"]
    f_strand = int(features_dict["strand"])
    f_name = features_dict["name"]
    f_label = features_dict["label"]
    f_color = features_dict["color"]

    # logging.info(f"Adding {f_name} at {f_start} - {f_end} on {f_strand} strand")
    feature = SeqFeature(FeatureLocation(f_start, f_end), strand=f_strand)
    track.add_feature(
        feature,
        name=f_name,
        label=f_label,
        sigil="ARROW",
        arrowshaft_height=1.0,
        color=f_color,
        # label_position="top",
        label_angle=0,
        label_size=6,
        height=5,
    )
    return


def parse_gff_record(line):
    segment = re.split("\t", line)
    attributesString = re.split("\n|;", segment[8])
    attributes = {}
    for s in attributesString:
        e = s.split("=")
        if len(e) >= 2:
            attributes.update({e[0]: urllib.parse.unquote(e[1])})
    if len(attributes) == 0:
        return None

    gene_suffix = attributes["ID"].split("_")[1]
    record = {
        "name": f"{segment[0]}_{gene_suffix}",
        "contig": segment[0],
        "feature": segment[2],
        "start": int(segment[3]) - 1,
        "end": int(segment[4]),
        "strand": segment[6],
    }
    return record


def read_gff_file(gff_file):
    with open(gff_file) as gff:
        for line in gff.readlines():
            if re.match("##FASTA", line):
                break
            if line[0] == "#":
                continue
            feature = parse_gff_record(line)

            if (feature["strand"] != "+") and (feature["strand"] != "-"):
                continue

            if feature is None:
                continue

            yield feature


def get_features_of_interest(
    gff_file,
    gene_of_interest,
    color_list,
    gene_names_list,
    window_size=None,
    cluster_size=None,
):
    # returns a list of dictionaries
    # [{
    #     name : str
    #     start : int,
    #     end : int,
    #     strand : +1/-1,
    # }]
    list_of_records = list()
    genes_in_cluster = cluster_size

    # First pass through the gff to identify the location of the gene of interest
    # and determine the window start and end coordinates
    feature_locii = set()
    features = read_gff_file(gff_file)
    contig_w_gene_of_interest = ""
    positions = list()
    past_features = dict()
    operon_strand = "0"
    for f, feature in enumerate(features):
        past_features.update({f: (feature["name"], feature["start"], feature["end"])})
        if feature["name"] == gene_of_interest:
            # if strand is "+" get next 2 from file
            # if strand is "-" get previous 2 from past_features
            operon_strand = feature["strand"]
            feature_locii.add(feature["name"])
            positions.append(feature["start"])
            positions.append(feature["end"])
            contig_w_gene_of_interest = feature["contig"]
            cluster_size = cluster_size - 1
            while cluster_size > 0:
                # If the selected gene is on a negative strand,
                # add the previous genes to the selected operon set.
                # else continue on to the next one.
                if operon_strand == "-":
                    name, start, end = past_features[(f - cluster_size)]
                elif operon_strand == "+":
                    try:
                        feature = next(features)
                        name = feature["name"]
                        start = feature["start"]
                        end = feature["end"]

                    except StopIteration as s:
                        logging.error(f"{feature['name']} \t {cluster_size}")
                        raise s

                feature_locii.add(name)
                positions.append(start)
                positions.append(end)
                cluster_size = cluster_size - 1

        if len(feature_locii) == genes_in_cluster:
            break

    start_pos_of_interest = min(positions)
    end_pos_of_interest = max(positions)

    assert (
        start_pos_of_interest != end_pos_of_interest
    ), f"start and end positions are the same: {end_pos_of_interest}"

    logging.info(
        f"Found {len(feature_locii)} gene(s) of interest at {start_pos_of_interest} - {end_pos_of_interest}"
    )

    window_start = int((start_pos_of_interest + end_pos_of_interest) / 2) - int(
        window_size / 2
    )
    window_end = int((start_pos_of_interest + end_pos_of_interest) / 2) + int(
        window_size / 2
    )

    # If the gene is on a negative strand,
    # reverse the color list so the color order is maintained.
    # when the track and features are adjusted to show everything in the same direction
    if operon_strand == "-":
        color_list = color_list[::-1]
        gene_names_list = gene_names_list[::-1]

    logging.info(f"Searching for neighbors between {window_start} - {window_end}")

    # Second pass through the GFF to aggregate all the genes within the indentified window
    i = 0
    features2 = read_gff_file(gff_file)
    for feature in features2:
        if feature["contig"] != contig_w_gene_of_interest:
            continue

        this_strand = feature["strand"]
        # If the feature is on the - strand, flip the feature strand
        # so the figure shows all tracks in the same direction
        if operon_strand == "-":
            this_strand = flip_strand(feature["strand"])

        if (feature["start"] >= window_start) and (feature["end"] <= window_end):
            should_label = False
            feature_color = "grey"
            feature_name = feature["name"]
            if feature["name"] in feature_locii:
                feature_color = color_list[i]
                feature_name = gene_names_list[i]
                should_label = True
                i += 1

            feature_start = feature["start"]
            feature_end = feature["end"]
            # If the feature is on the - strand, flip the feature, start, stop
            # so the figure shows all tracks in the same direction
            if operon_strand == "-":
                feature_start = -1 * feature["start"]
                feature_end = -1 * feature["end"]

            feature_start, feature_end = sorted([feature_start, feature_end])
            list_of_records.append(
                {
                    "name": feature_name,
                    "start": feature_start,
                    "end": feature_end,
                    "strand": "+1" if this_strand == "+" else "-1",
                    "color": feature_color,
                    "label": should_label,
                }
            )
        elif feature["end"] > window_end:
            break

    # If the feature is on the - strand, flip the track
    # so the figure shows all tracks in the same direction
    if operon_strand == "-":
        window_start, window_end = sorted([-window_start, -window_end])

    track_metadata = {
        "contig": "",
        "start": window_start,
        "end": window_end,
    }

    return track_metadata, list_of_records


def flip_strand(item):
    if item == "+":
        return "-"
    elif item == "-":
        return "+"


def main():
    args = usage()

    gff_list = args["seed"]
    prefix = args["prefix"]
    window_size = 20000
    cluster_size = 3

    diagram = GenomeDiagram.Diagram(f"{prefix}")
    num_genomes = 0
    # TODO @sunitj: make the random color generator generate cluster_size 'distinct' colors
    # color_list = [random_color() for i in range(cluster_size)]
    color_list = ["#bf3f53", "#3f85bf", "#bdbf3f"]  # red, blue, pale green
    gene_names_list = ["T-Cell antigen, SBP", "Gene B", "Gene C"]
    logging.info(f"Color pallette: {color_list}")
    with open(gff_list, "r") as gff_h:
        for idx, line in enumerate(gff_h):
            if line[0] == "#":
                continue
            gff, locus = line.strip().split("\t")
            logging.info(f"Creating new track for {gff}")
            track_name = os.path.basename(gff)
            track = define_genome_track(diagram, track_name, idx)
            track_metadata, feature_dict_list = get_features_of_interest(
                gff,
                locus,
                color_list,
                gene_names_list,
                window_size=window_size,
                cluster_size=cluster_size,
            )

            logging.info(f"\tAdding {len(feature_dict_list)} features.")
            _ = [
                add_feature_to_track(track, feature_dict, track_metadata)
                for feature_dict in feature_dict_list
            ]
            num_genomes = idx
    diagram.draw(
        format="linear",
        pagesize=(30 * inch, num_genomes * inch),
        fragments=1
        # start=track_metadata["start"] - buffer,
        # end=track_metadata["end"] + buffer,
        # track_size=window_size + buffer,
    )
    diagram.write(f"{prefix}.gene_neighborhood.pdf", "pdf")


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.DEBUG, format="%(asctime)s\t[%(levelname)s]:\t%(message)s",
    )
    main()
