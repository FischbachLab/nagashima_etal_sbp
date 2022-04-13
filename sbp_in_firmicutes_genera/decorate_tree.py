#!/usr/bin/env python3

# Accept a tree file, a stats file and optionally a map file
# if map file is present,
## aggregate the stats file such that
## each leaf of the tree corresponds to exactly one row in the
## aggregated stats file.
# Construct a tree with the additional tracks showing the stats


import logging
from ete3 import PhyloTree, faces, AttrFace, TreeStyle, NodeStyle
import pandas as pd
import numpy as np

def usage():
    

def main():
    args = usage()

if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s\t[%(levelname)s]:\t%(message)s",
    )
    main()