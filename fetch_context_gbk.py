#!/usr/bin/env python

"""
Scans a given GenBank file and searches in the CDS features for a specified
label (as gene id, locus_tag, etc.) and extracts the locus around that feature
"""

from pathlib import Path
import sys
import os
import argparse
from Bio import SeqIO

__author__ = "Jorge Navarro"
__version__ = "1.0"
__maintainer__ = "Jorge Navarro"
__email__ = "jorge.c.navarro.munoz@gmail.com"

def parameter_parser():
    def_extra = 20
    def_o = Path("./output")

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="GenBank file", type=Path, required=True)
    parser.add_argument("-l", "--label", help="Target label (gene_id, locus_tag, etc.)",
        type=str, required=True)
    parser.add_argument("-e", "--extra", help=f"Number of kbps at either side of \
        the target gene/protein/locus_tag to download. Default={def_extra}", \
        default=def_extra, type=int)
    parser.add_argument("-o", "--outputfolder", help=f"Where to put retrieved \
        files. Default={def_o}", type=Path, default=def_o)

    return parser.parse_args()


# throw here as many qualifiers as desired
def find_label(qualifiers, target):
    if "gene" in qualifiers:
        if qualifiers["gene"][0] == target:
            return True
    if "protein_id" in qualifiers:
        if qualifiers["protein_id"][0] == target:
            return True
    if "proteinId" in qualifiers:
        if qualifiers["proteinId"][0] == target:
            return True
    if "locus_tag" in qualifiers:
        if qualifiers["locus_tag"][0] == target:
            return True

    if "name" in qualifiers:
        if qualifiers["name"][0] == target:
            return True

    return False


# lazy way of getting new borders that include features. Too many comparisons
def check_borders_left(spans, pos):
    for span in spans:
        if pos in span:
            return span[0] - 10
        # if pos is in the middle of two feature spans, return orig. pos
        elif span[0] > pos:
                return pos
    return pos



def check_borders_right(spans, pos):
    for span in spans:
        if pos in span:
            return span[-1] + 10
        elif span[0] < pos:
            return pos
    # reached the and and there were not annotations. Return orig pos
    return pos



def scan_and_extract(gbk, target, extra, o):
    num_extraction = 1 # just in case there is more than 1 feature with target label
    out_filename_base = f"{gbk.stem}_{target}_"
    target_found = False

    try:
        records = list(SeqIO.parse(str(gbk), "genbank"))
    except ValueError as e:
        print("Error, not able to parse file {}: {}".format(str(gbk), str(e)))
    else:
        for record in records:
            features = [f for f in record.features if f.type == "CDS"]
            spans = [range(f.location.start, f.location.end) for f in features]

            for feature in features:
                if feature.type == "CDS":
                    if find_label(feature.qualifiers, target):
                        target_found = True
                        start_feature = int(feature.location.start)
                        end_feature = int(feature.location.end)

                        start_extraction = max(0, start_feature - 1000*extra)
                        end_extraction = min(end_feature + 1000*extra, len(record))

                        # check if we have features at borders. If so, include them
                        start_extraction = check_borders_left(spans, start_extraction)
                        end_extraction = check_borders_right(spans, end_extraction)

                        extraction = record[start_extraction:end_extraction]
                        rc = ""
                        if feature.location.strand != 1:
                            # extract and reverse
                            rc = "_rc"
                            extraction = extraction.reverse_complement(id=f"{extraction.id}{rc}", annotations=True)

                        with open(o / f"{out_filename_base}{num_extraction}{rc}.gbk", "w") as ef:
                            SeqIO.write(extraction, ef, "genbank")
                            num_extraction += 1



    if not target_found:
        print("Finished, but nothing found with target label...")
    else:
        print(f"Finished, {num_extraction-1} extraction(s) done")

    return


if __name__ == '__main__':
    args = parameter_parser()

    if not args.input.is_file():
        sys.exit(f"Error, {args.input} not a file")

    o = args.outputfolder
    if not o.is_dir():
        os.makedirs(o, exist_ok=True)

    print(f"Attempting to extract locus around feature with label '{args.label}'")

    scan_and_extract(args.input, args.label, args.extra, o)
