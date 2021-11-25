# Fetch genomic context from a GenBank file

## What it does

This script reads a GenBank file and, given a label, it will try to find a matching qualifier. If found, it will extract it, plus some extra kbps upstream/downstream (default: 20 kbps)


## Notes

- The label must be complete (no fuzzy search)
- The script only looks at CDS features
- Qualifiers where the label will be searched for: `gene`, `protein_id`, `proteinId`, `locus_tag` and `name`
- If the extension falls in the middle of a CDS feature, it will be stretched as far as needed to completely include the CDS
- The extracted locus reverse-complemented if the targeted feature is in the reverse strand
- Output is also a GenBank file


## Usage:

```
usage: fetch_context_gbk.py [-h] -i INPUT -l LABEL [-e EXTRA] [-o OUTPUTFOLDER]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        GenBank file
  -l LABEL, --label LABEL
                        Target label (gene_id, locus_tag, etc.)
  -e EXTRA, --extra EXTRA
                        Number of kbps at either side of the target gene/protein/locus_tag to download. Default=20
  -o OUTPUTFOLDER, --outputfolder OUTPUTFOLDER
                        Where to put retrieved files. Default=output
```

## Requirements

* Python 3
* biopython
Provide a GenBank, a label (gene id, locus_tag, protein id) and extract another GenBank with a specified number of kilobases