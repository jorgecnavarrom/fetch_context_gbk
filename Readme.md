# Fetch genomic context from a GenBank file

## What it does

This script reads a GenBank file and, given a label, it will try to find a matching qualifier. If found, it will extract it, plus some extra **bps** upstream/downstream (default: 20 kbps)


## Notes

- The label must be complete (no fuzzy search)
- The script only looks at CDS features
- Qualifiers where the label will be searched for: `gene`, `protein_id`, `proteinId`, `locus_tag` and `name`
- If the extension falls in the middle of a CDS feature, it will be stretched as far as needed to completely include the CDS
- The extracted locus is reverse-complemented if the targeted feature is in the reverse strand
- Output is also a GenBank file
- For even more control at each side, use the `--downstream` and `--upstream` paramteres


## Usage:

```
usage: fetch_context_gbk.py [-h] -i INPUT -l LABEL [-e EXTRA] [-u UPSTREAM] [-d DOWNSTREAM] [-o OUTPUTFOLDER]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        GenBank file
  -l LABEL, --label LABEL
                        Target label (gene_id, locus_tag, etc.)
  -e EXTRA, --extra EXTRA
                        Number of bps at either side of the target gene/protein/locus_tag to extract. Default=20000
  -u UPSTREAM, --upstream UPSTREAM
                        Override extension upstream of target with desired bps
  -d DOWNSTREAM, --downstream DOWNSTREAM
                        Override extension downstream of target by desired bps
  -o OUTPUTFOLDER, --outputfolder OUTPUTFOLDER
                        Where to put retrieved files. Default=output
```

## Requirements

* Python 3
* biopython
