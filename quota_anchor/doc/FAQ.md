# FAQ

## GFF format

1. The software uses the regex: `r'^(\S+)\t(\S+)\tCDS\t(\d+)\t+(\d+)\t+(\S+)\t+(\S+)\t+(\S+)\t+.*Parent=([\s\S]+?);'` and `r'^(\S+)\t(\S+)\tCDS\t(\d+)\t+(\d+)\t+(\S+)\t+(\S+)\t+(\S+)\t+.*Parent=([\S\s]+?)$'` to parse the CDS feature lines in the gff file.
2. The software uses the regex: `r'^(\S+)\t(\S+)\t\S+\t(\d+)\t+(\d+)\t+(\S+)\t+(\S+)\t+(\S+)\t+.*ID=([\s\S]+?);Parent=([\s\S]+?);'` and `r'^(\S+)\t(\S+)\t\S+\t(\d+)\t+(\d+)\t+(\S+)\t+(\S+)\t+(\S+)\t+.*ID=([\s\S]+?);Parent=([\s\S]+?)$'` to parse the mRNA feature lines in the gff file.
