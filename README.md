# Genetic barcoding uncovers the clonal makeup of solid and liquid biopsies and their ability to capture intra tumoral heterogeneity

Code associated with the manuscript *Genetic barcoding uncovers the clonal makeup of solid and liquid biopsies and their ability to capture intra tumoral heterogeneity*  (2025, under review).

- Contact: Antonin Serrano, antonin.serrano@unimelb.edu.au
- Contact: Tom Weber, weber.ts@wehi.edu.au
- Contact: Delphine Merino, delphine.merino@onjcri.org.au

## Helpful links

- [bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2025.06.12.659267v1)
- Barcode-seq reads: Available from SRA at PRJNA1368999 once the paper has been accepted.

## generate barcode counts

To generate a barcode count matrix from a barcode-seq file using the edgeR processAmplicon function. 

> Note, fastq files must be uncompressed for this function.

```R
library(edgeR)

counts <-  processAmplicons(
  readfile = "path/to/read_R1.fastq",
  readfile2 = "path/to/read_R2.fastq",
  barcodefile = "index_tables/{run}_TrueIndexes.txt", # tsv of sample specific barcodes
  hairpinfile = "LongBarcodes.txt", # barcode sequences
  verbose = TRUE
)
```


