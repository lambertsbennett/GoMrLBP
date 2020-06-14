# Multi-resolution local binary patterns in Go

This package allows you to compute local binary patterns from sequence data.
Inspired by [Gouchaki et al. (2019)](https://www.nature.com/articles/s41598-018-38197-9).

GoMrLBP carries out several functions to transform assembled contigs into local binary pattern histograms:
- Integer representation of sequences.
- Calculation of LBP codes.
- Construction of LBP histogram.
- Singular value decomposition of LBP histograms.

Usage example:
```
gomrlbp -n NUM_PROCESSORS -file CONTIG_FILE -o OUTPUT_FILE
```
This implementation is an improvement over the C++ implementation in that all operations are carried out in parallel, allowing entire genomes/transcriptomes to be processed in approx. 1 minute. 

INPUTS: Contig fasta file (Can be gzipped).
OUTPUTS: CSV file with:
- Sequence header (read ID)
- Species ID (if known)
- LBP histogram
- SVD results.
