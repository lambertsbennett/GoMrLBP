 [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Multi-resolution local binary patterns in Go

This package allows you to compute local binary pattern histograms from sequence data.
Inspired by [Kouchaki et al. (2019)](https://www.nature.com/articles/s41598-018-38197-9).

GoMrLBP carries out several functions to transform assembled contigs into local binary pattern histograms:
- Integer representation of sequences.
- Calculation of LBP codes.
- Construction of LBP histogram.
- Singular value decomposition of LBP histograms.

Usage example:
```
gomrlbp -n NUM_PROCESSORS -file CONTIG_FILE -o OUTPUT_FILE -max-win MAX WINDOW LBP -single USE SINGLE WINDOW
```
This implementation is an improvement over the C++ implementation in that all operations (except SVD) are carried out in parallel, allowing genomes/transcriptomes to be processed rapidly. 

----

**INPUTS:** Contig fasta file (Can be gzipped).

**OUTPUTS:** CSV file with:
- Sequence header (read ID)
- Species ID (if known)
- LBP histogram
- SVD results.

----

**Still under development.**

To add:
- Option to output as parquet file.
- Consistent error handling.
- Unit tests.

----

This package makes extensive use of gonum (https://github.com/gonum), gorgonia (https://github.com/gorgonia/gorgonia), and the truncated SVD implementation from James Bowman's nlp package (https://github.com/james-bowman/nlp).
