# nanopore_sequencing

The goal of this project is to generate the reference sequence needed to perform [Adaptive Sampling](https://nanoporetech.com/document/adaptive-sampling#sample-preparation-and-analysis-advanced) on the Oxford Nanopore.

Adaptive Sampling requires:
1. A Reference Genome
2. A BED file containing Regions of Interest (ROIs) + buffer sequences. 
These buffer sequences help ensure that reads that only partially cover the region of interest are still included.
