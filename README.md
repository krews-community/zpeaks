# ZPeaks

A Dockerized application for finding peaks from alignment files.

## Steps

- Piles up alignments

## Requirements

- Memory: ZPeaks stores all pile-up data for all chromosomes given in dense arrays, in memory, at once. This means
you need
- CPUs: ZPeaks can parallelize for compute intensive tasks at a sub-chromosome level. It will run efficiently with 
any number of cores.