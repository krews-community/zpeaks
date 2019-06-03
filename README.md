# ZPeaks

A Dockerized application for finding peaks from alignment files.

## Steps

![ZPeaks Flow](img/zpeaks_flow.png)

## Requirements

- Memory: ZPeaks stores all pile-up data for all chromosomes given in dense arrays, in memory, at once. As a rule 
of thumb, you should use at least 5x the number of base pairs in the genome + 1GB. For human, this means 16GB.
- CPUs: ZPeaks can parallelize for compute intensive tasks at a sub-chromosome level. It will run efficiently with 
any number of cores.

## Running

`java -jar `