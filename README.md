# PREMISE: A probabilistic framework for source assignment of Illumina reads

This repository contains an Expectation-Maximization (EM)-based metagenomic classifier.

## Installation

To install you must first have cargo and rustup installed:
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

After installing the above command you can run the following to install premise:
```bash
cargo install --git=https://github.com/sriram98v/premise
```

Alternatively, you can install premise by cloning this repository and building it locally:
```bash
git clone https://github.com/sriram98v/premise
cd premise
cargo install --path=./
```

## Usage
This tool aligns reads to references by first building a suffix array, then aligning the reads to the references by finding seed matches and extending them. After finding all alignments, an EM algorithm is employed to identify the most-likely true source for each read, and to estimate population proportions of the detected references.

### Building an Index
In order to start classifying a set of reads, you need to build an index from you set of reference sequences. You can build an index using the following command.
```bash
premise build -s <reference sequence file path>
```
This will create an index file with the extension ```.fmidx```.

### Classifying reads
You can now begin classifying read using the following command.
```bash
premise query -s <reference sequence file path> -p <Percent mismatch> -r <Reads> -t <Threads [Default: 2]> -o <Output file [Default: out.matches]>
```

You can refer the man page for premise for more details by running
```bash
premise -h
```
