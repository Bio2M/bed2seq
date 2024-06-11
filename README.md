# bed2seq


From a BED file, return the sequences according to the genome supplied

## Installation 

### Solution 1 (Preferred)

Install with pip

```
pip install bed2seq
```

### Solution 2

Installation from github:

```
git clone https://github.com/Bio2M/bed2seq.git
```

The `pyfaidx` python package is required, install it with `pip`, `apt` or  your preferred method.

## usage

```
positional arguments:
  bed                   bed file

options:
  -h, --help            show this help message and exit
  -g genome, --genome genome
                        genome as fasta file
  -a APPEND, --append APPEND
                        enlarge the sequence ('-a 20' append 20 bp on each side)
  -r, --remove          only with '--append' option, keep only appended part
  -n, --nostrand        don't reverse complement when strand is '-'
  -o OUTPUT, --output OUTPUT
                        Output file (default: <input_file>-bed2seq.fa)
  -v, --version         show program's version number and exit
```