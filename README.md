# bed2seq


From a BED file, return the sequences according to the genome supplied

## Installation 

### Solution 1 (Preferred)

In a virtualenv, install with pip

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
  -g, --genome genome   genome as fasta file
  -e, --extend EXTEND   extend the sequence ('-a 20' append 20 bp on each side)
  -r, --remove          only with '-e/--extend' option, keep only axtended part
  -n, --nostrand, --nostranded
                        don't reverse complement when strand is '-'
  -a, --add-columns ADD_COLUMNS [ADD_COLUMNS ...]
                        Add one or more columns to header (ex: '-a 3 AA' will 
                        add columns 3 and 27). The first column is '1' (or 'A')
  -d, --delimiter DELIMITER
                        with -a/--add-columns and a fasta format output, 
                        specifies a delimiter (default: '|')
  -o, --output OUTPUT   Output file (default: <input_file>-bed2seq.tsv)
  -v, --version         show program's version number and exit

```
