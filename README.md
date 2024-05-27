# bed2seq


From a BED file, return the sequences according to the genome supplied

WARNING
This is an alpha version, it could contain issues.



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
                        Output file (default: <input_file>-bed2seq.tsv)
  -v, --version         show program's version number and exit
```