#!/usr/bin/env python3

"""
From a BED file, return the sequences according to the genome supplied

WARNING
This is an alpha version, it could contain issues.
contact me in case of bug at benoit.guibert@inserm.fr
"""

import sys
import os
import argparse
try:
    import pyfaidx
except ModuleNotFoundError as err:
    sys.exit(err)

import info


def main():
    """ Function doc """
    args = usage()
    ctrl_args(args)
    try:
        chr_dict = pyfaidx.Fasta(args.genome) # if fai file doesn't exists, it will be automatically created
    except pyfaidx.FastaNotFoundError as err:
        sys.exit(f"FastaNotFoundError: {err}")
    except OSError as err:
        sys.exit(f"\n{COL.RED}WriteError: directory {os.path.dirname(args.genome)!r} may not be "
                  "writable.\nIf you can't change the rights, you can create a symlink and target "
                  f"it. For example:\n  ln -s {args.genome} $HOME\n{COL.END}")
    bed_ctrl(args, chr_dict)
    res = compute(args, chr_dict)
    # ~ print(*res, sep='\n---\n')
    write(args, res)


def ctrl_args(args):
    if args.remove and not args.append:
        sys.exit(f"{COL.RED}option '--remove' needs '--append' option.")


def bed_ctrl(args, chr_dict):
    with open(args.bed.name) as fh:
        for row in fh:
            if row.startswith('#'):
                continue
            chr, start, end, name, score, strand, *rest = row.rstrip('\n').split('\t')

            ### Check some commonly issues
            if chr not in chr_dict:
                sys.exit(f"{COL.RED}ErrorChr: Chromosomes are not named in the same way in the "
                          "query and the genome file. Below the first chromosome found: \n"
                         f" your query: {chr}\n"
                         f" genome: {next(iter(chr_dict.keys()))}\n"
                         f"Please, correct your request (or modify the file '{args.genome}.fai').")
            break


def compute(args, chr_dict):

    res = []

    for row in args.bed:
        chr, start, end, name, score, strand, *ext = row.rstrip().split('\t')
        start = int(start) - args.append
        end = int(end) +  args.append

        seq = chr_dict[chr][start:end]

        ### Handle strand
        seq = seq.complement.reverse.seq if strand == '-' and not args.nostrand else seq.seq

        ### Handle remove
        if args.remove:
            seq = seq[:args.append] + seq[-args.append:]

        ### push in results
        res.append(f">{name}\n{seq}\n")

    return res


def write(args, res):
    """ Function doc """
    ### define output file
    if not args.output:
        name, ext = os.path.splitext(os.path.basename(args.bed.name))
        args.output = f"{name}-bed2seq.fa"
    ### write results in file
    with open(args.output, 'w') as fh:
        if not res:
            return
        for row in res:
            fh.write(row)



class COL:
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'


def usage():
    doc_sep = '=' * min(80, os.get_terminal_size(2)[0])
    parser = argparse.ArgumentParser(description= f'{doc_sep}{__doc__}{doc_sep}',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,)
    parser.add_argument("bed",
                        help="bed file",
                        type=argparse.FileType('r'),
                       )
    parser.add_argument("-g", "--genome",
                        help="genome as fasta file",
                        metavar="genome",
                        required=True,
                       )
    parser.add_argument('-a', '--append',
                        type=int,
                        help="enlarge the sequence ('-a 20' append 20 bp on each side)",
                        default=0,
                       )
    parser.add_argument('-r', '--remove',
                        action="store_true",
                        help="only with '--append' option, keep only appended part",
                        )
    parser.add_argument('-n', '--nostrand',
                        action="store_true",
                        help="don't reverse complement when strand is '-'",
                        )
    parser.add_argument("-o", "--output",
                        type=str,
                        help=f"Output file (default: <input_file>-{info.APPNAME}.tsv)",
                        )
    parser.add_argument('-v', '--version',
                        action='version',
                        version=f"{parser.prog} v{info.VERSION}",
                       )
    ### Go to "usage()" without arguments or stdin
    if len(sys.argv) == 1 and sys.stdin.isatty():
        parser.print_help()
        sys.exit()
    return parser.parse_args()


if __name__ == "__main__":
    main()

import sys
import argparse