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
import textwrap
import ascii
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

    resp = compute(args, chr_dict)
    # ~ print(*res, sep='\n---\n')
    write(args, resp)


def ctrl_args(args):
    if args.remove and not args.extend:
        sys.exit(f"{COL.RED}option '--remove' needs '--append' option.")


def _tab_length(rows):
    """ Function doc """
    tab_len = 0
    for row in rows:
        if not row.startswith('#'):
            tab_len = len(row.split('\t'))
            break
    return tab_len


def _input_ok(args, rows, resp, chr_dict, cols_id):
    ### find first no commented line and check it
    for i,row in enumerate(rows):
        if row.startswith('#'):
            continue
        fields = row.rstrip('\n').split('\t')
        nfields = len(fields)
        try:
            chr, start, end = fields[:3]
        except ValueError:
            resp["error"] = f"not enough columns at line {i+1} (check your bed file)."
            resp["is_ok"] = False
            return False

        if nfields < 6 and not args.nostrand:
            resp["warning"].append("strand column missing: strands cannot be evaluated.")

        ### Try to accoding chr column
        if chr not in chr_dict:
            ### chromosome begin by 'chr' in user file, not in genome --> need to remove chr
            if chr.lstrip('chr') in chr_dict:
                args.chr = '-chr'
            ### genome begin by 'chr', but not in user file --> need to add chr
            if f"chr{chr}" in chr_dict:
                args.chr = '+chr'

        ### check if --add-columns is compatible with number of columns
        if args.add_columns and max(cols_id) > nfields:
            resp["error"] = (f"BED file has {nfields} columns, but you asked for "
                  f"{max(cols_id)}.")
            resp["is_ok"] = False
            return False

        '''
        ### Check some commonly issues
        if chr not in chr_dict:
            resp["error"] = ("chromosomes are not named in the same way in the "
                      "query and the genome file. Below the first chromosome found: \n"
                     f"     your query: {chr}\n"
                     f"     genome: {next(iter(chr_dict.keys()))}\n"
                     f"   Please, correct your request.")
            resp["is_ok"] = False
            return False
        '''
        break
    return True



def compute(args, chr_dict):

    resp = {
        "is_ok": True,
        "result": [],
        "warning": [],
        "error": None
    }

    ### convert input as row list
    if isinstance(args.input, str):
        rows = args.input.splitlines()
    else:
        rows = args.input.read().splitlines()

    ### How many columns in the BED file
    tab_len = _tab_length(rows)

    ### columns characters are converted as index, ex: AA -> 27
    cols_id = ascii.get_index(args.add_columns)

    ### check input syntax
    if not _input_ok(args, rows, resp, chr_dict, cols_id):
        return resp

    for i,row in enumerate(rows):
        fields = row.rstrip('\n').split('\t')
        # ~ nfields = len(fields)
        if tab_len >= 6:
            chr, start, end, id_seq, score, strand = fields[:6]
        elif tab_len >= 4:
            chr, start, end, id_seq = fields[:4]
        else:
            chr, start, end = fields[:3]

        if tab_len < 4:
            id_seq = f"sequence_{i+1}"

        start = int(start) - args.extend
        end = int(end) +  args.extend

        ### handle chromosome
        if hasattr(args, 'chr'):
            if args.chr == '-chr':
                chr = chr.lstrip('chr')
            else:
                chr = f"chr{chr}"

        ### get sequence
        seq = chr_dict[chr][start:end]

        ### Handle strand
        seq = seq.complement.reverse.seq if tab_len >=6 and strand == '-' and not args.nostrand else seq.seq

        ### Handle remove
        if args.remove:
            seq = seq[:args.extend] + seq[-args.extend:]

        ### append additional selected columns to the header
        if cols_id:
            added_cols = f"{args.delimiter}{args.delimiter.join([fields[num-1] for num in cols_id])}"
            id_seq += added_cols

        ### push in results
        resp["result"].append(f">{id_seq}\n{textwrap.fill(seq, width=100)}")

    return resp


def write(args, resp):

    ### define output file
    if not args.output:
        name, ext = os.path.splitext(os.path.basename(args.input.name))
        args.output = f"{name}-bed2seq.fa"

    if resp["is_ok"]:
        ## write results in file
        if resp["result"]:
            with open(args.output, 'w') as fh:
                for result in resp["result"]:
                    fh.write(f"{result}\n")
            print(f"\n🧬 {args.output} succefully created.\n")
        ### WARNINGS
        if resp["warning"]:
            print(f"{COL.PURPLE}⚠️  Warnings:")
            for warning in resp["warning"]:
                print(f"   - {warning}")
            print(COL.END)
    else:
        print(f"\n☠️  {COL.RED}Error: {resp['error']}\n")


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
    parser.add_argument("input",
                        help="bed file",
                        type=argparse.FileType('r'),
                       )
    parser.add_argument("-g", "--genome",
                        help="genome as fasta file",
                        metavar="genome",
                        required=True,
                       )
    parser.add_argument('-e', '--extend',
                        type=int,
                        help="extend the sequence ('-a 20' append 20 bp on each side)",
                        default=0,
                       )
    parser.add_argument('-r', '--remove',
                        action="store_true",
                        help="only with '-e/--extend' option, keep only axtended part",
                        )
    parser.add_argument('-n', '--nostrand', '--nostranded',
                        action="store_true",
                        help="don't reverse complement when strand is '-'",
                        )
    parser.add_argument("-a", "--add-columns",
                        help="Add one or more columns to header (ex: '-a 3 AA' will add columns "
                             "3 and 27). The first column is '1' (or 'A')",
                        nargs= '+',
                        )
    parser.add_argument("-d", "--delimiter",
                        help="with -a/--add-columns and a fasta format output, specifies a delimiter (default: space)",
                        default= ' ',
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
