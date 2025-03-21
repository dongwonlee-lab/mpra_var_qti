#!/bin/env python
"""
    make_mpra_oligos.py: make a full set of oligo sequences for an MPRA experiment

    Copyright (C) 2025 Dongwon Lee

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import sys
import argparse

def read_cre_file(filename):
    creids = []
    creseqs = []
    fp = open(filename, 'r')
    for line in fp:
        f = line.strip().split('\t')
        creids.append(f[0])
        creseqs.append(f[1])
    fp.close()
    return creids, creseqs

def read_barcode_file(filename):
    barcodes = []
    barcode_set = []
    fp = open(filename, 'r')
    for line in fp:
        f = line.strip().split('\t')
        creno = int(f[0])
        barcodeno = int(f[1])

        # add a new barcode set
        if barcodeno == 1:
            barcodes.append([])

        barcodes[creno-1].append(f[2])
    fp.close()

    return barcodes

def main():
    desc_txt = "Generate a full set of oligo sequences for an MPRA experiment. Sequence oligos will be built as follows: ForwardPrimer-CRE-RE1-spacer-RE2-Barcode-ReversePrimer"

    parser = argparse.ArgumentParser(description=desc_txt, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("cre_file", help="list of cis-regulatory elements to test") 
    parser.add_argument("barcode_file", help="list of barcodes generated by select_mpra_barcodes.py") 
    parser.add_argument("-o", "--output", default='mpra_oligo_output.txt', help="output file name")
    parser.add_argument("-f", "--forward-primer", type=str, default="CACTGAGGACCGGATCAACT", help="Forward primer sequence")
    parser.add_argument("-r", "--reverse-primer", type=str, default="CATTGCGTGAACCGAGACCT", help="Reverse primer sequence")
    parser.add_argument("-e", "--re1", type=str, default="CCTGCAGG", help="restriction enzyme site 1")
    parser.add_argument("-E", "--re2", type=str, default="GAATTC", help="restriction enzyme site 2")
    parser.add_argument("-s", "--spacer", type=str, default="AAAA", help="sequence for spacer between re1 and re2")

    args = parser.parse_args()

    creids, creseqs = read_cre_file(args.cre_file)
    barcodes = read_barcode_file(args.barcode_file)
    
    if len(creids) < len(barcodes):
        sys.stderr.write("WARNING: the number of CREs is smaller than the number of barcode sets. Remaining barcodes will be created as empty vector.\n")
        for i in range(len(barcodes) - len(creids)):
            creids.append("Empty" + str(i+1))
            creseqs.append("")

    if len(creids) > len(barcodes):
        sys.stderr.write("WARNING: the number of CREs is larger than the number of barcode sets. Exit.\n")
        sys.exit()

    fout = open(args.output, 'w')
    for creid, creseq, barcodes in zip(creids, creseqs, barcodes):
        for barcode in barcodes:
            oligo_id = '_'.join([creid, barcode])
            oligo = args.forward_primer + creseq + \
                    args.re1 + args.spacer + args.re2 + \
                    barcode + args.reverse_primer

            fout.write("%s\t%s\n" % (oligo_id, oligo))
    fout.close()

if __name__=='__main__': main()
