#!/bin/env python
"""
    select_mpra_barcodes.py: sample barcodes for MPRA experiments

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
import random
import operator
try:
    from itertools import imap
except ImportError:
    # Python 3...
    imap=map


def read_barcode_file(filename):
    fp = open(filename, 'r')
    barcodes = []
    for line in fp:
        barcodes.append(line.strip())
    fp.close()
    return barcodes

def hamming(str1, str2):
    ne = operator.ne
    return sum(imap(ne, str1, str2))

def main():
    desc_txt="select barcodes for each element that have at least pairwise N (default: 5) Hamming distances"

    parser = argparse.ArgumentParser(description=desc_txt, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("barcode_file", help="list of barcodes; the number of barcodes should be more than (num_elements * num_barcodes)")
    parser.add_argument("num_elements", type=int, help="number of enhancer elements to test")
    parser.add_argument("num_barcodes", type=int, help="number of barcodes per enhancer element")
    parser.add_argument("-o", "--output", default="selected_mpra_barcode_output.txt", help="output filename")
    parser.add_argument("-r", "--random-seed", type=int, default=1, help="set random seed for repproduciblity")
    parser.add_argument("-m", "--min-hamdist", type=int, default=5, help="set a minimum Hamming distance (i.e. number of mismatches) within barcode sets")
    parser.add_argument("-n", "--max-try", type=int, default=1000000, help="set a maximum number of tries")

    args = parser.parse_args()

    barcodes = read_barcode_file(args.barcode_file)
    nelem = args.num_elements
    nbarcode = args.num_barcodes

    random.seed(args.random_seed)

    elem_barcodes_matrix = []
    for i in range(args.num_elements):
        elem_barcodes = []
        for j in range(args.num_barcodes):
            failed = True
            attempts = 0
            while failed:
                if attempts > args.max_try:
                    sys.stderr.write("ERROR: You reached the maximum number of tries.  You should generate more barcodes first.\n")
                    exit()

                attempts += 1
                barcode_ind = random.randint(0, len(barcodes)-1)
                ibarcode = barcodes[barcode_ind]
                failed = False
                for ebarcode in elem_barcodes:
                    if hamming(ebarcode, ibarcode) < args.min_hamdist:
                        failed = True
                        break
            
                if not failed:
                    elem_barcodes.append(ibarcode)
                    barcodes.pop(barcode_ind)
	
        elem_barcodes_matrix.append(elem_barcodes)

    outfp = open(args.output, 'w')
    for i in range(len(elem_barcodes_matrix)):
        for j in range(len(elem_barcodes_matrix[i])):
            outfp.write("%d\t%d\t%s\n" % (i+1, j+1, elem_barcodes_matrix[i][j]))
    outfp.close()
    
if __name__=='__main__': main()
