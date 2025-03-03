#!/bin/env python
"""
    extract_mpra_elements.py: fetch sequences of the given variants for MPRA experiments

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
import os
import os.path
import argparse
import gzip

def read_snpfile(fn, add_chr):
    if fn[-2:] == 'gz':
        f = gzip.open(fn, 'r')
    else:
        f = open(fn, 'r')

    positions = []

    for line in f:
        if line[0] == "#": 
            continue

        if line.find('\t') > 0:
            l = line.strip().split('\t')
        else:
            l = line.strip().split(' ')

        chrom = l[0]
        #add chr to the chromosome name of the variants
        if add_chr:
            chrom = "chr" + chrom

            if chrom == "chrMT":
                chrom = "chrM"

        snppos = int(l[1]) - 1 # 1-based converted to 0-based
        rsid = l[2]
        ref_allele = l[3]
        alt_allele = l[4]

        # handle multi-allelic variants
        if ',' in alt_allele:
            alt_allele = l[4].split(',')[0]
            sys.stderr.write("WARNING: %s is multi-allelic (%s).  The first allele (%s) will be used.\n" % (rsid, l[4], alt_allele))

        positions.append((chrom, snppos, rsid, ref_allele, alt_allele)) # chrom, snp, rsid, ref_allele, alt_allele

    f.close()
    return positions


def read_chromfile(filename):
    try:
        f = open(filename, 'r')
        seq = ''.join(map(lambda s: s.rstrip('\n').upper(), f.readlines()[1:]))
        f.close()

    except IOError as err:
        print("I/O error:", err , filename)
        sys.exit(0)

    return seq


def read_chromfile_gz(filename):
    try:
        f = gzip.open(filename, 'rt')
        seq = ''.join(map(lambda s: s.rstrip('\n').upper(), f.readlines()[1:]))
        f.close()

    except IOError as err:
        print("I/O error:", err , filename)
        sys.exit(0)
    return seq


def read_genomefile_chrom(genomefn, idxfn, chrom):
    chrom2lenoff = dict()
    try:
        fp = open(idxfn, 'r')
        for line in fp:
            f = line.strip().split('\t')
            chrom2lenoff[ f[0] ] = (int(f[1]), int(f[2]))
        fp.close()

        fp = open(genomefn, 'r')
        (chr_len, chr_offset) = chrom2lenoff[chrom]
        fp.seek(chr_offset)
        seqs = fp.read(chr_len)
        seq = ''.join(map(lambda s: s.strip().upper(), seqs.splitlines()))
        fp.close()

        return seq
    except IOError as err:
        print("I/O error:", err , filename)
        sys.exit(0)


def main(argv = sys.argv):
    desc_txt = "extract genomic DNA sequences from variant data for an MPRA experiment"

    parser = argparse.ArgumentParser(description=desc_txt, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("snpinfo_file", help="variant information file (VCF file format; only the first five columns are necessary)")
    parser.add_argument("element_length", type=int, help="desired (maximum) length of MPRA element design")
    parser.add_argument("-o", "--output", default='mpra_elements_output.txt', help="output file name")
    parser.add_argument("-d", "--genome-dir", default='.', help="set the directory of reference sequences. It searches for {CHROM}.fa[.gz] or genome.fa in the given dir")
    parser.add_argument("-m", "--max-indel", type=int, default=5, help="set the maximum length of indel allowed")
    parser.add_argument("-c", "--add-chr", action='store_true', help="if set, add chr in front of variant chromosome")

    args = parser.parse_args()

    snps = read_snpfile(args.snpinfo_file, args.add_chr)

    chrnames = sorted(set(map(lambda p: p[0], snps)))
    
    fout = open(args.output, 'w')
    for chrom in chrnames:
        sys.stderr.write("processing %s...\n" % (chrom))

        refseq_fn = os.path.join(args.genome_dir, '.'.join((chrom, 'fa')))
        refseq_gz_fn = os.path.join(args.genome_dir, '.'.join((chrom, 'fa', 'gz')))
        genome_fn = os.path.join(args.genome_dir, 'genome.fa')
        genome_idx_fn = os.path.join(args.genome_dir, 'genome.fa.fai')
        if os.path.isfile(refseq_fn):
            refseq = read_chromfile(refseq_fn)
        elif os.path.isfile(refseq_gz_fn): 
            refseq = read_chromfile_gz(refseq_gz_fn)
        elif os.path.isfile(genome_fn) and os.path.isfile(genome_idx_fn):
            refseq = read_genomefile_chrom(genome_fn, genome_idx_fn, chrom)
        else:
            sys.stderr.write("ERROR: cannot find chromosome/genome sequence file in the directory [%s].\n" % (args.genome_dir))
            sys.exit()

        for snp in sorted(filter(lambda p: p[0] == chrom, snps), key=lambda x: x[1]):
            ref_allele_len = len(snp[3])
            alt_allele_len = len(snp[4])

            # filter out INDEL variants that is longer than predefined threshold (default: 5bp)
            if ref_allele_len > (args.max_indel+1) or alt_allele_len > (args.max_indel+1):
                sys.stderr.write("WARNING: %s INDEL variant is longer than %d. skip.\n" % (snp[2], args.max_indel))
                continue

            # determine the extenion based on length of alleles and the maximum length of MPRA element
            # fixed length constructs for indels
            lextend_ref = int((args.element_length - ref_allele_len)/2)
            rextend_ref = args.element_length - ref_allele_len - lextend_ref
            lextend_alt = int((args.element_length - alt_allele_len)/2)
            rextend_alt = args.element_length - alt_allele_len - lextend_alt

            pos_start_ref = snp[1] - lextend_ref
            pos_end_ref = snp[1] + ref_allele_len + rextend_ref
            ref_allele_seq = refseq[pos_start_ref:pos_end_ref]
            if ref_allele_seq[lextend_ref:(lextend_ref+ref_allele_len)] != snp[3]:
                sys.stderr.write("WARNING: %s ref allele is not consistent. (%s =/= %s)\n" % (snp[2], ref_allele_seq[lextend_ref], snp[3]))

            pos_start_alt = snp[1] - lextend_alt
            pos_end_alt = snp[1] + ref_allele_len + rextend_alt
            alt_allele_seq = refseq[pos_start_alt:pos_end_alt]
            alt_allele_seq = alt_allele_seq[0:lextend_alt] + snp[4] + alt_allele_seq[(lextend_alt+ref_allele_len):]

            ref_allele_id = '_'.join([snp[0], str(snp[1]), snp[2], snp[3], snp[4], "ref"])
            alt_allele_id = '_'.join([snp[0], str(snp[1]), snp[2], snp[3], snp[4], "alt"])

            fout.write("%s\t%s\n" % (ref_allele_id, ref_allele_seq))
            fout.write("%s\t%s\n" % (alt_allele_id, alt_allele_seq))

    fout.close()

if __name__ == '__main__': main()
