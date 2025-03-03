#!/bin/env python

"""
    gen_mpra_barcodes.py: generate barcodes for MPRA experiments

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
import re
import argparse

#Trie structure for efficient finding of existing barcodes within certain Hamming distance
class Node:
    def __init__(self, letter=None):
        self.letter = letter
        self.children = dict()

    def addChild(self, key):
        self.children[key] = Node(key)

    def __getitem__(self, key):
        return self.children[key]

class Trie:
    def __init__(self):
        self.head = Node()

    def __getitem__(self, key):
        return self.head.children[key]

    def add(self, word):
        current_node = self.head
        word_finished = True

        for i in range(len(word)):
            if word[i] in current_node.children:
                current_node = current_node.children[word[i]]
            else:
                word_finished = False
                break

        if not word_finished:
            while i < len(word):
                current_node.addChild(word[i])
                current_node = current_node.children[word[i]]
                i += 1

    def DFS_recursive(self, cnode, word, level, nmismatches, hamming_dist):
        #print word, level, nmismatches
        if level == len(word):
            return True

        for letter in cnode.children:
            if letter == word[level]:
                if self.DFS_recursive(cnode.children[letter], word, level+1, nmismatches, hamming_dist):
                    return True
            else:
                if nmismatches == hamming_dist:
                    continue
                else:
                    if self.DFS_recursive(cnode.children[letter], word, level+1, nmismatches+1, hamming_dist):
                        return True
        return False


    def has_word_within_hamming_dist(self, word, hamming_dist=0):
        return self.DFS_recursive(self.head, word, 0, 0, hamming_dist)


def revcomp(seq):
    rc = {'A':'T', 'G':'C', 'C':'G', 'T':'A'}
    return ''.join([rc[seq[i]] for i in range(len(seq)-1, -1, -1)])


def id2kmer(kmerid, k):
    kmer = ''
    nts = ['A', 'C', 'G', 'T']
    for i in range(k):
        kmer = nts[(kmerid % 4)] + kmer
        kmerid = int(kmerid/4)

    return kmer


def read_mirna_file(filename, species_ids):
    mirna_seeds = set()
    fp = open(filename, 'rt')
    header = fp.readline() #remove header
    for line in fp:
        f = line.strip().split()
        if int(f[2]) in species_ids:
            #mirna_seeds.add(f[1].replace("U", "T"))
            seed_rc = revcomp(f[1].replace("U", "T"))
            mirna_seeds.add(seed_rc)
    fp.close()

    sys.stderr.write("%d distinct miR seeds loaded.\n" % len(mirna_seeds)) 
    return mirna_seeds

def main():
    desc_txt="generate all possible barcodes for MPRA experiments that meet the following criteria:\
            a) barcodes must contain all four nucleotides;\
            b) GC contents must be between 40% and 60%;\
            c) barcodes must NOT contain the predefined restriction enzyme sites (default: SbfI, EcoRI) and \
            homopolymer with predefined length (default: 3bp) and longer;\
            e) barcodes must NOT contain 7bp miRNA seeds in the barcode with 6bp flanking sequences;\
            f) barcodes with 1bp deletion + 1bp from the 5' flanking side\
               must not match any of the previously selected barcodes (greedy search);\
            g) Hamming distance between barcodes must be at least 2 (greedy search)"

    parser = argparse.ArgumentParser(description=desc_txt, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("barcode_length", type=int, help="desired barcode length")
    parser.add_argument("-o", "--output", default="mpra_barcode_output.txt", help="output filename")
    parser.add_argument("-g", "--min_gc", type=float, default=0.4, help="minimum GC contents allowed") 
    parser.add_argument("-G", "--max_gc", type=float, default=0.6, help="maximum GC contents allowed") 
    parser.add_argument("-r", "--re-seq", default="CCTGCAGG,GAATTC", help="restriction enzyme sites to exclude. SbfI and EcoRI RE sites are defaults values")
    parser.add_argument("-l", "--length-homopolymer", type=int, default=3, help="length of homopolymer to exclude")
    parser.add_argument("-m", "--mirna-file", default="miR_Family_Info.txt", help="miR database file from targetscan.org") 
    parser.add_argument("-u", "--upstream-seq", default="GAATTC", help="5' upstream flanking sequence of barcodes for scanning of bad sequences (RE & miR)")
    parser.add_argument("-d", "--downstream-seq", default="CATTGC", help="3' downstream flanking sequence of barcodes for scanning of bad sequences (RE & miR)")
    parser.add_argument("-s", "--species_ids", default="9606,10090", help="species IDs of miRNA for scanning") 

    args = parser.parse_args()

    #exclude_str = "AAA|CCC|GGG|TTT|GGTACC|TCTAGA"
    exclude_str = "A{%d}" % (args.length_homopolymer)
    exclude_str += "|C{%d}" % (args.length_homopolymer)
    exclude_str += "|G{%d}" % (args.length_homopolymer)
    exclude_str += "|T{%d}" % (args.length_homopolymer)
    exclude_str += ('|'+'|'.join(args.re_seq.split(',')))

    flanking5p = args.upstream_seq
    flanking3p = args.downstream_seq
    species_ids = list(map(int, args.species_ids.split(',')))
    gcc_min = args.min_gc
    gcc_max = args.max_gc

    kmerlen = args.barcode_length

    mirna_seeds = read_mirna_file(args.mirna_file, species_ids)
    mirna_seeds_trie = Trie()
    for seed in mirna_seeds:
        mirna_seeds_trie.add(seed)

    patternA1 = re.compile("A")
    patternC1 = re.compile("C")
    patternG1 = re.compile("G")
    patternT1 = re.compile("T")
    pattern_include = (patternA1, patternC1, patternG1, patternT1)
    patterns_exclude = re.compile(exclude_str)

    barcodes = set()
    barcodes_trie = Trie()
    found = 0

    outfp = open(args.output, 'w')

    for kid in range(4**kmerlen):
        if kid % 10000 == 0:
            sys.stderr.write("\r%d found.. %d/%d (%.2f%%) examined.." % (found, kid, 4**kmerlen, kid*100/float(4**kmerlen)))

        kmer = id2kmer(kid, kmerlen)
        failed = False

        #1. must contain all four nucleotides
        for p in pattern_include:
            if not p.search(kmer):
                failed = True
                break
        if failed:
            continue

        #2. GC contents must be between 40% and 60%"
        gcc = (kmer.count("C") + kmer.count("G"))/float(len(kmer))
        if gcc<gcc_min or gcc>gcc_max:
            continue

        #3. must NOT contain the predefined patterns
        if patterns_exclude.search(kmer):
            continue

        #4. must NOT contain 7bp miRNA seeds in the barcode with 6bp flanking sequences
        kmer_ext = flanking5p + kmer + flanking3p
        for i in range(len(kmer_ext)-7+1):
            if mirna_seeds_trie.has_word_within_hamming_dist(kmer_ext[i:(i+7)], 0):
                failed = True
                break
        if failed:
            continue

        #5. barcodes with single base deletion + one base from the 5-prime flanking side must not match existing barcodes (greedy search)"
        for i in range(len(kmer)):
            kmer1del_f5p = flanking5p[-1:] + kmer[:i] + kmer[(i+1):]
            if kmer1del_f5p in barcodes:
                failed = True
                break
        if failed:
            continue

        #6. Hamming distance with existing barcodes should be at least 2
        if barcodes_trie.has_word_within_hamming_dist(kmer, 1):
            continue

        barcodes.add(kmer)
        barcodes_trie.add(kmer)

        found += 1
        outfp.write(kmer + "\n")

    sys.stderr.write("\r%d found.. %d/%d (%.2f%%) examined.. done!\n" % (found, kid+1, 4**kmerlen, (kid+1)*100/float(4**kmerlen)))
    outfp.close()

if __name__=='__main__': main()
