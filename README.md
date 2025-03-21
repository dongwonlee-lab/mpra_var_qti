
# MPRA for QT-interval associated variants

This is the repository for codes to replicate the MPRA experiment and analyses
published in the following manuscript:

Dongwon Lee, Lavanya Gunamalai, Jeerthi Kannan, Kyla Vickery, Or Yaacov,
Ana C. Onuchic-Whitford, Aravinda Chakravarti, and Ashish Kapoor.
Massively parallel reporter assays identify functional enhancer variants at QT interval GWAS loci.
*in prep*.

## Requirements

We have Python3 scripts for generating MPRA oligo sequences (MPRA Variant
Disigner), and Jupyter Notebooks for MPRAalyze analyses in R.

 * Python >= 3.7
 * R >= 4.x
 * R Packages
     * tidyverse
     * MPRAnalyze
     * ggplot2

## MPRA Variant Designer

There are four scripts in our MPRA Variant Designer tool developed in this
manuscript.

 1. `extract_mpra_elements.py`
 2. `gen_mpra_barcodes.py`
 3. `select_mpra_barcodes.py`
 4. `make_mpra_oligo.py`

These scripts are used to generate the MPRA sequence oligos for testing
variant effects on enhancer activity, and can be found in
`mpra_variant_designer` directory. We will explain how to use them
step-by-step in the next section.

### Extract Genomic sequences for test variants

```
Usage: extract_mpra_elements.py [options] <snpinfo_file> <element_length>

extract genomic DNA sequences from variant data for MPRA experiments

positional arguments:
  snpinfo_file          variant information file (VCF file format; only
                        the first five columns are necessary)
  element_length        desired (maximum) length of MPRA element design

options:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        output file name (default: mpra_elements_output.txt)
  -d GENOME_DIR, --genome-dir GENOME_DIR
                        set the directory of reference sequences. It
                        searches for {CHROM}.fa[.gz] or genome.fa in the
                        given dir (default: .)
  -m MAX_INDEL, --max-indel MAX_INDEL
                        set the maximum length of indel allowed (default: 5)
  -c, --add-chr         if set, add chr in front of variant chromosome
                        (default: False)
```

The first step is to extract genomic sequences for test variants by running
`extract_mpra_elements.py`. It also creates sequences with alternative alleles.

Note that this script also requires genome sequences in FASTA format. You can
use either `genome.fa` or FASTA files splitted by chromosomes, such as
`chr1.fa.[gz]`. These files are available from UCSC Genome Browser or NCBI.

* hg19: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/
* hg38: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/

Let's try to run this script to generate sequences for MPRA test elements using
the example file in the `data` directory. The genome assembly of the variants in
the example file is hg19. Let's also use 129 nt as the length of a test element,
so that we will end up with 200 bp of oligo length, which is the same as 
reported in the manuscript.


```
$ python mpra_variant_designer/extract_mpra_elements.py -d hg19/chromosomes data/test_variants.hg19.txt 129
```

This will generate `mpra_elements_output.txt`, which will then be used for
making the final oligo sequences later. The output file is tab-delimited,
and contains only two columns, `SEQ_ID` and `Sequence` without the header.
`SEQ_ID` is formatted as `[chrom]_[position]_[rsid]_[ref-allele]_[alt-allele]_[ref/alt]`.

##### Example Output: `mpra_elements_output.txt`

```
chr1_6147297_rs11583631_C_T_ref GGAAGGCATCTCCTCCGGACTCCGGTCACACTAGGCCCCCCACTGTCCAACTCCCACTGCCTGACTCACAGTTACTGCCTCGGTGGCTGCCTCCTGGCTCCATGAGGAGCCCACGCACCCCGGCTGGGC
chr1_6147297_rs11583631_C_T_alt GGAAGGCATCTCCTCCGGACTCCGGTCACACTAGGCCCCCCACTGTCCAACTCCCACTGCCTGATTCACAGTTACTGCCTCGGTGGCTGCCTCCTGGCTCCATGAGGAGCCCACGCACCCCGGCTGGGC
chr1_6147340_rs11584419_A_C_ref TGTCCAACTCCCACTGCCTGACTCACAGTTACTGCCTCGGTGGCTGCCTCCTGGCTCCATGAGGAGCCCACGCACCCCGGCTGGGCACGTGGTACACATCATATGATGGAGAAGTGCCCATTTCATGCC
chr1_6147340_rs11584419_A_C_alt TGTCCAACTCCCACTGCCTGACTCACAGTTACTGCCTCGGTGGCTGCCTCCTGGCTCCATGAGGCGCCCACGCACCCCGGCTGGGCACGTGGTACACATCATATGATGGAGAAGTGCCCATTTCATGCC
...
```

### Generate barcodes for MPRA elements

The second step is to create and assign barcodes to the MPRA test elements
generated above. This step is acheived by running `gen_mpra_barcodes.py` and
`select_mpra_barcodes.py`.

Note that the generated barcodes should meet the following criteria to ensure
that they introduce minimal technical biases:

* They contain all four nucleotides
* They have a GC content between 40 and 60 percent
* They do not contain restriction enzyme sites
* They do not contain homopolymers
* They do not contain miRNA binding sites

#### Generate all possible barcodes

```
Usage: gen_mpra_barcodes.py [options] <barcode_length>

generate all possible barcodes for MPRA experiments that meet the following
criteria: a) barcodes must contain all four nucleotides; b) GC contents must 
be between 40% and 60%; c) barcodes must NOT contain the predefined restriction
enzyme sites (default: SbfI, EcoRI) and homopolymer with predefined length
(default: 3bp) and longer; e) barcodes must NOT contain 7bp miRNA seeds in the
barcode with 6bp flanking sequences; f) barcodes with 1bp deletion + 1bp from
the 5' flanking side must not match any of the previously selected barcodes
(greedy search); and g) Hamming distance between barcodes must be at least 2
(greedy search)

positional arguments:
  barcode_length           desired barcode length

options:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        output filename (default: mpra_barcode_output.txt)
  -g MIN_GC, --min_gc MIN_GC
                        minimum GC contents allowed (default: 0.4)
  -G MAX_GC, --max_gc MAX_GC
                        maximum GC contents allowed (default: 0.6)
  -r RE_SEQ, --re-seq RE_SEQ
                        restriction enzyme sites to exclude SbfI and EcoRI
                        RE sites are defaults values (default: CCTGCAGG,GAATTC)
  -l LENGTH_HOMOPOLYMER, --length-homopolymer LENGTH_HOMOPOLYMER
                        length of homopolymer to exclude (default: 3)
  -m MIRNA_FILE, --mirna-file MIRNA_FILE
                        miR database file from targetscan.org
                        (default: miR_Family_Info.txt)
  -u UPSTREAM_SEQ, --upstream-seq UPSTREAM_SEQ
                        5' upstream flanking sequence of barcodes for
                        scanning of bad sequences (RE & miR) (default: GAATTC)
  -d DOWNSTREAM_SEQ, --downstream-seq DOWNSTREAM_SEQ
                        3' downstream flanking sequence of barcodes for
                        scanning of bad sequences (RE & miR) (default: CATTGC)
  -s SPECIES_IDS, --species_ids SPECIES_IDS
                        species IDs of miRNA for scanning (default: 9606,10090)
```

`gen_mpra_barcodes.py` generates all possible barcodes with a specific length
that meet the several criteria as described above. In order to exclude the
barcodes matching miRNA binding site, this script requires miRNA seeds as an
input, which can be found here:
https://www.targetscan.org/vert_80/vert_80_data_download/miR_Family_Info.txt.zip

After downloading this, try the following command to generate barcodes with
length of 13 nt. It will generate `mpra_barcode_output.txt` containing barcodes.

```
$ python mpra_variant_designer/gen_mpra_barcodes.py -m miR_Family_Info.txt 13
```

This will take about ~10 mins, as all possible k-mers are evaluated.

##### Example Output: `mpra_barcode_output.txt`

```
AAGAACAACCGTC
AAGAACAACGTCC
AAGAACAACTCCG
AAGAACAATCGCC
AAGAACACAGTCG
...
```

#### Assign barcodes to MPRA test elements

```
Usage: select_mpra_barcodes.py [options] <barcode_file> <num_elements> <num_barcodes>

select barcodes for each element that have at least pairwise N (default: 5)
Hamming distances

positional arguments:
  barcode_file          list of barcodes; the number of barcodes should be
                        more than (num_elements * num_barcodes)
  num_elements          number of enhancer elements to test
  num_barcodes          number of barcodes per enhancer element

options:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        output filename
                        (default: selected_mpra_barcode_output.txt)
  -r RANDOM_SEED, --random-seed RANDOM_SEED
                        set random seed for repproduciblity (default: 1)
  -m MIN_HAMDIST, --min-hamdist MIN_HAMDIST
                        set a minimum Hamming distance (i.e. number of
                        mismatches) within barcode sets (default: 5)
  -n MAX_TRY, --max-try MAX_TRY
                        set a maximum number of tries (default: 1000000)
```

In MPRA, multiple barcodes should be assigned to each test element for robust
estimation of enhancer activit. This barcode assignment process is done by
running `select_mpra_barcodes.py`. This script takes three arguments: (1) full
list of barcodes generated from the previous step, (2) number of test MPRA
elements (sequences with aleteative alleles are consider as distinct MPRA test
elements), and (3) the number of barcodes to assign per element. Try this script
with the following command. It will generate `selected_mpra_barcode_output.txt`
containing three columns: element index, barcode index in a element, and
selected barcodes.

```
$ python mpra_variant_designer/select_mpra_barcodes.py mpra_barcode_output.txt 20 50 
```


##### Example Output: `selected_mpra_barcode_output.txt`

```
1       1       ACGTGGACTCGTA
1       2       GCCGGTATATCTA
1       3       TAACACTTCGCCG
1       4       GTGATGCGTCGTA
1       5       GTCATCACATCCG
...
```

### Generate MPRA oligo construct

```
Usage: make_mpra_oligo.py [options] cre_file barcode_file

Generate a full set of oligo sequences for an MPRA experiment. Sequence oligos
will be built as follows: ForwardPrimer-CRE-RE1-spacer-RE2-Barcode-ReversePrimer

options:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        output file name (default: mpra_oligo_output.txt)
  -f FORWARD_PRIMER, --forward-primer FORWARD_PRIMER
                        Forward primer sequence (default: CACTGAGGACCGGATCAACT)
  -r REVERSE_PRIMER, --reverse-primer REVERSE_PRIMER
                        Reverse primer sequence (default: CATTGCGTGAACCGAGACCT)
  -e RE1, --re1 RE1     restriction enzyme site 1 (default: CCTGCAGG)
  -E RE2, --re2 RE2     restriction enzyme site 2 (default: GAATTC)
  -s SPACER, --spacer SPACER
                        sequence for spacer between re1 and re2 (default: AAAA)
```

Finally, the full set of sequences for oligo synthesis can be made by running
`make_mpra_oligo.py`. It takes two arguments, the MPRA test elements and the
selected barcodes, both of which were generated above. Note that all the given
sequences, including test elements, barcodes, primers, and restriction enzyme
site sequences, will be concatenated as is. Try the following is an example
command to make the final output, `mpra_oligo_output.txt`. The output contains
two columns, `OLIGO_ID` and `OLIGO_SEQ`. `OLIGO_ID` is formatted as
`[SEQ_ID]_[BARCODE]`.

```
$ python make_mpra_oligo.py mpra_elements_output.txt selected_mpra_barcode_output.txt
```

##### Example Output: `mpra_oligo_output.txt`

```
chr1_6147297_rs11583631_C_T_ref_ACGTGGACTCGTA   CACTGAGGACCGGATCAACT-[sequence]-CCTGCAGG-AAAA-GAATTC-ACGTGGACTCGTA-CATTGCGTGAACCGAGACCT
chr1_6147297_rs11583631_C_T_ref_GCCGGTATATCTA   CACTGAGGACCGGATCAACT-[sequence]-CCTGCAGG-AAAA-GAATTC-GCCGGTATATCTA-CATTGCGTGAACCGAGACCT
...
```

Your MPRA sequences for oligo synthesis are now ready!


** Please email Dongwon Lee (dongwon.lee AT childrens DOT harvard DOT edu) if you have any questions. **
