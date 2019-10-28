#! /bin/bash

echo "Simulate HiC"
# python Simulate_HiC.py

echo ""
echo ""
echo ""
echo ""
echo "BAM2count"
python ../scripts/waffle-bam2count.py --bam data/fake.bam -r 10000 -C 8 -o tmp -b data/biases.pickle

echo ""
echo ""
echo ""
echo ""
echo "metaplot intra A vs B"
python ../scripts/waffle-peaks.py --genomic_matrix tmp/tmp_bam_10kb.tsv --peaks data/peaks_protA.bed data/peaks_protB.bed -o tmp/tmp.pickle -s 5 -w intra

echo ""
echo ""
echo ""
echo ""
echo "metaplot inter A vs B"
python ../scripts/waffle-peaks.py --genomic_matrix tmp/tmp_bam_10kb.tsv --peaks data/peaks_protA.bed data/peaks_protB.bed -o tmp/tmp.pickle -s 5 -w inter

echo ""
echo ""
echo ""
echo ""
echo "metaplot intra A+B"
python ../scripts/waffle-peaks.py --genomic_matrix tmp/tmp_bam_10kb.tsv --peaks data/peaks_prot.bed -o tmp/tmp.pickle -s 2 -w intra
