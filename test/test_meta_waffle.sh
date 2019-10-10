#! /bin/bash

echo "Simulate HiC"
python Simulate_HiC.py

echo ""
echo ""
echo ""
echo ""
echo "BAM2count"
python ../scripts/bam2count.py --bam data/fake.bam -r 10000 -C 8 -o tmp -b data/biases.pickle

echo ""
echo ""
echo ""
echo ""
echo "metaplot intra A vs B"
python ../scripts/peaks_x_interaction.py --genomic_matrix tmp/tmp_bam_10kb.tsv --bam data/fake.bam --biases data/biases.pickle --peaks data/peaks_protA.bed data/peaks_protB.bed -r 10000 -o tmp -s 5 -w intra

echo ""
echo ""
echo ""
echo ""
echo "metaplot inter A vs B"
python ../scripts/peaks_x_interaction.py --genomic_matrix tmp/tmp_bam_10kb.tsv --bam data/fake.bam --biases data/biases.pickle --peaks data/peaks_protA.bed data/peaks_protB.bed -r 10000 -o tmp -s 5 -w inter

echo ""
echo ""
echo ""
echo ""
echo "metaplot intra A+B"
python ../scripts/peaks_x_interaction.py --genomic_matrix tmp/tmp_bam_10kb.tsv --bam data/fake.bam --biases data/biases.pickle --peaks data/peaks_prot.bed -r 10000 -o tmp -s 2 -w intra
