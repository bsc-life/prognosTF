#! /bin/bash

echo "Simulate HiC"
# python Simulate_HiC.py

ret=`python -c 'import sys; print("%i" % (sys.hexversion<0x03000000))'`

echo ""
echo ""
echo ""
echo ""
echo "BAM2count"
if [ $ret -eq 0 ]; then
    python ../scripts/waffle-bam2count.py --bam data/fake.bam -r 10000 -C 8 -o data/data_bam_10kb.tsv -b data/biases3.pickle
else
    python ../scripts/waffle-bam2count.py --bam data/fake.bam -r 10000 -C 8 -o data/data_bam_10kb.tsv -b data/biases.pickle
fi

echo ""
echo ""
echo ""
echo ""
echo "metaplot intra A vs B"
python ../scripts/waffle-peaks.py --genomic_matrix data/data_bam_10kb.tsv --peaks data/peaks_protA.bed data/peaks_protB.bed -o tmp/tmp.tsv -s 5 -w intra

echo ""
echo ""
echo ""
echo ""
echo "metaplot inter A vs B"
python ../scripts/waffle-peaks.py --genomic_matrix data/data_bam_10kb.tsv --peaks data/peaks_protA.bed data/peaks_protB.bed -o tmp/tmp.pickle -s 5 -w inter --output_format pickle

echo ""
echo ""
echo ""
echo ""
echo "metaplot intra A+B"
python ../scripts/waffle-peaks.py --genomic_matrix data/data_bam_10kb.tsv --peaks data/peaks_prot.bed -o tmp/tmp.pickle -s 2 -w intra
