#!/bin/bash

# Bin ID
BIN=$1
# Bin width
W=$2

K=`expr '(' $BIN '-' 1 ')' "*" $W + 1`

FPKM=$3
VCF=$4
PEAK=$5
OUT=$6

PHMDIR=`printenv PHMDIR`

for I in `sh $PHMDIR/script/extract.sh $PEAK $BIN $W | awk '{print $1";"$2";"$3}'`
do
        IFS=';'
        set -- $I
	
	PEAKCENTRE=`expr '(' $2 '+' $3 ')' / 2`

	$PHMDIR/bin/bayeslm \
		--vcf $VCF \
		--normalised-count $FPKM \
		--feature-bed $PEAK \
		--feature-id $K \
		--window-chromosome $1 \
		--window-centre $PEAKCENTRE \
		--window-size 1000000 \
		--output $OUT
        
	K=`expr $K + 1`
done



