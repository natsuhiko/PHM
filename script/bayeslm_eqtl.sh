#!/bin/bash

# Bin ID
BIN=$1
# Bin width
W=$2

K=`expr '(' $BIN '-' 1 ')' "*" $W + 1`

FPKM=$3
VCF=$4
GENE=$5
OUT=$6

PHMDIR=`printenv PHMDIR`

for I in `sh $PHMDIR/script/extract.sh $GENE $BIN $W | awk '{print $1";"$2";"$3";"$5}'`
do
        IFS=';'
        set -- $I
	
	if [[ $4 == "+" ]]
	then
		TSS=$2
	else
		TSS=$3
	fi

	$PHMDIR/bin/bayeslm \
		--vcf $VCF \
		--normalised-count $FPKM \
		--feature-id $K \
		--window-chromosome $1 \
		--window-centre $TSS \
		--window-size 1000000 \
		--output $OUT
        
	K=`expr $K + 1`
done



