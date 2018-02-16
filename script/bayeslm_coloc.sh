#!/bin/bash

# Bin ID
BIN=$1
# Bin width
W=$2

K=`expr '(' $BIN '-' 1 ')' "*" $W + 1`

FPKM_EQTL=$3
FPKM=$4
VCF_EQTL=$5
VCF=$6
GENE=$7
PEAK=$8
VL0=$9
PI1_EQTL=${10}
PI1_CAQTL=${11}
OUT=${12}

PHMDIR=`printenv PHMDIR`

for I in `sh $PHMDIR/script/extract.sh $GENE $BIN $W | awk '{print $1";"$2";"$3";"$4}'`
do
        IFS=';'
        set -- $I
	
	# TSS
	if [[ $4 == "+" ]]
	then
		TSS=$2
	else
		TSS=$3
	fi
	
	# Window
	A=`expr $TSS '-' 500000`
	B=`expr $TSS '+' 500000`
	if [ $A -le 0 ]; then A=1; fi

	# N peaks
	PEAKID=0
	PEAKID=`cat $PEAK | gunzip | awk -v CHR=$1 -v A=$A -v B=$B '$1==CHR && int($2+$3)/2>A && int($2+$3)/2<B {print NR}' | head -n 1` 

	$PHMDIR/bin/bayeslm \
		--vcf $VCF_EQTL --vcf2 $VCF \
		--normalised-count $FPKM_EQTL --normalised-count2 $FPKM \
		--feature-id $K --feature-id2 $PEAKID \
		--feature-bed $PEAK \
		--window-chromosome $1 \
		--window-centre $TSS \
		--window-size 1000000 \
		--variant-level $VL0 --variant-level2 $VL0 \
		--feature-level $PI1_EQTL --feature-level2 $PI1_CAQTL \
		--output $OUT

	K=`expr $K + 1`
done



