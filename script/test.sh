#!/bin/bash

PHMDIR=`printenv PHMDIR`

if [ -z $PHMDIR ]
then
	echo ERROR: Please set PHMDIR enviroment variable.
	echo ERROR: Aborted.
	exit
fi

#################### USER DEFINITION #####################
# Bayes factors calculated by script/bayeslm1.sh
IN1=$PHMDIR/data/input1.gz
#
# Output directory of hm
OUT1=$PHMDIR/data/Stage1
#
# Bayes factors calculated by script/bayeslm2.sh
IN2=$PHMDIR/data/input2.gz
#
# Output directory of phm
OUT2=$PHMDIR/data/Stage2
#
# Log normalised read counts (binary double array)
FPKM=$PHMDIR/data/log.fpkm.chr22.bin
#
# VCF file (tabix indexed)
VCF=$PHMDIR/data/chr22.vcf.gz
#
# Peak annotation (tabix indexed)
PEAK=$PHMDIR/data/peaks.chr22.bed.gz
#
# Hyper-parameter estimate of variant level prior
VL=$PHMDIR/data/Stage1/variant_level.bin
#
# Prior probability that each peak is a QTL
PI1=$PHMDIR/data/Stage1/Pi1.bin
##########################################################

# First stage
sh $PHMDIR/script/bayeslm1.sh 1 2219 $FPKM $VCF $PEAK $IN1
R1=`cat $IN1 | gunzip | wc -l | awk '{print $1}'`
NF=`cat $PEAK | gunzip | wc -l | awk '{print $1}'`
$PHMDIR/bin/hm -i $IN1 -c I,S,S,S,S,S,S,S,C2,C3,S,S,B -r $R1 -f $NF -p -o $OUT1

# Second stage
sh $PHMDIR/script/bayeslm2.sh 1 2219 $FPKM $VCF $PEAK $VL $PI1 $IN2
R2=`cat $IN2 | gunzip | wc -l | awk '{print $1}'`
$PHMDIR/bin/phm -i $IN2 -c J,K,N4,B10 -r $R2 -f $NF -o $OUT2

