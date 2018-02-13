#!/bin/bash

PHMDIR=`printenv PHMDIR`

if [ -z $PHMDIR ]
then
	echo ERROR: Please set PHMDIR enviroment variable.
	echo ERROR: Aborted.
	exit
fi

#################### USER DEFINITION #####################
# Bayes factors for caQTLs
IN_CAQTL=$PHMDIR/data/input1.gz
if [ ! -e $IN_CAQTL ]; then echo ERROR: Run $PHMDIR/script/test.sh first.; exit; fi
#
# Bayes factors for eQTLs
IN_EQTL=$PHMDIR/data/input_eqtl.gz
#
# Output directory of hm for eQTL mapping
OUT_EQTL=$PHMDIR/data/Eqtl
#
# Output directory of hm for caQTL mapping
OUT_CAQTL=$PHMDIR/data/Caqtl
#
# Bayes factors calculated by script/bayeslm2.sh
IN_COLOC=$PHMDIR/data/input_coloc.gz
#
# Output directory of coloc
OUT_COLOC=$PHMDIR/data/Coloc
#
# Log normalised read counts of RNA-seq from gEUVADIS 372 European samples (binary double array)
FPKM_EQTL=$PHMDIR/data/log.fpkm.geu.372.chr22.bin
#
# Log normalised read counts of in-house ATAC-seq (binary double array)
FPKM=$PHMDIR/data/log.fpkm.chr22.bin
#
# VCF file for gEUVADIS 372 European samples (tabix indexed)
VCF_EQTL=$PHMDIR/data/chr22.geu.372.vcf.gz
#
# VCF file for GBR 100 LCL samples (tabix indexed)
VCF=$PHMDIR/data/chr22.vcf.gz
#
# Peak annotation (tabix indexed)
PEAK=$PHMDIR/data/peaks.chr22.bed.gz
#
# Gene annotation (tabix indexed)
GENE=$PHMDIR/data/genes.chr22.bed.gz
#
# Vector of 0s suggesting variant level prior if noninformative
VL0=$PHMDIR/data/variant_level0.bin
#
# Prior probability that each peak is a caQTL
PI1_CAQTL=$PHMDIR/data/Caqtl/Pi1.bin
#
# Prior probability that each gene is an eQTL
PI1_EQTL=$PHMDIR/data/Eqtl/Pi1.bin
##########################################################

# First stage (eQTL)
NF1=`cat $GENE | gunzip | wc -l | awk '{print $1}'`
sh $PHMDIR/script/bayeslm_eqtl.sh 1 $NF1 $FPKM_EQTL $VCF_EQTL $GENE $IN_EQTL
R1=`cat $IN_EQTL | gunzip | wc -l | awk '{print $1}'`
$PHMDIR/bin/hm -i $IN_EQTL  -c I,S,S,S,S,S,S,S,S,S,S,S,B -r $R1 -f $NF1 -o $OUT_EQTL -p

# First stage (caQTL)
NF2=`cat $PEAK | gunzip | wc -l | awk '{print $1}'`
R2=`cat $IN_CAQTL | gunzip | wc -l | awk '{print $1}'`
$PHMDIR/bin/hm -i $IN_CAQTL -c I,S,S,S,S,S,S,S,S,S,S,S,B -r $R2 -f $NF2 -o $OUT_CAQTL -p

# Coloc
sh $PHMDIR/script/bayeslm_coloc.sh 1 $NF1 $FPKM_EQTL $FPKM $VCF_EQTL $VCF $GENE $PEAK $VL0 $PI1_EQTL $PI1_CAQTL $IN_COLOC
R3=`cat $IN_COLOC | gunzip | wc -l | awk '{print $1}'`
$PHMDIR/bin/phm -i $IN_COLOC -c J,K,B6 -r $R3 -o $OUT_COLOC




