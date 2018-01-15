#!/bin/bash

# File
GZFILE=$1
# Bin ID
BIN=$2
# N rows in each bin
W=$3


A=`expr $BIN "*" $W`
B=$W
# N rows of GZFILE
N=`cat $1 | gunzip | wc -l | awk '{print $1}'`

if [ $A -ge $N ]
then
	B=`expr $B '-' $A '+' $N`
	A=$N
fi

cat $1 | gunzip | head -n $A | tail -n $B
