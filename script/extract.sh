
A=`expr $2 "*" $3`
B=$3
N=`zwc $1 | cut -f 1`
if [ $A -ge $N ]
then
	B=`expr $B - $A + $N`
	A=$N
fi
zcat $1 | head -n $A | tail -n $B
