#!/bin/bash
ens=$1
path=$2
file=$3

gcc -g motif_analysis.c packages/motifs_packages/*.c packages/C_packages/*.c -lm -o motif_analysis.out
if [ $? -eq 0 ];then
	echo "Compile complete!"
else
	echo "Compile failed!"
	exit 9
fi

./motif_analysis.out $ens $path $file &

WORK_PID=`jobs -l | awk '{print $2}'`
wait $WORK_PID

