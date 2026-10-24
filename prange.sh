#!/bin/bash
 
debut=$( hostname | tr -cd '[0-9]' )

for (( num=debut; num < 179; num+=2 ))
do
	file=data/$num-3-3-7.zip
	if [ ! -f $file ]; then
		gzip -kd $file.gz
	fi
	ls -lart  $file
	./prange.exe  -t$num  -w
done
