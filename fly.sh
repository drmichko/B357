#!/usr/bin/bash
	running=0
	task=178
	core=6
	for(( j = 1; j < task; j++ )) ;  do
		./onthefly.exe -t $j  -p72
	done
wait

