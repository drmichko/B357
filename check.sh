#!/usr/bin/bash
	running=0
	task=178
	core=6
	for(( j = 1; j < task; j++ )) ;  do
		./check.exe  $j   &
        	(( running++ ))
    		if (( running >= core )); then
       		  	wait -n  
       		 (( running-- ))
    		fi
	done
wait

