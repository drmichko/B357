#!/usr/bin/bash
	running=0
	task=178
	core=2
	for(( j = 1; j < task; j++ )) ;  do
		./zipcheck.exe  $j   &
        	(( running++ ))
    		if (( running >= core )); then
       		  	wait -n  
       		 (( running-- ))
    		fi
	done
wait

