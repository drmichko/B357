#!/bin/bash

base=$1

taille=$( bc <<< "2^32" )
small=0
large=0
total=0
cs=0
class=12
while read sharp num stab orb last etc
do
	#echo $num $stab $orb $last
	orb=${orb##*=}
	size=$( bc <<< "$orb*$base"  )
	if [ $size -lt $taille ]; then
		if [ ! -f "data/$num-3-3-7.zip" ] ; then 
			echo num=$num orb=$orb size=$size 
		fi
		let small++
		let cs+=orb
		let total+=size
	else
		echo $num $size large
		let large++
	fi
	class=$( bc <<< "$class+$orb" )
done < <( grep num= info-3-5-7.txt )
echo small : $small  count=$cs  size=$total
echo large : $large  size=$( bc <<< "$large * $taille")
bc <<< "178 * $taille"
echo classe $class
