#!/bin/bash

# Removing files if they exists
# DOES NOT WORK
for file in "C.xyz_", "O.xyz_", "C.vxyz_", "O.vxyz_"
do
    echo "Checking for "$file
    if [ -e $file ]
    then
	echo "File found, removing..."
	rm $file
    fi
done
# DOES WORK
rm *_

echo "Parsing..."
typeset -i count=1
# Writting files
while read line
do
    if [ $count -gt 1 ]
    then
	string="96"
	typeset -i count2=1
	while [ $count2 -lt 97  ]
	do
	    string=$string" "$count2
	    count2=$count2+1
	done
	echo $string > "velocities.vxyz_"
    fi
    if [ $count -gt 2 ]
    then
	# Parsing GEOMETRY file
	x=$(echo $line | cut -d " " -f2 )
	y=$(echo $line | cut -d " " -f3 )	
	z=$(echo $line | cut -d " " -f4 )
	vx=$(echo $line | cut -d " " -f5 )
	vy=$(echo $line | cut -d " " -f6 )	
	vz=$(echo $line | cut -d " " -f7 )
	echo $vx" "$vy" "$vz >> "velocities.vxyz2_"
	echo $vx" "$vy" "$vz >> "velocities.vxyz_"
	# 3-35 = C , 36-98 = O 
	if [ $count -lt 35 ]
	then
	    echo $x" "$y" "$z >> "C.xyz_"
	else
	    echo $x" "$y" "$z >> "O.xyz_"
	fi
    fi
    count=$count+1
done < GEOMETRY.xyz
echo "DONE."
