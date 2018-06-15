#!/bin/bash

#=======================#
# Equation Of State     #
#  Data Analysis Script #
#=======================#

#-------------
# Description
#-----------------------------------------------
#
#
#------------------------------------------------

#-------
# Input
#------------------------------------------------------------------------------------------
# -p [Path/To/Data/Folder]      : (Optionnal) Path to the where data are
# -r [Path/Where/To/Store/Data] : (Optionnal) Folder where to store the data for analysis
#------------------------------------------------------------------------------------------

#--------
# Output
#----------------------------------------------------------------------------
# - For each structrure:
#   - Data of the enthalpie as a function of the pressure
#   - The Bulk Modulus and its derivative
#   - The Volume
# - Graphs different functionnals against each other
# - Graphs all structure done with the same functionnal on the pessure range
#-----------------------------------------------------------------------------

# Start
echo -e "\n----------------------------------------------------"
echo -e " Equation of State - Data Analysis Script starting"
echo -e "----------------------------------------------------\n"

# Loading Usage Function

# Input Check

# Location
pwd=$(pwd)

# House for data files
data=$pwd"/data"
# If already exists, cleans it (new is already better)
if [ -e $data ]
then
    rm -r $data
fi
# Creating data 
mkdir $data

# Aggregating data file in a single data folder
#--------------------------------------------------
# Recursive exploration function
function explore {
   # Listing all files
    nohup ls > temp_explore 2>/dev/null
    # Listing all files ending with .out if no, creating empty file
    nohup ls *.out > temp_out 2>/dev/null || touch temp_out 
    # If there at .out files...
    if [ -s temp_out ]
    then
	# We move them all to $data
	while read line
	do
	    # New is always better...
	    if [ -e $data"/"$line ]
	    then
		rm $data"/"$line
	    fi
	    cp $line $data
	done < temp_out
    else
    # Otherwise we just go further up the tree, if we can.
	while read line
	do 
	    if [ -d $line ]
	    then
		cd $line
		explore
	    fi
	done < temp_explore
    fi
    # If we reach this point, we start down the tree
    if [[ $pwd != $(pwd) ]]
    then
	cd ..
    fi
    # And we remove temp files as we go
    if [ -e temp.out ]
    then
	rm temp.out
    fi
    if [ -e temp_explore ]
    then
	rm temp_explore
    fi
}

# Launching the recursion
explore

# Extracting Data
#----------------------

# Moving to data file
cd $data

# Removing non necessary files
rm temp* || echo "loutre" > /dev/null

# Listing files
ls > temp

# Identifiying the different structures and creating folder
while read line
do
    struct=$(echo $line | cut -d "." -f1)
    if [ -e $struct ]
    then
	rm -r $struct
    fi
    mkdir $struct
    cp $struct* $struct"/"
done < temp

# Getting energy and volume from all files
# - Loop on the files
while read line
do
    if [[ $line != "temp" ]]
    then
	grep bfgs $line > step
	grep volume $line > volume
	grep energy $line > energy
    fi
done < temp



# Removing temp file
rm temp