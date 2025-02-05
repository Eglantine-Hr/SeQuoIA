#!/bin/bash

#set -u
#set -o pipefail

echo "STARTING BASH SCRIPT"



#inputPath=$(dirname $1)
inputPath=$(echo $1)
outputPath=$(echo $2)
DATASET=$( echo $3) 
log_parsimony=$( echo $4)
#Num_process=$( echo $5)

currentdir=$(pwd)
cd ${currentdir}



echo ${currentdir}
echo ${DATASET}





CloneName=$(echo $(basename $inputPath)) 
inputPath=$(echo $(dirname $inputPath))

echo "all VAR"
echo ${CloneName}
echo ${inputPath}
echo ${outputPath}

#Getting clone size (for computational time optimization) 
CloneSize=$(echo $CloneName | grep -o -P '(?<=_n).*(?=_Id)')
CloneSize=$(expr $CloneSize + 0)

export QT_QPA_PLATFORM=offscreen

for currentChain_fasta in $(find ${inputPath}/${CloneName} -name '*Sequences.fasta' -printf "%P\n"  ) #getting fasta of heavy and light chains 
do
	if test -f "${inputPath}/${CloneName}/${currentChain_fasta}"; then
		
		echo $currentChain_fasta
		ChainType=$(echo "${currentChain_fasta%%_*}")  #before first _, should be Heavy or Light 
	
		#File management
		mkdir -p ${outputPath}/${CloneName}/${ChainType}
		cp ${inputPath}/${CloneName}/${currentChain_fasta} ${outputPath}/${CloneName}/${ChainType}
		cp ${inputPath}/${CloneName}/${ChainType}_isotypemap.txt ${outputPath}/${CloneName}/${ChainType}/isotypemap.txt
		cd ${outputPath}/${CloneName}/${ChainType}
		echo "directories created"
	
	
		# we deduplicate the sequences and convert to phylip alignment format, and also determine sequence abundances 
		deduplicate ${currentChain_fasta} --root Ancestor --abundance_file abundances.csv --idmapfile idmap.txt > dedup.phylip
		echo "deduplication done "
	
		# We use PHYLIP’s dnapars program to generate a set of maximum parsimony trees (https://matsengrp.github.io/gctree/quickstart.html)
		if (( $CloneSize > 50 )); then
			#mkconfig --quick dedup.phylip dnapars > dnapars.cfg 
			#echo "quick argument"
			mkconfig dedup.phylip dnapars > dnapars.cfg 
		else
			mkconfig dedup.phylip dnapars > dnapars.cfg 
		fi
		
		#here we replace absolute path with relative due to dnapars character limitation 
		sed "s#${outputPath}/${CloneName}/${ChainType}/##" dnapars.cfg  > dnapars_modif.cfg   #s substitute / to seperate search and replacement, if in path, use #
		phylip dnapars < dnapars_modif.cfg > dnapars.log 
		wait
		echo "parsimony trees done"
	
		#generating .p forest and computing galton params => pb unwanted outputs as well, hence --outbase to differentiate them  
		gctree infer outfile abundances.csv --root Ancestor --frame 1 --idmap idmap.txt --isotype_mapfile isotypemap.txt --idlabel --summarize_forest --outbase default --verbose || true
		wait
		
		# Extracting galton parameters from the summary log (select line containing keyword params and deleting irrelevant characters)
		
		sed -n '/^Param/p' default.forest_summary.log | sed 's/[a-zA-Z]\|[*:]\|[[:space:]]//g' > default.forest_summary2.log 
		
		#Parsing forest and generating outputs 
		python ${currentdir}/03_Script/02_Phylogeny/01_Generating_Forest/02_Parsing_Forest_Outputs.py
		echo "parsing trees done"
			
	 	

	fi
	

	
	cd ${currentdir} # going back to root
		
done 

     

