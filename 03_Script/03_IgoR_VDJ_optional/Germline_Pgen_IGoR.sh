#!/bin/bash

set -u
set -o pipefail

currentdir=$(pwd)
echo ${currentdir}

#cd $(dirname $1)
echo $(dirname $1)
cd $1
#ls



outputPath=$(dirname $2)
echo ${outputPath}

organism=$3
echo $organism

#cd 05_Output/01_BCRseq_processing/03_UCA_MRCA_inference_annotation/{wildcards.DATASET}/{wildcards.SAMPLE} 


echo 'CloneName, Chain, Pgen, NCAseqAbsolute' > Pgen_VDJ_summary.txt

for CloneName in `find . -type d -mindepth 1 -maxdepth 1` #excluding current folder
do
	#cd ../ 
   	#cd $CloneName
   	echo ${CloneName}
   	#echo "step1"
   	cd ${CloneName}
   	CloneName=${CloneName#*./}
   	#NCAabsolute= $(${CloneName}/Heavy_Absolute_NCA.txt)
   	#NCAabsolute="Heavy_Absolute_NCA.txt"
#   	#echo "$NCAabsolute"
   	
   	ls
   	
   	if test -f "Heavy_Absolute_Germline.txt"; then 
		NCAabsolute=$(cat "Heavy_Absolute_Germline.txt") #fetching NCA sequence
		#echo ${NCAabsolute}
   		PGEN=$(igor-compute_pgen ${organism} heavy_naive ${NCAabsolute:1:-1}) #computing probas, we remove quotation marks
   	 	PGEN=${PGEN#*$'\n'} #${PGEN#*txt}#https://stackoverflow.com/questions/428109/extract-substring-in-bash           ##PGEN=${PGEN#*$'\n'} #pk marche?
   	 	#echo ${PGEN} 
   	 	echo "${CloneName}, Heavy, ${PGEN}" >> ${currentdir}/${outputPath}/Pgen_VDJ_summary.txt #adding data to summary table ##HOW TO IMRPOVE THIS STEP? $3_
		#break 
   	fi
   	#"Support for light chains in command line is not ready yet due to the lack of genomic templates and suitable model."
#  	if test -f "Light_Absolute_Germline.txt"; then
#   		NCAabsolute=$(cat "Light_Absolute_Germline.txt") #fetching NCA sequence 
#   		PGEN=$(igor-compute_pgen human light $NCAabsolute) #computing probas
#   	 	PGEN=${PGEN#*$'\n'} #https://stackoverflow.com/questions/428109/extract-substring-in-bash
#   	 	echo "${CloneName}, Light, ${PGEN}}" >> ../Pgen_VDJ_summary.txt #adding data to summary table  ##, ${NCAabsolute}
#   	fi
	cd ../
done
 
