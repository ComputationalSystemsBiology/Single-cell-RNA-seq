#!/bin/bash

## Depends on gzip, bz2, XPath and wget

function compress {
	gunzip $1*.fastq.gz
	bzip2 $1*.fastq
}

# get experiment : given an ID download experiment xml
function get_experiment {
	wget -O $2 https://www.ebi.ac.uk/arrayexpress/xml/v3/experiments/$1
}

# write info : given an experiment xml write info file
function write_info {
	acc=$( xpath -q -e "experiments/experiment/accession" $1 | sed 's/<accession>//g' | sed 's/<\/accession>//g')
	organism=$( xpath -q -e "experiments/experiment/organism" $1 | sed 's/<organism>//g' | sed 's/<\/organism>//g')
	type=$( xpath -q -e "experiments/experiment/experimenttype" $1 | sed 's/<experimenttype>//g' | sed 's/<\/experimenttype>//g' )
	title=$(  xpath -q -e "experiments/experiment/bibliography/title" $1 | sed 's/<title>//g' | sed 's/<\/title>//g' )
	author=$(  xpath -q -e "experiments/experiment/bibliography/authors" $1 | sed 's/<authors>//g' | sed 's/<\/authors>//g' )
	desc=$(  xpath -q -e "experiments/experiment/description/text" $1 | sed 's/<text>//g' | sed 's/<\/text>//g' )
	cell_type=$(xpath -q -e "experiments/experiment/samplecharacteristic[category='cell type']/value" $1 | sed 's/<value>//g' | sed 's/<\/value>//g' )
	cell_line=$(xpath -q -e "experiments/experiment/samplecharacteristic[category='cell line']/value" $1 | sed 's/<value>//g' | sed 's/<\/value>//g' )
	
	info_file="${acc}.info"
	echo "	${acc}
Type : ${type}
Organism : ${organism}
Cells : ${cell_type}, ${cell_line}

Description : ${desc}

Publication : ${title}
Authors : ${author}" > ${info_file}
}

# get samples : given an ID download samples xml for parsing
function get_samples {
	wget -O $2 https://www.ebi.ac.uk/arrayexpress/xml/v3/experiments/$1/sample > samples.tmp
}
# get samples url : given a file name, extract samples URL for downloading
# parse XML file until samples URIs are found
function get_samples_url {
	xpath -q -e "experiment/sample/descendant::*[name='FASTQ_URI']/value" $1 | sed 's/<value>//g' | sed 's/<\/value>//g'
}

# Apply functions
exp_file="experiment.tmp"
sample_file="samples.tmp"
outDir=./Reads

# Process experiment infos
get_experiment $1 ${exp_file}
write_info ${exp_file}
rm ${exp_file}

# Download Samples
get_samples $1 ${sample_file}
urls=$(get_samples_url ${sample_file})
rm ${sample_file}


	# Check if URLs were found
if [ $(echo ${urls} | wc -l ) == 0 ]
then
	echo "ERROR : No URL found"
	exit  1
fi

# Create output directory
if [ ! -d ${outDir} ];
then
	mkdir ${outDir};
fi

cd ${outDir}

for url in ${urls}
do
	# Treat file name
	sample=$(echo ${url} | sed 's/[a-zA-Z0-9\.:]*\///g' | sed 's/\.[a-zA-Z0-9]*//g')
	file=$(echo ${url} | sed 's/[a-zA-Z0-9\.:]*\///g')
	if [ -e  ./${sample}*.fastq.bz2 ]
	then
		echo "${sample}.fastq.bz2 already exists ! skipping it...."
	
	else
	
		echo "Downloading ${sample}
url : ${url} "
		# Download
		wget ${url}
		# Compress
		compress ${sample}
	fi;
done

cd ..
