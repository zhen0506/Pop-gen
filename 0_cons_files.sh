#!/usr/bin/bash

#concatenate fastq files from the same genotype
file=$1
sample=""
raw="/orange/seonghee/zhen/popgen/raw_read/"
tail -n +2  UC_samplelist.txt | while read SRR SAM ID NAME TAXA; do
	if [[ ${sample} != ${ID} ]]; then
		if [[ -e ${arr1[0]} ]]; then
			echo "${arr1[@]} ${arr2[@]}"
			cat ${arr1[@]} | gzip > ${raw}${ID}_1.fastq.gz
                	cat ${arr2[@]} | gzip > ${raw}${ID}_2.fastq.gz
		fi
		sample=$ID
		arr1=(${raw}${SRR}"_1.fastq")
		arr2=(${raw}${SRR}"_2.fastq")
	else 
                arr1+=(${raw}${SRR}"_1.fastq")
                arr2+=(${raw}${SRR}"_2.fastq")
	fi
	done

