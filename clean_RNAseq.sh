#!/bin/bash

#clean_RNAseq.sh
#script to clean RNA-seq reads prior to assembly
#Megan Supple
#10 Jan 2014

#usage:  clean_RNAseq.sh <clean.RNAseq.input>
#<clean.RNAseq.input> is a text file with entries as in example input file (see example.clean.RNAseq.input or example at the end of this script

#requires:  Trimmomatic, FLASh, FastQC

#produces output in the current directory:
	#clean.log			-- output from stdout and stderr
	#resource.log			-- log of computational resources used
	#*.fastqc,*fastq.gz		-- FastQC reports
        #<sample>.pair1.fastq.gz        --**cleaned reads, pair 1**
        #<sample>.pair1.fastq.gz        --**cleaned reads, pair 2**
        #<sample>.single.fastq.gz       --**cleaned reads, unpaired** 



#read in input file
source $1
#set output file for stdout
exec &>clean.log


echo "cleaning RNA-Seq data" > $resourcelog

#loop over each paired end sample
for ((a=0; a<${#id[@]}; a++))
	do
		echo -e "\n\ncleaning sample ${id[a]}"	
		echo -e "\n\ncleaning sample ${id[a]}" >> $resourcelog
	
		#adapter trim, merge completely overlapping reads, and quality filtering with Trimmomatic
		echo -e "\ntrimmomatic"
		echo -e "\ntrimmomatic" >> resource.log
		if [ $score = "-phred33" ]
		  then
			/usr/bin/time -v -o $resourcelog -a java -classpath $path_trimmomatic org.usadellab.trimmomatic.TrimmomaticPE -threads $threads $score $files/${file1[a]} $files/${file2[a]} ${id[a]}.trim.p1.fastq ${id[a]}.trim.u1.fastq ${id[a]}.trim.p2.fastq ${id[a]}.trim.u2.fastq ILLUMINACLIP:$adapter:$seedMismatches:$palindromeClipThreshold:$simpleClipThreshold HEADCROP:$headCropLen SLIDINGWINDOW:$windowSize:$windowQuality LEADING:$leadQuality TRAILING:$trailQuality
		elif [$score = "-phred64" ]
		  then
                        /usr/bin/time -v -o $resourcelog -a java -classpath $path_trimmomatic org.usadellab.trimmomatic.TrimmomaticPE -threads $threads $score $files/${file1[a]} $files/${file2[a]} ${id[a]}.trim.p1.fastq ${id[a]}.trim.u1.fastq ${id[a]}.trim.p2.fastq ${id[a]}.trim.u2.fastq ILLUMINACLIP:$adapter:$seedMismatches:$palindromeClipThreshold:$simpleClipThreshold HEADCROP:$headCropLen SLIDINGWINDOW:$windowSize:$windowQuality LEADING:$leadQuality TRAILING:$trailQuality TOPHRED33
		else
		  echo "$score is not a valid quality score for trimmomatic. Should be -phred33 or -phred64"; exit 1
                fi

		#merge overlapping reads with FLASH
		echo -e "\nflash"
		echo -e "\nflash" >> $resourcelog
		/usr/bin/time -v -o $resourcelog -a flash -o ${id[a]} -t $threads -m $minOverlap -x $maxMismatchDensity -p 33 -M $maxOverlap ${id[a]}.trim.p1.fastq ${id[a]}.trim.p2.fastq
		
		#merge single end fastqs
		echo -e "\ncat"
		echo -e "\ncat" >> $resourcelog
		/usr/bin/time -v -o $resourcelog -a cat ${id[a]}.trim.u1.fastq ${id[a]}.trim.u2.fastq ${id[a]}.extendedFrags.fastq > ${id[a]}.singlemerge.fastq
		
		#remove small reads with trimmomatic
		echo -e "\ntrimmomatic2"
		echo -e "\ntrimmomatic2" >> $resourcelog
		/usr/bin/time -v -o $resourcelog -a java -classpath $path_trimmomatic org.usadellab.trimmomatic.TrimmomaticPE -phred33 ${id[a]}.notCombined_1.fastq ${id[a]}.notCombined_2.fastq ${id[a]}.pair1.fastq ${id[a]}.single1.fastq ${id[a]}.pair2.fastq ${id[a]}.single2.fastq MINLEN:$minlen 
		/usr/bin/time -v -o $resourcelog -a java -classpath $path_trimmomatic org.usadellab.trimmomatic.TrimmomaticSE -phred33 ${id[a]}.singlemerge.fastq ${id[a]}.single3.fastq MINLEN:$minlen
		
		#merge single end fastqs, again
		echo -e "\ncat2"
		echo -e "\ncat2" >> $resourcelog
		/usr/bin/time -v -o $resourcelog -a cat ${id[a]}.single1.fastq ${id[a]}.single2.fastq ${id[a]}.single3.fastq > ${id[a]}.single.fastq
		
		#FastQC final files
                echo -e "\nFastQC"
		fastqc -t $threads ${id[a]}.pair1.fastq 
		fastqc -t $threads ${id[a]}.pair2.fastq 
		fastqc -t $threads ${id[a]}.single.fastq 

                #gzip final files and remove intermediates
                gzip ${id[a]}.pair1.fastq &
                gzip ${id[a]}.pair2.fastq &
                gzip ${id[a]}.single.fastq &

		rm ${id[a]}.trim.u1.fastq ${id[a]}.trim.u2.fastq ${id[a]}.trim.p1.fastq ${id[a]}.trim.p2.fastq ${id[a]}.notCombined_1.fastq ${id[a]}.notCombined_2.fastq ${id[a]}.extendedFrags.fastq ${id[a]}.hist ${id[a]}.histogram ${id[a]}.single1.fastq ${id[a]}.single2.fastq ${id[a]}.single3.fastq ${id[a]}.singlemerge.fastq 

		echo -e "\ndone cleaning sample $id[a]}\n\n"				

	done

echo DONE!!!
echo | mutt -s "RNA-seq cleaning from $1 complete" $email

