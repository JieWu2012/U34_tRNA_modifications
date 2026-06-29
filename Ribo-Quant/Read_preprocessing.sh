# Analysis of ribosome (monosome/disome) profiling data


# Remove adapter sequence at the 3' end using fastx toolkit. (For wild type and ncs2Delp6D)

for i in *.fastq

do 

	echo $i

	a=${i/\.fastq/}

	fastx_clipper -a CTGTAGGCACCATCAAT -l 15 -c -Q 33  -i $i   -o $a"_clipped.fastq"

done

# Trim addtional 4 nt at the 3' end. 
# If there are additional nt at the 5' end, the reads need to be further trimmed using "-f". See manual of fastx toolkit. 

for i in *_clipped.fastq

do 

	echo $i

	a=${i/\.fastq/}

	fastx_trimmer -m 15 -t 4 -i $i -o $a"_trimmer_mt15.fastq"

	rm $i

done


# Fastqc

for i in *_trimmer_mt15.fastq

do 
	echo $i

	fastqc $i

done

# Remove non protein-coding reads

mkdir ncRNA_rm

ncRNA_STAR=SacCer3_ncRNA_STAR

for i in data/*.fastq 

do 

	echo $i

	a=${i/\.fastq/}

	b=${a/*\//}

	STAR 	--runThreadN 20 \
			--genomeDir $ncRNA_STAR \
			--readFilesIn  \
			--outFilterMultimapNmax 30 \
			--outSAMmultNmax 1 \
			--outFilterMismatchNoverReadLmax 0.1  \
			--alignEndsType EndToEnd \
			--outFileNamePrefix ncRNA_rm/$b"_toncRNA_"  \
			--outReadsUnmapped Fastx \
			--outSAMtype BAM Unsorted 
done

# Change file name

for i in *mate1;do echo $i; a=${i/\mate1/fastq}; mv $i $a; done;
