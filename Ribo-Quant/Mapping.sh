# Map reads to CDS and generate periodicity plot. 

CDSbtIndex=SacCer3_CDS_ex21

for i in rRNA_rm/*.fastq

do 
	echo $i
	./rp_analysis.sh -M Mapping \
					 -Q $i  \
					 -T $CDSbtIndex   \
					 -O mapping/  \
					 -X 34 \
					 -I 16 \
					 -E 21
done