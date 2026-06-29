#!/bin/bash
#!/bin/awk -f



help_info(){
  echo "
 **********************************************************************
*   _____    _____      _____     _____      _____   __    ___   ____  *
*  |  __  \ |  __  \   |  __  \  /  __  \   /  ___| |  |  /  /  /  __| *
*  | |__\ | | |__\ |   | |__\ | |  /  \  | |  /     |  | /  /  |  /    *
*  |   ___| |  ____|   |   ___| | |    | | |  |     |  |/  /    \  \   *
*  |   \    | |        |   \    | |    | | |  |     |      \     \  \  *
*  | |\ \   | |        | |\ \   |  \__/  | |  |___  |  |\   \   __/  / *
*  |_| \_\  |_|        |_| \_\   \ ____ /   \ ____| |__| \___\ |____/  *
*                                                                      *
*                                                                      *
*          Copyright by Jie Wu   2017.01.23                            *
*                                                                      *
 **********************************************************************
  "
  echo "./rp_analysis.sh : tool for analyzing Ribosome profiling data and make codon plot"
  echo "Usage: ./rp_analysis.sh -M [Mode] [Options]"

  echo "
             [ -M|--mode ]                = Analyzation mode, choose one from \"Mapping\", \"CodonPlot\". 
             [ -H|--help|-h  ]            = Help information.

      1) When [Mode] is \"Mapping\", Options: 
      
             [ -Q|--fastq ]               = Input fastq file.
             [ -T|--transcriptom-index ]  = Prefix of bowtie index of cds transcriptome. 
             [ -O|--output ]              = Name of output folder. 
             [ -X|--max_length ]          = Maximum read length to check the frame information.
             [ -I|--min_length ]          = Minimum read length to check the frame information.
             [ -E|--extend ]              = Length extended flanking annotated CDS, which was built as the transcriptome index.
      2) When [Mode] is \"CodonPlot\", Options:

             [ -A|--fasta ]               = Fasta file that were used to build the transcriptome index.
             [ -D|--index ]               = Index file (7 columns seperated by tab including sam/bed6 file name, sample name, replicates, read length, frame, A site, class). eg. \"CodonPlot.index\".
             [ -O|--output ]              = Name of output folder.
             [ -E|--extend ]              = Length extended flanking annotated CDS, which was built as the transcriptome index.       
             [ -F|--fexcluded ]           = Codon length excluded at the 5' end of each transcript when plotting the codon plot. 
             [ -B|--bexcluded ]           = Codon length excluded at the 3' end of each transcript when plotting the codon plot. 
  " 
  echo "
  eg. 1) ./rp_analysis.sh -M Mapping -Q fastq_file -T index_prefix -O output -X 31 -I 28 -E 18
      2) ./rp_analysis.sh -M CodonPlot -A  transcriptome_fasta_file -D index_file -O output -E 18 -F 15 -B 15
" 
  if ! [ -x "$(command -v R)" ];then
     echo -e "Warnings: \n\tR is not installed in this server so pdf files can't be made. Please install it or use the files from the results such as \"all_codon_Asite_nor_wt.txt\", \"all_codon_Asite_offset.txt\", to make the plot in Rstudio."

  fi

}


if [ $# -eq 0 ]
then
  help_info
fi

TEMP=`getopt -o M:Q:T:O:L:A:X:I:D:E:S:C:F:B:Hh --long mode:,fastq:,transcriptome_index:,output:,length:,fasta:,max_length:,min_length:,index:,extend:,fexcluded:,bexcluded:,offset:,codons:,help  -- "$@"`
eval set -- "$TEMP"


while true
do
        case "$1" in
    -M|--mode)
          ARG_mode=$2
          shift 2
        ;;
    -Q|--fastq)
          ARG_fastq=$2
          shift 2
        ;;
    -T|--transcriptom-index)
          ARG_transcriptome_index=$2
          shift 2
        ;;
    -O|--output)
          ARG_output=$2
          shift 2
        ;;
    -A|--fasta)
          ARG_fasta=$2
          shift 2
        ;;
    -X|--max_length)
          ARG_max_length=$2
          shift 2
        ;;
    -I|--min_length)
          ARG_min_length=$2
          shift 2
        ;;
    -D|--index)
         ARG_index=$2
         shift 2
        ;;
    -E|--extend)
          ARG_extend=$2
          shift 2
        ;;
    -F|--fexcluded)
          ARG_exclude_5end=$2
          shift 2
        ;;
    -B|--bexcluded)
          ARG_exclude_3end=$2
          shift 2
        ;;
    -H|-h|--help)
          help_info
          shift
        ;;
    --)
          shift
          break
        ;;
    *)
          echo "Internal error!"
          #help_info
          exit 1
        ;;
        esac
done
#echo $ARG_fastq
if [[ $ARG_output != "" ]]; then

   if [[ ! -d $ARG_output ]]; then
  
      mkdir $ARG_output

   fi
fi
chmod +x codon_plot.r
chmod +x wave_plot.r
chmod +x frame_plot.r

if [[ $ARG_mode == "Mapping" ]]; then
    
    Start_mapping=`date +%s`
    
    echo -e "\nMapping......"
    name=${ARG_fastq/\.fastq}
    name=${name/*\//}
   bowtie -p 16 -v 1 -m 1 --no-unal  --norc  --best --strata -t $ARG_transcriptome_index -q $ARG_fastq -S $ARG_output/$name"_mapped.sam"  2>$ARG_output/$name"_mapped.log"
   

   awk  'BEGIN{FS="\t";OFS="\t"} $2==0{
   
     sub(/M/,"",$6)
     if($0~/MD:Z:0[A-Z]/){              
         p=substr($10,2,length($10))
         print $3,$4-1+1,$4-1+$6,"@"$1"|"p,"0","+"    # reads with the first unmapped nucleotide.
       }else{
         print $3,$4-1,$4-1+$6,"@"$1"|"$10,"0","+"
     }
   
   }' $ARG_output/$name"_mapped.sam" > $ARG_output/$name".bed6" #produce bed6 file for codon plot.  
    
    echo -e "\nProducing $name frame file..."
    awk 'BEGIN{OFS="\t";FS="\t"}$2==0{                                                                                                                                                                                                                                                                                                                                                                                                                                                    

      sub(/M/,"",$6)
      if($6>="'$ARG_min_length'" && $6<="'$ARG_max_length'"){
          print $4-"'$ARG_extend'"-1,$6
      }
    
    }' $ARG_output/$name"_mapped.sam" | awk 'BEGIN{OFS="\t";FS="\t"}{hash[$0]++}END{for(i in hash){print i,hash[i]}}'  > $ARG_output/$name"_frame.txt" # read counts in frame
    #echo $ARG_fastq
    if [[  -f $ARG_output/$name"_frame_star.txt" ]]; then
      
       rm $ARG_output/$name"_frame_star.txt"
    
    fi

    for i in $(eval echo "{$ARG_min_length..$ARG_max_length}"); do   # traverse all the read lengths provided. 
      echo -e "\n Frame information for $i nt reads" >>  $ARG_output/$name"_frame_star.txt"
      awk 'BEGIN{OFS="\t";FS="\t"}$1>=-20 && $1<=50 && $2=="'$i'" && NR==FNR{
        if(hash[$2]<$3){
          hash[$2]=$3
        }
      }$1>-20 && $1<50 && $2=="'$i'" && NR!=FNR{
        star=sprintf("%*s",$3*100/hash[$2],"*")
        gsub(/ /,"*",star)
        print $0"\t"star
      }' $ARG_output/$name"_frame.txt"  $ARG_output/$name"_frame.txt" | sort -k1,1 -n  >> $ARG_output/$name"_frame_star.txt"  #show the frame around start codon

    done
    if ! [ -x "$(command -v R)" ];then
       echo -e "Warnings: \n\tR is not installed in this server so pdf files can't be made. Please install it or use the files from the results such as files with suffix \"_frame.txt\" to make the plot in Rstudio."

    else
       R --vanilla --slave --args $ARG_output/$name"_frame.txt" $ARG_max_length $ARG_min_length <frame_plot.r  # making the frame plot. see the pdf file. 
    fi
    End_mapping=`date +%s`
    Time_mapping=$((End_mapping-Start_mapping))
    echo "Took $((Time_mapping/60)) mins to map and produce frame plot"
fi


if [[ $ARG_mode =~ "CodonPlot" ]]; then
 Start_plotting=`date +%s`
 if [[  -f $ARG_output/all_codon_Asite.txt ]]; then
  
   rm $ARG_output/all_codon_Asite.txt

 fi
 if [[  -f $ARG_output/all_codon_Asite_offset.txt ]]; then
  
   rm $ARG_output/all_codon_Asite_offset.txt

 fi
 echo -e "\n Assigning ribosome occupancy of A site...... "
 IFS=$'\n'
 fasta_filename=${ARG_fasta/\.fa*/}
 fasta_filename=${fasta_filename/*\//}

 for i in `cat $ARG_index`;do

   OLD_IFS="$IFS"
   IFS=$'\t'
   array=($i)
   IFS="$OLD_IFS"
   if [[ ${array[0]} =~ ".sam" ]]; then
         #statements
      file_name=${array[0]/\.sam/}
      file_name1=${file_name/*\//}
      awk  'BEGIN{OFS="\t";FS="\t"} $2==0{
      
        sub(/M/,"",$6)
        if($0~/MD:Z:0[A-Z]/){              
          p=substr($10,2,length($10))
          print $3,$4-1+1,$4-1+$6,"@"$1"|"p,"0","+"    # reads with the first unmapped nucleotide.
        }else{
          print $3,$4-1,$4-1+$6,"@"$1"|"$10,"0","+"
        } 
      
      }' $file_name".sam" > $file_name".bed6"  # producing bed file. 

   elif [[ ${array[0]} =~ ".bed6" ]]; then
     #statements   
       file_name=${array[0]/\.bed6/}  
       file_name1=${file_name/*\//}

   fi
   

   echo -e "\n\tProcessing ${array[0]} ......"
    awk 'BEGIN{
      OFS="\t"
      split("'${array[3]}'",length_array,",")
      for(len=1;len<=length(length_array);len++){
        hash_len[length_array[len]]=1
      } # read length in hash
      split("'${array[4]}'",frame_array,",")
      split("'${array[5]}'",A_site_array,",")
      for(i=1;i<=length(frame_array);i++){
        split(frame_array[i],frame_length,":")
        split(A_site_array[i],A_site_length,":")
        for(j=1;j<=length(frame_length);j++){
          frame_length_hash[length_array[i]"\t"frame_length[j]]=A_site_length[j]
          
        }
      } # read length and frame in the hash
    }NR==FNR{
          if($1~/>/){
            sub(/>/,"",$1)
            line[NR]=$1
          }else{
            if(line[NR-1]){
              length_gene[line[NR-1]]=length($0)
              split($0,nt,"")
              for(i=1;i<=length($0);i++){
                hash_codon[line[NR-1]"\t"i]=nt[i]
              }
            }
          }   # length of genes
           
    }NR!=FNR{
    
      if(hash_len[$3-$2]){

             x=$2-"'$ARG_extend'" # location in the extended gene
            if(x<0){
    
              if(x%3+3==3){frame=0}else if(x%3+3==1){frame=1}else if(x%3+3==2){frame=2}
            }else if(x>=0){
              if(x%3==0){frame=0}else if(x%3==1){frame=1}else if(x%3==2){frame=2}
            }       # decide which frame of the 5 end. 
            split($4,a,"|")
            if($3-$2<60){
              for(j=$3+1;j<=$2+60;j++){
                a[2]=a[2]""hash_codon[$1"\t"j]
              }
            }
            
            if(frame==0 && frame_length_hash[$3-$2"\t0"]!=""){
              if(x+frame_length_hash[$3-$2"\t0"]-1>= 3*"'$ARG_exclude_5end'" && x+frame_length_hash[$3-$2"\t0"]-1<=length_gene[$1]-2*"'$ARG_extend'"-1-3*"'$ARG_exclude_3end'"){
                print $0,substr(a[2],frame_length_hash[$3-$2"\t0"],3),substr(a[2],frame_length_hash[$3-$2"\t0"]+3,3),substr(a[2],frame_length_hash[$3-$2"\t0"]+6,3),substr(a[2],frame_length_hash[$3-$2"\t0"]+9,3),x+frame_length_hash[$3-$2"\t0"]-1
              }
            }else if(frame==2 && frame_length_hash[$3-$2"\t2"]!=""){
              if(x+frame_length_hash[$3-$2"\t2"]-1>=3*"'$ARG_exclude_5end'" && x+frame_length_hash[$3-$2"\t2"]-1<=length_gene[$1]-2*"'$ARG_extend'"-1-3*"'$ARG_exclude_3end'"){
                print $0,substr(a[2],frame_length_hash[$3-$2"\t2"],3),substr(a[2],frame_length_hash[$3-$2"\t2"]+3,3),substr(a[2],frame_length_hash[$3-$2"\t2"]+6,3),substr(a[2],frame_length_hash[$3-$2"\t2"]+9,3),x+frame_length_hash[$3-$2"\t2"]-1
              }
            }else if(frame==1 && frame_length_hash[$3-$2"\t1"]!=""){
              if(x+frame_length_hash[$3-$2"\t1"]-1>=3*"'$ARG_exclude_5end'" && x+frame_length_hash[$3-$2"\t1"]-1<=length_gene[$1]-2*"'$ARG_extend'"-1-3*"'$ARG_exclude_3end'"){  
                print $0,substr(a[2],frame_length_hash[$3-$2"\t1"],3),substr(a[2],frame_length_hash[$3-$2"\t1"]+3,3),substr(a[2],frame_length_hash[$3-$2"\t1"]+6,3),substr(a[2],frame_length_hash[$3-$2"\t1"]+9,3),x+frame_length_hash[$3-$2"\t1"]-1  # only useful when 5 end locates at -11 in the frame plot.   
              }
            }
       } # producing codons in the A, +1, +2, +3 sites and the locus of A site in genes. 
    
    }' $ARG_fasta $file_name".bed6" > $ARG_output/$file_name1"_codon_Asite.txt" 


  awk 'BEGIN{FS="\t";OFS="\t"}{
  
      hash[$7]++;hash1[$8]++;hash2[$9]++;hash3[$10]++
  
  }END{
  
    for(i in hash){if(i !~/N/ && i!="TAG" && i!="TAA" && i!="TGA"){print i,hash[i],hash1[i],hash2[i],hash3[i]}} #calculate all the counts of echo codon. 
  
  }'  $ARG_output/$file_name1"_codon_Asite.txt" | awk 'BEGIN{OFS="\t";FS="\t";sumP=0;sumP1=0;sumP2=0;sumP3=0}{
  
    hash[$1]=$0
    sumP=sumP+$2
    sumP1=sumP1+$3
    sumP2=sumP2+$4
    sumP3=sumP3+$5
    hashP[$1]=$2
    hashP1[$1]=$3
    hashP2[$1]=$4
    hashP3[$1]=$5
  }END{
    for(i in hash){
      print i, hashP[i],hashP[i]/sumP,hashP1[i],hashP1[i]/sumP1,hashP2[i],hashP2[i]/sumP2,hashP3[i],hashP3[i]/sumP3  # fraction of each codon. 
    }
  }' | awk 'BEGIN{FS="\t";OFS="\t"}{

     print $0,$3*3/($5+$7+$9),"'${array[1]}'","'${array[6]}'" # normarized to the average of +1,+2,+3 sites. 

   }' | sort -k1,1 >> $ARG_output/all_codon_Asite.txt 

      
  awk 'BEGIN{OFS="\t"}{
    if($1~/>/){
      split($1,a,">")
      hash[NR]=a[2]  # length of the genes. 
    }else if(hash[NR-1]){
      len=length($0)
      for(i="'$ARG_extend'"+1+3*"'$ARG_exclude_5end'";i<=len-"'$ARG_extend'"-3*"'$ARG_exclude_3end'";i=i+3){
        codon=substr($0,i,3)
        gsub(/t/,"T",codon)
        gsub(/g/,"G",codon)
        gsub(/c/,"C",codon)
        gsub(/a/,"A",codon)
        print hash[NR-1],(i-"'$ARG_extend'"-1)/3+1,codon
      }
    }
  }' $ARG_fasta > $ARG_output/$fasta_filename"_codon_excluded.txt" # a file that contains the codons and their locations with both end excluded.  

 awk 'BEGIN{OFS="\t";FS="\t"}{
   hash[$1"\t"$11]++
 }END{for(i in hash){
         print i,hash[i]
       }
 }' $ARG_output/$file_name1"_codon_Asite.txt" > $ARG_output/$file_name1"_codon_Asite_offset.txt"  

  
 rm $ARG_output/$file_name1"_codon_Asite.txt"

done  
  
  echo -e "\n\tNormalizing to wt......"
  awk 'BEGIN{OFS="\t";FS="\t"}{
 
     hash[$1"\t"$11"\t"$12]=hash[$1"\t"$11"\t"$12]""$3","
     hash_nor[$1"\t"$11"\t"$12]=hash_nor[$1"\t"$11"\t"$12]""$10","
     count[$1"\t"$11"\t"$12]++
   }END{
     for(i in hash){
       split(hash[i],a,",")
       split(hash_nor[i],a_nor,",")
       sum=0
       sum_nor=0
       for(j=1;j<length(a);j++){
         sum=sum+a[j];sum_nor=sum_nor+a_nor[j]
       }

        mean=sum/(length(a)-1)
        mean_nor=sum_nor/(length(a_nor)-1)
        square_sum=0
        for(k=1;k<length(a);k++){
          square_sum=square_sum+(a_nor[k]-mean_nor)*(a_nor[k]-mean_nor)
        }
        sd_nor=sqrt(square_sum/(length(a)-1))
        print i,mean,mean_nor,sd_nor   # get the mean of codon occupancy before or after normalizing to average of +1, +2, +3 sites.  
     }
   }' $ARG_output/all_codon_Asite.txt | awk 'BEGIN{OFS="\t";FS="\t"}{hash[$1"\t"$2"\t"$3]=$0;hash2[$1"\t"$3]=$0}END{

     for(i in hash){
       split(i,a,"\t")
       split(hash[i],b,"\t")
       split(hash2[a[1]"\twt"],c,"\t")
       mean_nor_wt=b[5]/c[5]
       sd_nor_wt=b[6]/c[5]
       print hash[i]"\t"mean_nor_wt"\t"sd_nor_wt # mean and sd were normalized to wt. 
     }
   }'  | sort -k1,1 -k2,2  > $ARG_output/all_codon_Asite_nor_wt.txt
 
  if ! [ -x "$(command -v R)" ];then
    echo -e "Warnings: \n\tR is not installed in this server so pdf files can't be made. Please install it or use the files from the results such as \"all_codon_Asite_nor_wt.txt\", \"all_codon_Asite_offset.txt\", to make the plot in Rstudio."

  else 
    echo "Plotting CodonPlot..."           
    R --vanilla --slave --args $ARG_output/all_codon_Asite_nor_wt.txt $ARG_output < codon_plot_log2.r # codon plot including normalized to wt and non-normalized one.
  fi
  
  End_plotting=`date +%s` 
  Time_plotting=$((End_plotting-Start_plotting))
  echo "Took $((Time_plotting/60)) mins for making codon plot"

fi
