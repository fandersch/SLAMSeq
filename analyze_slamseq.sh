#!/bin/bash
#SBATCH -N 1      # nodes requested
#SBATCH -n 1     # tasks requested
#SBATCH -c 1      # cores requested
#SBATCH --qos=medium  # time
#SBATCH --mem=20G  # memory in Mb
#SBATCH --job-name analyze_slamseq
mode=$1
read_length=$2
spikein=$3
rename=$4
umi=$5

if [[ ${mode} = "mm" ]]; then

  echo "mouse"
  UTRbed="/groups/zuber/USERS/florian.andersch/data/genomes/mouse/3UTR/mm10_refseq_ensembl_3UTR_fillup.bed"
  if [[ ${spikein} = "spikein" ]]; then
    genome=""  
  else
    genome="/groups/zuber/USERS/florian.andersch/data/genomes/mouse/fa/GRCm38.fa"
  fi

elif [[ ${mode} = "hs" ]]; then

  echo "human"
  UTRbed="/groups/zuber/USERS/florian.andersch/data/genomes/human/3UTR/hg38_refseq_ensembl_3UTR_fillup.bed"
  if [[ ${spikein} = "spikein" ]]; then
    genome="/groups/zuber/USERS/florian.andersch/data/genomes/GRCh38_dmr6/GRCh38_flybase6.fa"  
  else
    genome="/groups/zuber/USERS/florian.andersch/data/genomes/human/fa/hg38_all.fa"
  fi

else
  echo "please specify mode: mm or hs"
  exit 1
fi

if [[ ! -z "$6" ]]; then
  UTRbed="$6"
fi


if [[ ! -z "$7" ]]; then
  genome="$7"
fi

echo "UTRbed: ${UTRbed}"
echo "UTRgtf: ${UTRgtf}"
echo "Genome: ${UTRgtf}"

wait_for_jobs(){
  sleep 60  # seconds, give time to the schedular to put up the task
  sleeptime=120  # ask every 2 mins, for the first 10 mins
  n=1
  while true; do
    if [ $(squeue | grep florian. | grep -c $1) -eq 0 ]; then
      break
    else
      echo sleep another $((sleeptime / 60)) minutes...
      sleep $sleeptime
    fi
    n=$((n + 1))
    if [ $n -eq 5 ]; then
      sleeptime=300  # if still running after 10 mins, ask every 5 mins
    fi
    if [ $n -eq 10 ]; then
      sleeptime=600  # if still running after 30 mins, ask every 10 mins
    fi
  done
}

source ~/.bashrc
conda activate slamseq_quantseq

#######################
# rename fastq files  #
#######################
if [[ ${rename} = "rename" ]]; then
  srun --time=30 --mem=20G Rscript workflow_scripts/rename_filter_NGS_files.R fastq_${mode} input/input_${mode}.txt
fi

#################################################################################
# check read quality with fastq                                                 #
# @param1: directory of raw reads                                               #
# @param2: output directory                                                     #
#################################################################################

bash workflow_scripts/fastqc.sh fastq_${mode} QC_${mode}

output_dir=""
if [[ ${umi} = "umi" ]]; then

  ###########################################################################################
  # extract umis from read and place into header                                            #
  # @param1: directory of merged reads                                                      #
  # @param2: output directory                                                               #
  ###########################################################################################

  bash workflow_scripts/umitools_extract.sh fastq_${mode}/ fastq_umi_extracted_${mode}/
  wait_for_jobs umi_extr

  ###########################################################
  # check read quality with fastq after umi extraction      #
  # @param1: directory of raw reads                         #
  # @param2: output directory                               #
  ###########################################################

  bash workflow_scripts/fastqc.sh fastq_umi_extracted_${mode}/ QC_umi_${mode}

  output_dir="umi_extracted_"

fi

############################################################################
# trimm adapter and first bases because of low quality,                    #
# @param1: directory of reads (*.fq.gz),                                   #
# @param2: output directory,                                               #
# @param3: quality threshold,                                              #
# @param4: bases to trim at the beginning,                                 #
# @param5: adapter stringency                                              #
############################################################################

bash workflow_scripts/trimming_fastq.sh fastq_${output_dir}${mode}/ trimmed_fastq_${output_dir}${mode}/ 27 4 3
wait_for_jobs trimming

########################3###########################################
# map reads to genome with slam-dunk pipeline (NextGenMap),        #
# @param1: directory of trimmed umi reads (*.fq.gz),               #
# @param2: reference genome fats file,                             #
# @param3: output directory,                                       #
# @param4: 3UTR bed-file                                           #
####################################################################
bash workflow_scripts/slamdunk.sh trimmed_fastq_${output_dir}${mode}/ $genome /$UTRbed $read_length dunk $umi slamdunk
wait_for_jobs slamdunk
srun --time=60 --mem=20G Rscript workflow_scripts/merge_slamseq_count_files_slamdunk.R input/input_${mode}.txt ${mode} slamdunk/counts_collapsed/


#################
# Spike-in      #
#################

if [[ ${spikein} = "spikein" ]]; then
  # Restructure
  rm -r slamdunk/filter_presplit
  rm -r slamdunk/filter_split

  if [[ ${umi} = "umi" ]]; then
    mv slamdunk/dedup slamdunk/filter_presplit
  else
    mv slamdunk/filter slamdunk/filter_presplit
  fi
  
  mv slamdunk/count slamdunk/count_presplit
  mv slamdunk/counts_collapsed slamdunk/counts_collapsed_presplit
  mkdir -p slamdunk/filter
  mkdir -p slamdunk/filter_split
  mkdir -p slamdunk/count
  mkdir -p slamdunk/counts_collapsed

  # De-spike
  workflow_scripts/fileutilities.py T slamdunk/filter_presplit --dir 'bam$' | workflow_scripts/fileutilities.py P --loop sbatch ,--qos=short ,-o /dev/null ,-e /dev/null ,--mem=20G ,-J slamspik workflow_scripts/Split_BAM_Spikein.py ,-a {abs} ,-o slamdunk/filter_split
  wait_for_jobs slamspik

  # Plug subject reads back into slamdunk
  workflow_scripts/fileutilities.py T slamdunk/filter_split --dir subject | perl -e 'while(<>){~s/_subject$//;print}' | workflow_scripts/fileutilities.py P --link slamdunk/filter/

  # Quantify features.
  # Yes I could split apart `slamdunk all` into its constituent steps and insert the de-spiking in between and quantify only once on the correct files,
  # but that would add failure opportunitiess that I don't want to deal with at the moment.
  workflow_scripts/fileutilities.py T slamdunk/filter --dir 'bam$' | workflow_scripts/fileutilities.py P --loop sbatch ,--time=60 ,-J samidx ,-o /dev/null ,-e /dev/null ,--wrap "'samtools index {abs}'"
  wait_for_jobs samidx
  
  #slamdunk count
  workflow_scripts/fileutilities.py T slamdunk/filter --dir 'bam$' | workflow_scripts/fileutilities.py P --loop sbatch ,--time=60 ,-o /dev/null ,-e /dev/null ,-J slamcnt slamdunk count ,-o slamdunk/count ,-s slamdunk/snp ,-r $genome ,-b $UTRbed ,-l $read_length ,-q 27 ,-t 1 {abs}
  wait_for_jobs slamcnt
  workflow_scripts/fileutilities.py T slamdunk/count --dir 'tcount.tsv$' | workflow_scripts/fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,-J alleyoop ,--wrap "'alleyoop collapse -o slamdunk/counts_collapsed {abs}'"
  wait_for_jobs alleyoop
  srun --time=60 --mem=20G Rscript workflow_scripts/merge_slamseq_count_files_slamdunk.R input/input_${mode}.txt ${mode} slamdunk/counts_collapsed/

  # Quantify library partitions. I can use samtools because slamdunk filters out multimappers. (therefore ambiguous should be empty)
  workflow_scripts/fileutilities.py T slamdunk/filter_split --dir 'spikein.bam$' | workflow_scripts/fileutilities.py P --loop printf '"%s\t"' '"{cor}"' '>>' slamdunk/filter_split/spike_counts.txt \&\& samtools view ,-F 4 ,-c {abs} '>>' slamdunk/filter_split/spike_counts.txt
  workflow_scripts/fileutilities.py T slamdunk/filter_split --dir 'ambiguous.bam$' | workflow_scripts/fileutilities.py P --loop printf '"%s\t"' '"{cor}"' '>>' slamdunk/filter_split/ambiguous_counts.txt \&\& samtools view ,-F 4 ,-c {abs} '>>' slamdunk/filter_split/ambiguous_counts.txt
  workflow_scripts/fileutilities.py T slamdunk/filter_split --dir 'subject.bam$' | workflow_scripts/fileutilities.py P --loop printf '"%s\t"' '"{cor}"' '>>' slamdunk/filter_split/subject_counts.txt \&\& samtools view ,-F 4 ,-c {abs} '>>' slamdunk/filter_split/subject_counts.txt
  workflow_scripts/fileutilities.py T slamdunk/filter_split/*counts.txt --appnd -i | perl -e 'while(<>){~s/_\|1//g;~s/_trimmed//;print}' > slamdunk/spike_summary_counts.txt
  wait_for_jobs slamcnt

  # Calculate spike-in scaling factors
  workflow_scripts/scaleFactors_from_spikeSummary.R slamdunk/spike_summary_counts.txt  
  
  # Apply scaling factors
  workflow_scripts/scaleFactors_apply.R -c slamdunk/counts_collapsed/counts.txt -f slamdunk/spike_summary_counts.factors.txt -i 3 -o slamdunk/counts_collapsed/counts
  workflow_scripts/scaleFactors_apply.R -c slamdunk/counts_collapsed/tcCounts.txt -f slamdunk/spike_summary_counts.factors.txt -i 3 -o slamdunk/counts_collapsed/tcCounts

fi

##############
# alleyoop   #
##############
bash workflow_scripts/slamdunk.sh trimmed_fastq_${output_dir}${mode}/ $genome /$UTRbed $read_length alleyoop $umi slamdunk
wait_for_jobs slamdunk

#############
# multiqc
#############
srun --time=60 alleyoop summary -o slamdunk/summary/summary.txt -t slamdunk/count slamdunk/filter/*.bam

bash workflow_scripts/multiqc.sh slamdunk multiqc_final

##########
# DESeq2 #
##########

if [[ ${spikein} = "spikein" ]]; then
  srun --time=60 --mem=20G Rscript workflow_scripts/DESeq2_slamseq.R input/input_${mode}.txt DEA_spikein_scaled_${mode}/ spikein
else
  srun --time=60 --mem=20G Rscript workflow_scripts/DESeq2_slamseq.R input/input_${mode}.txt DEA_${mode}/ nospikein
fi