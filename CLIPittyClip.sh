terminal_width=$(tput cols)
separator_line=$(printf "%${terminal_width}s" | tr ' ' '*')

## Set Default Values:
# Starting file:
EXP_ID=""
# Type of CLIP experiment for naming:
EXP_TYPE=""
# Keep intermeidate files:
KEEP="no"
# Quality filter threshod:
QUALITY=30
# Perform or skip demultiplexing:
DEMUX="no"
# Allowed mismatch for barcode demultiplexing:
barcodeMis=0
# Minimum read length after 5 and 3' trimming/clipping:
minLength=16
# Default 5'end nucleotide to clip:
adapter5=10
# Default 3'end adapter sequence (L32):
adapter3="GTGTCAGTCACTTCCAGCGG"
# Path to genome index:
genome_index=""
# Allowed mismatch
alignMis=0
# Seed length
seedLength=15
# Number of threads
threads=1
# Peak calling minimum distance for HOMER 
PeakMinDist=50
# Peak size for HOMER
PeakSize=20
# Fragment size for HOMER
FragLength=25

# Function to display information
function show_info {
  echo "$separator_line"
  echo "CLIPittyClip: Single-line CLIP data analysis pipeline"
  echo "$separator_line"
  echo "Version 2.0.0"
  echo "Author: Soon Yi"
  echo "Last updated: 2024-01-03"
  echo "$separator_line"
  echo "CLIPittyClip.sh"
  echo "$separator_line"
  echo "This is a standard CLIP analysis pipeline, from fastq.gz to peaks."
  echo "If demultiplexing is needed, set -d option to 'yes' and provide the proper barcode file."
  echo "The pipeline utilizes the following programs: "
  echo " - fastx_toolkit (collapser, barcode_splitter, trimmer, clipper)"
  echo " - ctk (stripBarcode, tag2collapse)"
  echo " - bowtie2"
  echo " - samtools (view, sort, index)"
  echo " - bedtools (bamtobed, genomecov, coverage)"
  echo " - Homer (makeTagDirectory, findPeaks)"
  echo "$separator_line"
}

# Function to display usage information
function show_usage {
  echo "Usage: CLIPittyClip.sh -i ? -y ? -5 ? -3 ? ... -x ? ..."
  echo "$separator_line"
  echo "- Experiment ID (-i), type (-y) and genome index (-x) are required:"
  echo "  Experiment ID (-i) represents a unique identifier for your experiment."
  echo "  Experiment type (-y) represents a short descripter for your experiment."
  echo "  Starting fastq.gz file should have a name convention as: id_type.fastq.gz"
  echo "  Example: JL1000_Input.fastq.gz"
  echo "$separator_line"
  echo "- If demultiplexing is required, make sure to set -d option to 'yes'."
  echo "  The Barcode file should have the same name convention as fastq.gz file, followed by '_BC.txt'."
  echo "  Example: JL1000_Input_BC.txt"
  echo "  Check our GitHub page for sample barcode file."
  echo "$separator_line"
  echo "- Path and file name to the index for bowtie2 mapping (-x) should follow bowtie2 convention."
  echo "  Path/file name will also be used for the .fa.fai file for bedtools genomecov."
  echo "  Make sure the bowtie2 genome index and bedtools genome files exists in your designated path."
  echo "$separator_line"
  echo "- The other options will use default options if not specified."
  echo "$separator_line"
  echo "Options:"
  echo "  -h: print usage information"
  echo "  -i: experiment ID"
  echo "  -y: experiment type (e.g., Input, Enrich, or Fraction)"
  echo "  -k: keep intermediate files                             (default: no)"
  echo "  -q: quality score threshold                             (default: 30)"
  echo "  -d: perform demultiplexing                              (default: no)"
  echo "  -b: number of mismatch allowed in barcode for demux     (default:  0)"
  echo "  -5: number of nucleotide to clip from the 5'-end        (default: 10)"
  echo "  -3: sequence of the 3'-end adapter                      (default: L32 sequence)"
  echo "  -l: minimum read length after 5/3 end trimming/clipping (default: 16)"
  echo "  -x: path to the genome index for bowtie2 and bedtools"
  echo "  -m: number of mismatch allowed for mapping              (default:  0)"
  echo "  -s: seed length for mapping                             (default: 15)"
  echo "  -t: number of threads to be utilized by bowtie/samtools (default:  1)"
  echo "  -p: minimum distance between peaks for homer            (default: 50)"
  echo "  -z: size of peaks for homer                             (default: 20)"
  echo "  -f: fragment length for homer                           (default: 25)"
  echo "$separator_line"
}

function on_exit {
  read -p "Press enter to exit."
}

# Check if arguments are provided:
if [[ $# -eq 0 ]]; then
  show_info
  echo "Error: No option arguments have been passed. use '-h' option for usage information."
  on_exit
  exit 1
fi

# Parse command-line options
while getopts "i:y:k:q:d:b:5:3:l:x:m:s:t:p:z:f:h" opt; do
  case $opt in
    i)
      EXP_ID="$OPTARG"
      ;;
    y)
      EXP_TYPE="$OPTARG"
      ;;
    k)
      KEEP="$OPTARG"
    ;;
    q)
      QUALITY="$OPTARG"
      ;;
    d)
      DEMUX="$OPTARG"
      ;;
    b)
      barcodeMis="$OPTARG"
      ;;
    5)
      adapter5="$OPTARG"
      ;;
    3)
      adapter3="$OPTARG"
      ;;
    l)
      minLength="$OPTARG"
      ;;
    x)
      genome_index="$OPTARG"
      ;;
    m)
      alignMis="$OPTARG"
      ;;
    s)
      seedLength="$OPTARG"
      ;;
    t)
      threads="$OPTARG"
      ;;
    p)
      PeakMinDist="$OPTARG"
      ;;
    z)
      PeakSize="$OPTARG"
      ;;
    f)
      FragLength="$OPTARG"
      ;;
    h)
      show_info
      show_usage
      on_exit
      exit 0
      ;;
    \?)
      echo "$separator_line"
      echo "Error: Invalid option -$OPTARG" >&2
      echo "$separator_line"
      show_usage
      on_exit
      exit 1
      ;;
    :)
      echo "$separator_line"
      echo "Error: Option -$OPTARG requires an argument." >&2
      echo "$separator_line"
      show_usage
      on_exit
      exit 1
      ;;
  esac
done

# Check if all required options are provided
if [[ -z "$EXP_ID" || -z "$EXP_TYPE" || -z "$genome_index" ]]; then
  echo "$separator_line"
  echo "Error: Options -i, -y, and -x are required."
  echo "$separator_line"
  show_usage
  on_exit
  exit 1
fi

## Start of pipeline:
#####################################################################################
BASE_NAME="${EXP_ID}_${EXP_TYPE}"

show_info | tee -a ${BASE_NAME}_logs.txt
echo "# Analysis started: $(date +'%Y/%m/%d %H:%M')" | tee -a ${BASE_NAME}_logs.txt
echo "# User Input: '$0 $@'\n" | tee -a ${BASE_NAME}_logs.txt

if [ "$KEEP" == "yes" ]; then
  echo "# All intermediate files will be kept. Note this may use up a lot of disk space." | tee -a ${BASE_NAME}_logs.txt
else
  echo "# All intermediate files will be removed as the pipeline proceeds." | tee -a ${BASE_NAME}_logs.txt
fi

## Check if specified fastq file exists:
if [ ! -e "${BASE_NAME}.fastq.gz" ]; then
  echo "# ${BASE_NAME}.fastq.gz file does not exist." | tee -a ${BASE_NAME}_logs.txt
  on_exit
  exit 1
fi

## Check if demultiplexing option is on:
if [ "$DEMUX" == "yes" ]; then
  # Check if the BC file exists in the current directory.
  if [ ! -e "${BASE_NAME}_BC.txt" ]; then
  # If BC file does not exists, exit.
    echo "# Barcode information not found! Make sure to provide proper barcode file." | tee -a ${BASE_NAME}_logs.txt
    echo "# Barcode file should have the same name convention as fastq.gz file, followed by '_BC.txt'.\n" | tee -a ${BASE_NAME}_logs.txt
    on_exit
    exit 1
  else
    echo "# Barcode information found! Demultiplexing will be performed based on the barcode information.\n" | tee -a ${BASE_NAME}_logs.txt
  fi
else
  echo "# Demultiplexing will be skipped. Option -d will be ignored.\n" | tee -a ${BASE_NAME}_logs.txt
fi

## Check if bedtools genome index file exists.
if [ ! -e "${genome_index}.fa.fai" ]; then
  # If bedtools index file does not exists, exit.
  echo "# Genome index for bedtools not found! Make sure to provide proper genome file for bedtools.\n" | tee -a ${BASE_NAME}_logs.txt
  on_exit
  exit 1
else
  echo "# Genome index for bedtools found!\n" | tee -a ${BASE_NAME}_logs.txt
fi

## Check if bowtie2 genome index file exists.
bt2_index=("1.bt2" "2.bt2" "3.bt2" "4.bt2" "rev.1.bt2" "rev.2.bt2")
for index in "${bt2_index[@]}"; do
  file_to_check="${genome_index}.${index}"
  if [ ! -e "${file_to_check}" ]; then
    # If bowtie2 index files doesn't exists, exit.
    echo "# ${index} index file not found! Make sure all bowtie2 genome index files are in the provided directory.\n" | tee -a ${BASE_NAME}_logs.txt
    on_exit
    exit 1
  else
    echo "# ${index} index file found!" | tee -a ${BASE_NAME}_logs.txt
  fi
done

## Make directories:
#####################################################################################
echo "\n# Create necessary directories: " | tee -a ${BASE_NAME}_logs.txt

mkdir ${BASE_NAME}
mkdir ${BASE_NAME}_groomed_reads
mkdir ${BASE_NAME}_mapped_reads
mkdir BAM
mkdir BED
mkdir collapsedBED
mkdir BEDGRAPH
mkdir peakCoverage
mkdir ${BASE_NAME}_peaks
if [ "$KEEP" == "yes" ]; then
  mkdir intermediate
fi


## Actual start of the pipeline: 
#####################################################################################
## Unzip fastq.gz:
echo "\n# Unzip fasta.gz file." | tee -a ${BASE_NAME}_logs.txt
gunzip -k ${BASE_NAME}.fastq.gz

## Filter reads based on quality:
echo "\n# Quality fiilter reads using fastx fastq_quality_filter:" | tee -a ${BASE_NAME}_logs.txt
echo "\n  ${BASE_NAME}.fastq --> ${BASE_NAME}_filtered.fastq" | tee -a ${BASE_NAME}_logs.txt
fastq_quality_filter -v -q ${QUALITY} -i ${BASE_NAME}.fastq -o ${BASE_NAME}_filtered.fastq | tee -a ${BASE_NAME}_logs.txt

## Collapse identical reads:
echo "\n# Collapse identical reads using fastx_collapser:" | tee -a ${BASE_NAME}_logs.txt
echo "\n  ${BASE_NAME}_filtered.fastq --> ${BASE_NAME}_collapsed.fa" | tee -a ${BASE_NAME}_logs.txt
fastx_collapser -v -i ${BASE_NAME}_filtered.fastq -o ${BASE_NAME}_collapsed.fa | tee -a ${BASE_NAME}_logs.txt
if [ "$KEEP" == "yes" ]; then
  mv *.fastq intermediate
else
  rm *.fastq
fi


## Strip UMI:
echo "\n# Strip UMIs from the reads using CTK stripBarcode.pl:" | tee -a ${BASE_NAME}_logs.txt
echo "\n  ${BASE_NAME}_collapsed.fa --> ${BASE_NAME}_collapsed_rmBC.fa" | tee -a ${BASE_NAME}_logs.txt
stripBarcode.pl -len 7 -v ${BASE_NAME}_collapsed.fa ${BASE_NAME}_collapsed_rmBC.fa
if [ "$KEEP" == "yes" ]; then
  mv *_collapsed.fa intermediate
else
  rm *_collapsed.fa
fi


## Perform demultiplexing if DEMUX option is yes:
if [ "$DEMUX" == "yes" ]; then
  ## Demultiplex Samples:
  echo "\n# Demultiplex samples with barcodes using fastx_barcode_splitter, allowing ${barcodeMis} mismatches:" | tee -a ${BASE_NAME}_logs.txt
  echo "\n  ${BASE_NAME}_collapsed_rmBC.fa --> demultiplexed .fa files.\n  Unmatched.fa file will be removed at the end." | tee -a ${BASE_NAME}_logs.txt
  cat ${BASE_NAME}_collapsed_rmBC.fa | fastx_barcode_splitter.pl --bcfile ${BASE_NAME}_BC.txt --bol --mismatches ${barcodeMis} --prefix "${BASE_NAME}_" --suffix ".fa" | tee -a ${BASE_NAME}_logs.txt
  if [ "$KEEP" == "yes" ]; then
    mv *_collapsed_rmBC.fa intermediate
    mv *unmatched.fa intermediate
  else
    rm *_collapsed_rmBC.fa
    rm *unmatched.fa
  fi
fi

## Barcode, Spacer and L32 trimming.
echo "\n# Trim ${adapter5}nt from 5'-end using fastx_trimmer, followed by 3'-end adapter ${adapter3} removal using fastx_clipper, and keeping reads with minimum length of ${minLength}:" | tee -a ${BASE_NAME}_logs.txt
for input_file in ${BASE_NAME}*.fa; do
	echo "\n  ${input_file}:" | tee -a ${BASE_NAME}_logs.txt
  output_file="$(basename "${input_file}" .fa)_groomed.fa"
  fastx_trimmer -f ${adapter5} -i "${input_file}" | fastx_clipper -l ${minLength} -a ${adapter3} -v -o ${output_file} | tee -a ${BASE_NAME}_logs.txt
  echo "  Remaining Reads: $(awk "BEGIN { printf \"%.2f\", $(wc -l < $(basename "${input_file}" .fa)_groomed.fa) / $(wc -l < "$input_file") * 100 }")"% | tee -a ${BASE_NAME}_logs.txt
  if [ "$KEEP" == "yes" ]; then
    mv ${input_file} intermediate
  else
    rm ${input_file}
  fi
done

## Map reads and sort/index bam file:
echo "\n# Map reads using bowtie2 (${alignMis} mismatches with ${seedLength}nt seed length) and create sorted/indexed bam/bed files:" | tee -a ${BASE_NAME}_logs.txt
for input_fasta in *_groomed.fa; do
  echo "\n  ${input_fasta}:" | tee -a ${BASE_NAME}_logs.txt
  bowtie2 -f -N ${alignMis} -L ${seedLength} -p ${threads} -x ${genome_index} -U ${input_fasta} -S $(basename "$input_fasta" _groomed.fa).sam 2>&1 | tee -a ${BASE_NAME}_logs.txt 
  if [ "$threads" == "1" ]; then
    samtools view -b -h $(basename "$input_fasta" _groomed.fa).sam -o $(basename "$input_fasta" _groomed.fa).bam
    samtools sort $(basename "$input_fasta" _groomed.fa).bam -o $(basename "$input_fasta" _groomed.fa).sorted.bam
    samtools index -b $(basename "$input_fasta" _groomed.fa).sorted.bam $(basename "$input_fasta" _groomed.fa).sorted.bam.bai
  else
    samtools view -b -@ ${threads} -h $(basename "$input_fasta" _groomed.fa).sam -o $(basename "$input_fasta" _groomed.fa).bam
    samtools sort -@ ${threads} $(basename "$input_fasta" _groomed.fa).bam -o $(basename "$input_fasta" _groomed.fa).sorted.bam
    samtools index -@ ${threads} -b $(basename "$input_fasta" _groomed.fa).sorted.bam $(basename "$input_fasta" _groomed.fa).sorted.bam.bai
  fi
  bedtools bamtobed -i $(basename "$input_fasta" _groomed.fa).sorted.bam > $(basename "$input_fasta" _groomed.fa).sorted.bed
  mv ${input_fasta} ${BASE_NAME}_groomed_reads 
  if [ "$KEEP" == "yes" ]; then
    mv $(basename "$input_fasta" _groomed.fa).sam intermediate
    mv $(basename "$input_fasta" _groomed.fa).bam intermediate
  else
    rm $(basename "$input_fasta" _groomed.fa).sam
    rm $(basename "$input_fasta" _groomed.fa).bam
  fi
done

## Move files to corresponding directories:
mv *sorted.bam BAM
mv *sorted.bam.bai BAM
mv BAM ${BASE_NAME}_mapped_reads

## Collapse PCR dups and make bedgraph:
echo "\n# Collapse PCR duplicates from the mapped reads using CTK tag2collapse.pl." | tee -a ${BASE_NAME}_logs.txt
echo "# From collapsed bed file, use bedtools genomecov to make bedgraph files." | tee -a ${BASE_NAME}_logs.txt
for input_bed in *.bed; do
  tag2collapse.pl --keep-tag-name --keep-max-score --random-barcode -EM -1 --seq-error-model em-local -weight --weight-in-name  ${input_bed} $(basename "$input_bed" .bed).collapsed.bed | tee -a ${BASE_NAME}_logs.txt
  mv ${input_bed} BED
  bedtools genomecov -i $(basename "$input_bed" .bed).collapsed.bed -strand + -g ${genome_index}.fa.fai -bg > $(basename "$input_bed" .bed).collapsed.pos.bedgraph
  bedtools genomecov -i $(basename "$input_bed" .bed).collapsed.bed -strand - -g ${genome_index}.fa.fai -bg > $(basename "$input_bed" .bed).collapsed.rev.bedgraph
  mv $(basename "$input_bed" .bed).collapsed.bed collapsedBED
  mv $(basename "$input_bed" .bed).collapsed.pos.bedgraph BEDGRAPH
  mv $(basename "$input_bed" .bed).collapsed.rev.bedgraph BEDGRAPH
done

## Move files to corresponding directories:
mv BED ${BASE_NAME}_mapped_reads
mv collapsedBED ${BASE_NAME}_mapped_reads
mv BEDGRAPH ${BASE_NAME}_mapped_reads

if [ "$DEMUX" == "yes" ]; then
  HOMER_INPUT_NAME=${BASE_NAME}_combined
else
  HOMER_INPUT_NAME=${BASE_NAME}
fi

## Combine all bed files for demultiplexed files:
if [ "$DEMUX" == "yes" ]; then
  echo "\n# Combine collapsed BED files to make global BED file for the experiment." | tee -a ${BASE_NAME}_logs.txt
  cat ${BASE_NAME}_mapped_reads/collapsedBED/*.bed > ${HOMER_INPUT_NAME}.bed
else
  cat ${BASE_NAME}_mapped_reads/collapsedBED/*.bed > ${BASE_NAME}.bed
fi

mv ${BASE_NAME}_groomed_reads ${BASE_NAME}
mv ${BASE_NAME}_mapped_reads ${BASE_NAME}

## Make Tag Directory and call peaks:
echo "\n# Use HOMER to call peaks, with minimum distance of ${PeakMinDist}, peak size of ${PeakSize}, and fragment length of ${FragLength}" | tee -a ${BASE_NAME}_logs.txt
makeTagDirectory ${BASE_NAME}_TagDir/ ${BASE_NAME}.bed -single -format bed 2>&1 | tee -a ${BASE_NAME}_logs.txt 
findPeaks ${BASE_NAME}_TagDir/ -o auto -style factor -L 2 -localSize 10000 -strand separate -minDist ${PeakMinDist} -size ${PeakSize} -fragLength ${FragLength} 2>&1 | tee -a ${BASE_NAME}_logs.txt 

## Process output peaks file to create peaks bed file:
echo "\n# Process peaks.txt to generate peaks BED file." | tee -a ${BASE_NAME}_logs.txt
sed '/^[[:blank:]]*#/d;s/#.*//' ${HOMER_INPUT_NAME}_TagDir/peaks.txt > peaksTemp.bed
awk 'OFS="\t" {print $2, $3, $4, $1, $6, $5}' peaksTemp.bed > peaks.bed
rm peaksTemp.bed 
sort -k 1,1 -k2,2n peaks.bed > peaks_Sorted.bed

mv ${HOMER_INPUT_NAME}.bed ${BASE_NAME}_peaks
mv ${HOMER_INPUT_NAME}_TagDir ${BASE_NAME}_peaks
mv peaks.bed ${BASE_NAME}_peaks

## 
echo "\n# Calculate coverages for the peaks for each sample using bedtools coverage." | tee -a ${BASE_NAME}_logs.txt
for id in ${BASE_NAME}/${BASE_NAME}_mapped_reads/collapsedBED/*.bed; do
  bedtools coverage -s -a peaks_Sorted.bed -b "${id}" > "coverage_$(basename "$id")"
done

## move columns to peaks_sorted.bed and make new file.
echo "\n# Make coverage table for the peaks spanning all samples." | tee -a ${BASE_NAME}_logs.txt
cp peaks_Sorted.bed ${BASE_NAME}_peakCoverage.txt
echo "chr\tstart\tend\tname\tscore\tstrand" > colnames.txt
cat ${BASE_NAME}_peakCoverage.txt >> colnames.txt
mv colnames.txt ${BASE_NAME}_peakCoverage.txt

for cov_file in coverage_*; do
  echo $(echo $(basename "${cov_file}" .sorted.collapsed.bed) | cut -d '_' -f 2-) > temp1.txt
  awk 'FNR>0 {print $7}' ${cov_file} > temp2.txt && cat temp1.txt temp2.txt > temp3.txt
  paste ${BASE_NAME}_peakCoverage.txt temp3.txt > temp4.txt && mv temp4.txt ${BASE_NAME}_peakCoverage.txt
  rm temp*.txt
done

echo "\n# Final output: ${BASE_NAME}_peakCoverage.txt" | tee -a ${BASE_NAME}_logs.txt

mv peaks_Sorted.bed ${BASE_NAME}_peaks
mv coverage_* peakCoverage
mv peakCoverage ${BASE_NAME}_peaks
mv ${BASE_NAME}_peakCoverage.txt ${BASE_NAME}_peaks
mv ${BASE_NAME}_peaks ${BASE_NAME}
if [ "$KEEP" == "yes" ]; then
  mv intermediate ${BASE_NAME}
fi


echo "\n# Analysis Finished: $(date +'%Y/%m/%d %H:%M')" | tee -a ${BASE_NAME}_logs.txt
echo "\n# All outputs are in ${BASE_NAME} folder." | tee -a ${BASE_NAME}_logs.txt
if [ "$KEEP" == "yes" ]; then
  echo "\n# All intermediate files are in ${BASE_NAME}/intermediate folder." | tee -a ${BASE_NAME}_logs.txt
fi
mv ${BASE_NAME}_logs.txt ${BASE_NAME}
