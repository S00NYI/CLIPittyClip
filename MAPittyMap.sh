terminal_width=$(tput cols)
separator_line=$(printf "%${terminal_width}s" | tr ' ' '*')

## Set Default Values:
# Starting file:
EXP_ID=""
# Type of CLIP experiment for naming:
EXP_TYPE=""
# Keep intermeidate files:
KEEP="no"
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
  echo "MAPittyMap.sh"
  echo "$separator_line"
  echo "This is the mapping-and-onward part of CLIPittyClip."
  echo "Use this sub-pipeline to map processed fasta from CLIPittyClip.sh to a different genome."
  echo "The pipeline utilizes the following programs: "
  echo " - bowtie2"
  echo " - samtools (view, sort, index)"
  echo " - bedtools (bamtobed, genomecov, coverage)"
  echo " - Homer (makeTagDirectory, findPeaks)"
  echo "$separator_line"
}

# Function to display usage information
function show_usage {
  echo "Usage: MapittyMap.sh -i ? -y ? -5 ? -3 ? ... -x ? ..."
  echo "$separator_line"
  echo "- Experiment ID (-i), type (-y) and genome index (-x) are required:"
  echo "  Experiment ID (-i) represents unique identifier for your experiment."
  echo "  Experiment type (-y) represents short descripter for your experiment."
  echo "  The Starting fasta file should have name convention as: id_type.fasta"
  echo "  Example: JL1000_Input.fasta"
  echo "$separator_line"
  echo "- Path and file name to the index for bowtie2 mapping (-x) should follow bowtie2 convention."
  echo "  Path/file name will also be used for the .fa.fai file for bedtools genomecov."
  echo "  Make sure the bowtie2 genome index and bedtools genome files exists in your designated path."
  echo "$separator_line"
  echo "- If using groomed.fasta output from CLIPittyClip.sh:"
  echo "  Copy *groomed.fasta to a separate location (not required, but good for organization.)"
  echo "  In the directory containing the *groomed.fasta, run MapittyMap.sh."
  echo "  For option -i, use the same input as the original CLIPittyClip.sh run."
  echo "  For option -y, add '_collapsed_rmBC_groomed' at the end of the original input."
  echo "  Example: CLIPittyClip.sh -i JL100 -y Input -x path-to-genome-index"
  echo "           MapittyMap.sh -i JL100 -y Input_collapsed_rmBC_groomed -x path-to-diff-genome-index"
  echo "$separator_line"
  echo "- The other options will use default options if not specified."
  echo "$separator_line"
  echo "Options:"
  echo "  -h: print usage information"
  echo "  -i: experiment ID"
  echo "  -y: experiment type (e.g., Input, Enrich, or Fraction)"
  echo "  -k: keep intermediate files                             (default: no)"
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
while getopts "i:y:k:x:m:s:t:p:z:f:h" opt; do
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

## Check if specified fasta file exists
if [ ! -e "${BASE_NAME}.fa" ]; then
  echo "# ${BASE_NAME}.fa file does not exist." | tee -a ${BASE_NAME}_logs.txt
  on_exit
  exit 1
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
## Map reads and sort/index bam file:
echo "\n# Map reads using bowtie2, allowing ${alignMis} mismatches and with ${seedLength}nt seed length:" | tee -a ${BASE_NAME}_logs.txt
input_fasta=${BASE_NAME}.fa 
echo "\n  ${input_fasta}:" | tee -a ${BASE_NAME}_logs.txt
bowtie2 -f -N ${alignMis} -L ${seedLength} -p ${threads} -x ${genome_index} -U ${input_fasta} -S $(basename "$input_fasta" .fa).sam 2>&1 | tee -a ${BASE_NAME}_logs.txt 
samtools view -b -h $(basename "$input_fasta" .fa).sam -o $(basename "$input_fasta" .fa).bam
samtools sort $(basename "$input_fasta" .fa).bam -o $(basename "$input_fasta" .fa).sorted.bam
samtools index -b $(basename "$input_fasta" .fa).sorted.bam $(basename "$input_fasta" .fa).sorted.bam.bai
bedtools bamtobed -i $(basename "$input_fasta" .fa).sorted.bam > $(basename "$input_fasta" .fa).sorted.bed
if [ "$KEEP" == "yes" ]; then
  mv $(basename "$input_fasta" .fa).sam intermediate
  mv $(basename "$input_fasta" .fa).bam intermediate
else
  rm $(basename "$input_fasta" .fa).sam
  rm $(basename "$input_fasta" .fa).bam
fi

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

## Combine all bed files for demultiplexed files:
cat ${BASE_NAME}_mapped_reads/collapsedBED/*.bed > ${BASE_NAME}.bed

mv ${BASE_NAME}.fa ${BASE_NAME}
mv ${BASE_NAME}_mapped_reads ${BASE_NAME}

## Make Tag Directory and call peaks:
echo "\n# Use HOMER to call peaks, with minimum distance of ${PeakMinDist}, peak size of ${PeakSize}, and fragment length of ${FragLength}" | tee -a ${BASE_NAME}_logs.txt
makeTagDirectory ${BASE_NAME}_TagDir/ ${BASE_NAME}.bed -single -format bed 2>&1 | tee -a ${BASE_NAME}_logs.txt 
findPeaks ${BASE_NAME}_TagDir/ -o auto -style factor -L 2 -localSize 10000 -strand separate -minDist ${PeakMinDist} -size ${PeakSize} -fragLength ${FragLength} 2>&1 | tee -a ${BASE_NAME}_logs.txt 

## Process output peaks file to create peaks bed file:
echo "\n# Process peaks.txt to generate peaks BED file." | tee -a ${BASE_NAME}_logs.txt
sed '/^[[:blank:]]*#/d;s/#.*//' ${BASE_NAME}_TagDir/peaks.txt > peaksTemp.bed
awk 'OFS="\t" {print $2, $3, $4, $1, $6, $5}' peaksTemp.bed > peaks.bed
rm peaksTemp.bed 
sort -k 1,1 -k2,2n peaks.bed > peaks_Sorted.bed

mv ${BASE_NAME}.bed ${BASE_NAME}_peaks
mv ${BASE_NAME}_TagDir ${BASE_NAME}_peaks
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
