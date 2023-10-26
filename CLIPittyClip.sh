## Set Default Values:
# Starting file:
EXP_ID=""
# Type of CLIP experiment for naming:
EXP_TYPE=""
# Quality filter threshod:
QUALITY=30
# Perform or skip demultiplexing:
DEMUX="no"
# Allowed mismatch for barcode demultiplexing:
barcodeMis=0
# Minimum read length after 5 and 3' trimming/clipping:
minLength=16
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
  echo "-----------------------------------------------------------------------------------------------"
  echo "CLIPittyClip: Single-line CLIP data analysis pipeline"
  echo "-----------------------------------------------------------------------------------------------"
  echo "Version 1.0.0"
  echo "Author: Soon Yi"
  echo "Last update: 2023-10-25"
  echo "-----------------------------------------------------------------------------------------------"
  echo "This is a standard CLIP analysis pipeline, from fastq.gz to peaks."
  echo "If input fastq requires demultiplexing, use -d option and provide proper barcode file."
  echo "The pipeline utilizes following programs: "
  echo " - fastx_toolkit (collapser, barcode_splitter, trimmer, clipper)"
  echo " - ctk (stripBarcode, tag2collapse)"
  echo " - bowtie2"
  echo " - samtools (view, sort, index)"
  echo " - bedtools (bamtobed, genomecov, coverage)"
  echo " - Homer (makeTagDirectory, findPeaks)"
  echo "Make sure conda environment has all the necessary programs installed."
  echo "-----------------------------------------------------------------------------------------------"
}

# Function to display usage information
function show_usage {
  echo "-----------------------------------------------------------------------------------------------"
  echo "Usage: CLIPittyClip.sh -i ? -y ? ... -x ? ..."
  echo "- Experiment ID (-i), type (-y) and genome index (-x) are required:"
  echo "  Experiment ID (-i) represents unique identifier for your experiment."
  echo "  Experiment type (-y) represents short descripter for your experiment."
  echo "  Starting fastq.gz file should have name convention as: id_type.fastq.gz"
  echo "  Example: JL1000_Input.fastq.gz"
  echo "-----------------------------------------------------------------------------------------------"
  echo "- If demultiplexing is required, make sure to set -d option to yes."
  echo "  Barcode file should have the same name convention as fastq.gz file, followed by '_BC.txt'."
  echo "  Example: JL1000_Input_BC.txt"
  echo "- Check our GitHub page for sample barcode file."
  echo "-----------------------------------------------------------------------------------------------"
  echo "- Path and file name to genome index for bowtie2 mapping (-x) should follow bowtie2 convention."
  echo "  Path/file name will also be used for .fa.fai file for bedtools genomecov."
  echo "  Makre sure bowtie2 genome index and bedtools genome files exists in your designated path."
  echo "-----------------------------------------------------------------------------------------------"
  echo "- The other options will use default options if not specified."
  echo "-----------------------------------------------------------------------------------------------"
  echo "Options:"
  echo "  -h: print usage information"
  echo "  -i: experiment ID"
  echo "  -y: experiment type (e.g., Input, Enrich, or Fraction)"
  echo "  -q: quality score threshold                             (default: 30)"
  echo "  -d: barcode file for demultiplexing.                    (default: no)"
  echo "  -b: allowed barcode mismatch for demultiplexing         (default:  0)"
  echo "  -l: minimum read length after 5/3 end trimming/clipping (default: 16)"
  echo "  -x: path to genome index for bowtie2 mapping"
  echo "  -m: allowed mismatch for mapping                        (default:  0)"
  echo "  -s: seed length for mapping                             (default: 15)"
  echo "  -t: number of threads to utilize for mapping            (default:  1)"
  echo "  -p: minimum distance between peaks for homer            (default: 50)"
  echo "  -z: size of peaks for homer                             (default: 20)"
  echo "  -f: fragment length for homer                           (default: 25)"
  echo "-----------------------------------------------------------------------------------------------"
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
while getopts "i:y:q:d:b:l:x:m:s:t:p:z:f:h" opt; do
  case $opt in
    i)
      EXP_ID="$OPTARG"
      ;;
    y)
      EXP_TYPE="$OPTARG"
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
      echo "Error: Invalid option -$OPTARG" >&2
      show_usage
      on_exit
      exit 1
      ;;
    :)
      echo "Error: Option -$OPTARG requires an argument." >&2
      show_usage
      on_exit
      exit 1
      ;;
  esac
done

# Check if all required options are provided
if [[ -z "$EXP_ID" || -z "$EXP_TYPE" || -z "$genome_index" ]]; then
  echo "Error: Options -i, -y, and -x are required."
  show_usage
  on_exit
  exit 1
fi

# Start pipeline
BASE_NAME="${EXP_ID}_${EXP_TYPE}"

show_info | tee -a ${BASE_NAME}_logs.txt
echo "# Analysis started: $(date +'%Y:%m:%d:%H:%M')\n" | tee -a ${BASE_NAME}_logs.txt

## Check if specified fastq file exists
if [ ! -e "${BASE_NAME}.fastq.gz" ]; then
  echo "# ${BASE_NAME}.fastq.gz file does not exist." | tee -a ${BASE_NAME}_logs.txt
  on_exit
  exit 1
fi

## Check if demultiplexing option is on.
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
  echo "# Demultiplexing will be skipped.\n" | tee -a ${BASE_NAME}_logs.txt
fi

## Check if bedtools genome index file exists.
if [ ! -e "${genome_index}.fa.fai" ]; then
# If BC file does not exists, exit.
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
    echo "# ${index} index file not found! Make sure all bowtie2 genome index files are in the provided directory.\n" | tee -a ${BASE_NAME}_logs.txt
    on_exit
    exit 1
  else
    echo "# ${index} index file found!" | tee -a ${BASE_NAME}_logs.txt
  fi
done

## Make directories:
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
rm *.fastq

## Strip UMI:
echo "\n# Strip UMIs from the reads using CTK stripBarcode.pl:" | tee -a ${BASE_NAME}_logs.txt
echo "\n  ${BASE_NAME}_collapsed.fa --> ${BASE_NAME}_collapsed_rmBC.fa" | tee -a ${BASE_NAME}_logs.txt
stripBarcode.pl -len 7 -v ${BASE_NAME}_collapsed.fa ${BASE_NAME}_collapsed_rmBC.fa
rm *_collapsed.fa

## Perform demultiplexing if Demux option is yes:
if [ "$DEMUX" == "yes" ]; then
  ## Demultiplex Samples:
  echo "\n# Demultiplex samples with barcodes using fastx_barcode_splitter, allowing ${barcodeMis} mismatches:" | tee -a ${BASE_NAME}_logs.txt
  echo "\n  ${BASE_NAME}_collapsed_rmBC.fa --> demultiplexed .fa files.\n  Unmatched.fa file will be removed at the end." | tee -a ${BASE_NAME}_logs.txt
  cat ${BASE_NAME}_collapsed_rmBC.fa | fastx_barcode_splitter.pl --bcfile ${BASE_NAME}_BC.txt --bol --mismatches ${barcodeMis} --prefix "${BASE_NAME}_" --suffix ".fa" | tee -a ${BASE_NAME}_logs.txt
  rm *_collapsed_rmBC.fa
  rm *unmatched.fa
fi

## Barcode, Spacer and L32 trimming.
echo "\n# Trim barcode and spacer (GGG) using fastx_trimmer, followed by L32 removal using fastx_clipper, and keeping reads with minimum length of ${minLength}:" | tee -a ${BASE_NAME}_logs.txt
for input_file in ${BASE_NAME}*.fa; do
	echo "\n  ${input_file}:" | tee -a ${BASE_NAME}_logs.txt
    output_file="$(basename "${input_file}" .fa)_groomed.fa"
    fastx_trimmer -f 10 -i "${input_file}" | fastx_clipper -l ${minLength} -a GTGTCAGTCACTTCCAGCGG -v -o ${output_file} | tee -a ${BASE_NAME}_logs.txt
    echo "  Remaining Reads: $(awk "BEGIN { printf \"%.2f\", $(wc -l < $(basename "${input_file}" .fa)_groomed.fa) / $(wc -l < "$input_file") * 100 }")"% | tee -a ${BASE_NAME}_logs.txt
    rm ${input_file}
done

## Map reads and sort/index bam file:
echo "\n# Map reads using bowtie2, allowing ${alignMis} mismatches and with ${seedLength}nt seed length:" | tee -a ${BASE_NAME}_logs.txt
for input_fasta in *_groomed.fa; do
    echo "\n  ${input_fasta}:" | tee -a ${BASE_NAME}_logs.txt
    bowtie2 -f -N ${alignMis} -L ${seedLength} -p ${threads} -x ${genome_index} -U ${input_fasta} -S $(basename "$input_fasta" _groomed.fa).sam 2>&1 | tee -a ${BASE_NAME}_logs.txt 
    samtools view -b -h $(basename "$input_fasta" _groomed.fa).sam -o $(basename "$input_fasta" _groomed.fa).bam
    samtools sort $(basename "$input_fasta" _groomed.fa).bam -o $(basename "$input_fasta" _groomed.fa).sorted.bam
    samtools index -b $(basename "$input_fasta" _groomed.fa).sorted.bam $(basename "$input_fasta" _groomed.fa).sorted.bam.bai
    bedtools bamtobed -i $(basename "$input_fasta" _groomed.fa).sorted.bam > $(basename "$input_fasta" _groomed.fa).sorted.bed
    mv ${input_fasta} ${BASE_NAME}_groomed_reads
    rm $(basename "$input_fasta" _groomed.fa).sam
    rm $(basename "$input_fasta" _groomed.fa).bam
done

## Make directory for mapped reads:
mv *sorted.bam BAM
mv *sorted.bam.bai BAM
mv BAM ${BASE_NAME}_mapped_reads

## Collapse PCR dups
echo "\n# Collapse PCR duplicates from the mapped reads using CTK tag2collapse.pl." | tee -a ${BASE_NAME}_logs.txt
echo "\n# From collapsed bed file, use bedtools genomecov to make bedgraph files." | tee -a ${BASE_NAME}_logs.txt
for input_bed in *.bed; do
    tag2collapse.pl --keep-tag-name --keep-max-score --random-barcode -EM 30 --seq-error-model em-local -weight --weight-in-name  ${input_bed} $(basename "$input_bed" .bed).collapsed.bed
    mv ${input_bed} BED
    bedtools genomecov -i $(basename "$input_bed" .bed).collapsed.bed -strand + -g ${genome_index}.fa.fai -bg > $(basename "$input_bed" .bed).collapsed.pos.bedgraph
    bedtools genomecov -i $(basename "$input_bed" .bed).collapsed.bed -strand - -g ${genome_index}.fa.fai -bg > $(basename "$input_bed" .bed).collapsed.rev.bedgraph
    mv $(basename "$input_bed" .bed).collapsed.bed collapsedBED
    mv $(basename "$input_bed" .bed).collapsed.pos.bedgraph BEDGRAPH
    mv $(basename "$input_bed" .bed).collapsed.rev.bedgraph BEDGRAPH
done

mv BED ${BASE_NAME}_mapped_reads
mv collapsedBED ${BASE_NAME}_mapped_reads
mv BEDGRAPH ${BASE_NAME}_mapped_reads

## Combine all bed files:
echo "\n# Combine collapsed BED files to make global BED file for the experiment." | tee -a ${BASE_NAME}_logs.txt
cat ${BASE_NAME}_mapped_reads/collapsedBED/*.bed > ${BASE_NAME}_combined.bed 

mv ${BASE_NAME}_groomed_reads ${BASE_NAME}
mv ${BASE_NAME}_mapped_reads ${BASE_NAME}

## Make Tag Directory and call peaks:
echo "\n# Use HOMER to call peaks, with minimum distance of ${PeakMinDist}, peak size of ${PeakSize}, and fragment length of ${FragLength}" | tee -a ${BASE_NAME}_logs.txt
makeTagDirectory ${BASE_NAME}_combined_TagDir/ ${BASE_NAME}_combined.bed -single -format bed 2>&1 | tee -a ${BASE_NAME}_logs.txt 
findPeaks ${BASE_NAME}_combined_TagDir/ -o auto -style factor -L 2 -localSize 10000 -strand separate -minDist ${PeakMinDist} -size ${PeakSize} -fragLength ${FragLength} 2>&1 | tee -a ${BASE_NAME}_logs.txt 

## Process output peaks file to create peaks bed file:
echo "\n# Process peaks.txt to generate peaks BED file." | tee -a ${BASE_NAME}_logs.txt
sed '/^[[:blank:]]*#/d;s/#.*//' ${BASE_NAME}_combined_TagDir/peaks.txt > peaksTemp.bed
awk 'OFS="\t" {print $2, $3, $4, $1, $6, $5}' peaksTemp.bed > peaks.bed
rm peaksTemp.bed 
sort -k 1,1 -k2,2n peaks.bed > peaks_Sorted.bed
grep -E '^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)\s' < peaks_Sorted.bed > peaks_Sorted_Filtered.bed

mv ${BASE_NAME}_combined.bed ${BASE_NAME}_peaks
mv ${BASE_NAME}_combined_TagDir ${BASE_NAME}_peaks
mv peaks.bed ${BASE_NAME}_peaks
mv peaks_Sorted.bed ${BASE_NAME}_peaks

## 
echo "\n# Calculate coverages for the peaks for each sample using bedtools coverage." | tee -a ${BASE_NAME}_logs.txt
for id in ${BASE_NAME}/${BASE_NAME}_mapped_reads/collapsedBED/*.bed; do
    bedtools coverage -s -a peaks_Sorted_Filtered.bed -b "${id}" > "coverage_$(basename "$id")"
done

## move columns to peaks_sorted.bed and make new file.
echo "\n# Make coverage table for the peaks spanning all samples." | tee -a ${BASE_NAME}_logs.txt
cp peaks_Sorted_Filtered.bed ${BASE_NAME}_peakCoverage.txt
echo "chr\tstart\tend\tname\tscore\tstrand" > colnames.txt
cat ${BASE_NAME}_peakCoverage.txt >> colnames.txt
mv colnames.txt ${BASE_NAME}_peakCoverage.txt

for cov_file in coverage_*; do
    echo $(echo $(basename "${cov_file}" .sorted.collapsed.bed) | cut -d '_' -f 2-) > temp1.txt
    awk 'FNR>0 {print $7}' ${cov_file} > temp2.txt && cat temp1.txt temp2.txt > temp3.txt
    paste ${BASE_NAME}_peakCoverage.txt temp3.txt > temp4.txt && mv temp4.txt ${BASE_NAME}_peakCoverage.txt
    rm temp*.txt
done

mv peaks_Sorted_Filtered.bed ${BASE_NAME}_peaks
mv coverage_* peakCoverage
mv peakCoverage ${BASE_NAME}_peaks
mv ${BASE_NAME}_peakCoverage.txt ${BASE_NAME}_peaks
mv ${BASE_NAME}_peaks ${BASE_NAME}

echo "\nAnalysis Finished: $(date +'%Y:%m:%d:%H:%M')" | tee -a ${BASE_NAME}_logs.txt
mv ${BASE_NAME}_logs.txt ${BASE_NAME}
