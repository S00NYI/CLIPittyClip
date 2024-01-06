terminal_width=$(tput cols)
separator_line=$(printf "%${terminal_width}s" | tr ' ' '*')

# Peak calling minimum distance for HOMER 
PeakMinDist=50
# Peak size for HOMER
PeakSize=20
# Fragment size for HOMER
FragLength=25
# Base name for file
BASE_NAME=Combined

# Function to display information
function show_info {
  echo "$separator_line"
  echo "CLIPittyClip: Single-line CLIP data analysis pipeline"
  echo "$separator_line"
  echo "Version 2.0.0"
  echo "Author: Soon Yi"
  echo "Last updated: 2024-01-05"
  echo "$separator_line"
  echo "PEAKittyPeak.sh"
  echo "$separator_line"
  echo "This is the combined peak calling part of CLIPittyClip."
  echo "Use this sub-pipeline to call peaks on a combined bed file (e.g., from demultiplexed samples)."
  echo "The pipeline utilizes the following programs: "
  echo " - bedtools (coverage)"
  echo " - Homer (makeTagDirectory, findPeaks)"
  echo "$separator_line"
}

function show_usage {
  echo "Usage: PEAKittyPeak.sh -p ? -z ? -f ?"
  echo "$separator_line"
  echo "- If no option is selected, the program will run with default values (see below)."
  echo "- Make a folder named \"BED\" that contains all bed files that you want to call peaks on."
  echo "  Run this program inside the directory that contains \"BED\" folder."
  echo "  The program will then make a bed file that combines all the provided bed files."
  echo "  Peak calling will be performed using the combined bed file."
  echo "$separator_line"
  echo "Options:"
  echo "  -h: print usage information"
  echo "  -p: minimum distance between peaks for homer            (default: 50)"
  echo "  -z: size of peaks for homer                             (default: 20)"
  echo "  -f: fragment length for homer                           (default: 25)"
  echo "$separator_line"
}

function on_exit {
  read -p "Press enter to exit."
}

while getopts "p:z:f:h" opt; do
  case $opt in
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

# Count the number of .bed files in the current directory
bedFileCount=$(ls BED |wc -l)

# If no .bed files are found, print a message and exit
if [[ $bedFileCount -eq 0 ]]; then
  echo "No .bed files found in the BED folder. Exiting."
  on_exit
  exit 1
fi

# Start Pipeline
show_info | tee -a ${BASE_NAME}_peaks_logs.txt
echo "# Analysis started: $(date +'%Y/%m/%d %H:%M')" | tee -a ${BASE_NAME}_peaks_logs.txt
echo "# User Input: '$0 $@'\n" | tee -a ${BASE_NAME}_peaks_logs.txt

# Check if arguments are provided:
if [[ $# -eq 0 ]]; then
  echo "# No option arguments have been passed. Default options will be used."
fi

## Combine all bed files:
echo "\n# Combine BED files to make global BED file for the experiment." | tee -a ${BASE_NAME}_peaks_logs.txt
echo "  Combining $(ls BED |wc -l) files to a single bed file..." | tee -a ${BASE_NAME}_peaks_logs.txt
ls BED | tee -a ${BASE_NAME}_peaks_logs.txt > /dev/null
cat BED/*.bed > ${BASE_NAME}.bed 

## Make Tag Directory and call peaks:
echo "\n# Use HOMER to call peaks, with minimum distance of ${PeakMinDist}, peak size of ${PeakSize}, and fragment length of ${FragLength}" | tee -a ${BASE_NAME}_peaks_logs.txt
makeTagDirectory ${BASE_NAME}_TagDir/ ${BASE_NAME}.bed -single -format bed 2>&1 | tee -a ${BASE_NAME}_peaks_logs.txt 
findPeaks ${BASE_NAME}_TagDir/ -o auto -style factor -L 2 -localSize 10000 -strand separate -minDist ${PeakMinDist} -size ${PeakSize} -fragLength ${FragLength} 2>&1 | tee -a ${BASE_NAME}_peaks_logs.txt 

## Process output peaks file to create peaks bed file:
echo "# Process peaks.txt to generate peaks BED file." | tee -a ${BASE_NAME}_peaks_logs.txt 
sed '/^[[:blank:]]*#/d;s/#.*//' ${BASE_NAME}_TagDir/peaks.txt > peaksTemp.bed
awk 'OFS="\t" {print $2, $3, $4, $1, $6, $5}' peaksTemp.bed > peaks.bed
rm peaksTemp.bed 
sort -k 1,1 -k2,2n peaks.bed > peaks_Sorted.bed
# grep -E '^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)\s' < peaks_Sorted.bed > peaks_Sorted_Filtered.bed

mkdir ${BASE_NAME}_peaks
mkdir peakCoverage
mv ${BASE_NAME}.bed ${BASE_NAME}_peaks
mv ${BASE_NAME}_TagDir ${BASE_NAME}_peaks
mv peaks.bed ${BASE_NAME}_peaks
# mv peaks_Sorted.bed ${BASE_NAME}_peaks

## 
echo "\n# Calculate coverages for the peaks for each sample using bedtools coverage." | tee -a ${BASE_NAME}_peaks_logs.txt 
for id in BED/*.bed; do
    bedtools coverage -s -a peaks_Sorted.bed -b "${id}" > "coverage_$(basename "$id")"
done

## move columns to peaks_sorted.bed and make new file.
echo "\n# Make coverage table for the peaks spanning all samples." | tee -a ${BASE_NAME}_peaks_logs.txt 
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

echo "\n# Final output: ${BASE_NAME}_peakCoverage.txt" | tee -a ${BASE_NAME}_peaks_logs.txt 

mv peaks_Sorted.bed ${BASE_NAME}_peaks
mv coverage_* peakCoverage
mv peakCoverage ${BASE_NAME}_peaks
mv ${BASE_NAME}_peakCoverage.txt ${BASE_NAME}_peaks
mv BED ${BASE_NAME}_peaks

echo "\n# Analysis Finished: $(date +'%Y/%m/%d %H:%M')" | tee -a ${BASE_NAME}_peaks_logs.txt
echo "\n# All outputs are in ${BASE_NAME}_peaks folder." | tee -a ${BASE_NAME}_peaks_logs.txt
mv ${BASE_NAME}_peaks_logs.txt ${BASE_NAME}_peaks