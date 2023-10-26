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
  echo "-----------------------------------------------------------------------------------------------"
  echo "CLIPittyClip: Single-line CLIP data analysis pipeline"
  echo "-----------------------------------------------------------------------------------------------"
  echo "Version 1.0.0"
  echo "Author: Soon Yi"
  echo "Last update: 2023-10-25"
  echo "-----------------------------------------------------------------------------------------------"
  echo "This is a part of the CLIPittyCLIP analysis pipeline, from bed files to peaks."
  echo "This script will combine all bed files and call peaks using HOMER."
  echo "The pipeline utilizes following programs: "
  echo " - bedtools (coverage)"
  echo " - Homer (makeTagDirectory, findPeaks)"
  echo "Make sure conda environment has all the necessary programs installed."
  echo "-----------------------------------------------------------------------------------------------"
}

function show_usage {
  echo "-----------------------------------------------------------------------------------------------"
  echo "Usage: CombinedPeakCalling.sh -p ? -z ? -f ?"
  echo "- If no option is selected, the program will run with default values (see below)."
  echo "- Make a folder named "BED" that contains all bed files that you want to call peaks on."
  echo "- Run this program inside the directory that contains "BED" folder."
  echo "-----------------------------------------------------------------------------------------------"
  echo "Options:"
  echo "  -h: print usage information"
  echo "  -p: minimum distance between peaks for homer            (default: 50)"
  echo "  -z: size of peaks for homer                             (default: 20)"
  echo "  -f: fragment length for homer                           (default: 25)"
  echo "-----------------------------------------------------------------------------------------------"
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
  echo "No .bed files found in BED folder. Exiting."
  on_exit
  exit 1
fi

## Combine all bed files:
echo -e "\n# Combine BED files to make global BED file for the experiment." | tee -a ${BASE_NAME}_logs.txt
echo -e "  Combining $(ls BED |wc -l) files to a single bed file..." | tee -a ${BASE_NAME}_logs.txt
ls -l BED | tee -a ${BASE_NAME}_logs.txt > /dev/null
cat BED/*.bed > ${BASE_NAME}.bed 

## Make Tag Directory and call peaks:
echo -e "\n# Use HOMER to call peaks, with minimum distance of ${PeakMinDist}, peak size of ${PeakSize}, and fragment length of ${FragLength}" | tee -a ${BASE_NAME}_logs.txt
makeTagDirectory ${BASE_NAME}_TagDir/ ${BASE_NAME}.bed -single -format bed 2>&1 | tee -a ${BASE_NAME}_logs.txt 
findPeaks ${BASE_NAME}_TagDir/ -o auto -style factor -L 2 -localSize 10000 -strand separate -minDist ${PeakMinDist} -size ${PeakSize} -fragLength ${FragLength} 2>&1 | tee -a ${BASE_NAME}_logs.txt 

## Process output peaks file to create peaks bed file:
echo "# Process peaks.txt to generate peaks BED file." | tee -a ${BASE_NAME}_logs.txt 
sed '/^[[:blank:]]*#/d;s/#.*//' ${BASE_NAME}_TagDir/peaks.txt > peaksTemp.bed
awk 'OFS="\t" {print $2, $3, $4, $1, $6, $5}' peaksTemp.bed > peaks.bed
rm peaksTemp.bed 
sort -k 1,1 -k2,2n peaks.bed > peaks_Sorted.bed
grep -E '^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)\s' < peaks_Sorted.bed > peaks_Sorted_Filtered.bed

mkdir ${BASE_NAME}_peaks
mkdir peakCoverage
mv ${BASE_NAME}.bed ${BASE_NAME}_peaks
mv ${BASE_NAME}_TagDir ${BASE_NAME}_peaks
mv peaks.bed ${BASE_NAME}_peaks
mv peaks_Sorted.bed ${BASE_NAME}_peaks

## 
echo -e "\n# Calculate coverages for the peaks for each sample using bedtools coverage." | tee -a ${BASE_NAME}_logs.txt 
for id in collapsedBED/*.bed; do
    bedtools coverage -s -a peaks_Sorted_Filtered.bed -b "${id}" > "coverage_$(basename "$id")"
done

## move columns to peaks_sorted.bed and make new file.
echo -e "\n# Make coverage table for the peaks spanning all samples." | tee -a ${BASE_NAME}_logs.txt 
cp peaks_Sorted_Filtered.bed ${BASE_NAME}_peakCoverage.txt
echo -e "chr\tstart\tend\tname\tscore\tstrand" > colnames.txt
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