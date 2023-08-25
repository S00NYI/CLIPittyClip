PeakMinDist=50
PeakSize=20
FragLength=25
BASE_NAME=Combined

## Combine all bed files:
echo -e "\n# Combine collapsed BED files to make global BED file for the experiment." | tee -a ${BASE_NAME}_logs.txt
echo -e "  Combining $(ls collapsedBED |wc -l) files to a single bed file..." | tee -a ${BASE_NAME}_logs.txt
ls -l collapsedBED | tee -a ${BASE_NAME}_logs.txt > /dev/null
cat collapsedBED/*.bed > ${BASE_NAME}.bed 

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