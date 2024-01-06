# CLIPittyClip: Single-line CLIP data analysis pipeline
- Version 2.0
- Author: Soon Yi
- Last updated: 2024-01-05
-------------------------------------------------------------------------------------------------------------------
CLIPittyClip is a single command line CLIP analysis pipeline, from fastq.gz to peaks.

The pipeline utilizes following programs: 
 - fastx_toolkit (fastq_quality_filter, collapser, barcode_splitter, trimmer, clipper)
 - ctk (stripBarcode, tag2collapse)
 - bowtie2
 - samtools (view, sort, index)
 - bedtools (bamtobed, genomecov, coverage)
 - Homer (makeTagDirectory, findPeaks)
  
-------------------------------------------------------------------------------------------------------------------
CLIPittyClip has been tested in:
- Apple MacBook Pro (2016) Intel Core i7 Montery 12.6
- Apple MacBook Pro M2 Montery 12.6
- Apple MacBook Air M2 Ventura 13.5
- Ubuntu 20.04.6 LTS 64bit
-------------------------------------------------------------------------------------------------------------------
## Installation:

Use the CLIPittyClip_environment.yml to install necessary packages:
```
conda env create -n CLIPittyClip_env -f install_environment.yml
```
For MacBook with Apple silicon, you may want to try the following if above conda call does not work:
```
CONDA_SUBDIR=osx-64 conda env create -n CLIPittyClip_env -f install_environment.yml
```
***Note***: Using mamba instead of conda will be much faster.

Restart terminal and activate the CLIPittyClip_env: 
```
conda activate CLIPittyClip_env
```

Install Homer by following instructions on Homer website: http://homer.ucsd.edu/homer/introduction/install.html 

CLIPittyClip_env environment have all the necessary prerequisites to install Homer.

Also, make sure to add Homer to PATH variable.  

Lastly, add directory in which CLIPittyClip.sh is located in to PATH variable by running install_zshrc.sh (or install_bashrc.sh) file.  
```
bash /install_bashrc.sh
```
OR
```
bash /install_zshrc.sh
```
Restart terminal and activate the CLIPittyClip_env. You are now ready to use CLIPittyClip programs.  
```
conda activate CLIPittyClip_env
```
 
-------------------------------------------------------------------------------------------------------------------
## CLIPittyClip.sh 
### Options:

Options in **bold** are required.

| Options  |Descriptions                                                |Default|
| :-------:|------------------------------------------------------------|:-----:|
| h        | print usage information                                    |       |
| **i**    | **experiment ID**                                          |       |
| **y**    | **experiment type (e.g., Input, Enrich, or Fraction)**     |       |
| k        | keep intermediate files (yes/no)                           |no     |
| q        | quality score threshold (for fastx_quality_filter)         |30     |
| d        | perform demultiplexing (yes/no)                            |no     |
| b        | allowed number of mismatches in barcode for demultiplexing |0      |
| 5        | number of nucleotide to clip from the 5'-end               |10     |
| 3        | sequence of the 3'-end adapter                             |L32: GTGTCAGTCACTTCCAGCGG|
| l        | minimum read length after 5/3 end trimming/clipping        |16     |
| **x**    | **path to genome index**                                   |       |
| m        | allowed mismatch for mapping                               |0      |
| s        | seed length for mapping                                    |15     |
| t        | number of threads to be utilized by bowtie/samtools        |1      |
| p        | minimum distance between peaks for homer                   |50     |
| z        | size of peaks for homer                                    |20     |
| f        | fragment length for homer                                  |25     |

-------------------------------------------------------------------------------------------------------------------
### Usage:
On terminal, traverse to directory containing your fastq.gz file.

CLIPittyClip.sh can be ran as following:
```
CLIPittyClip.sh -i ID -y TYPE -x /PATH/TO/GENOME/ANNOTATION/ANNOTATION_FILE_NAME
```
Options **-i**, **-y**, and **-x** are required.  

-i and -y is used to define the input fasta.gz name.
- ID (*i*) is your experimental ID (e.g. SY1004).
- Type (*y*) is your short handle describing the experiemnt (e.g. HuR_CLIP).
- Then, input fastq.gz should have filename as: id_type.fastq.gz
 - For the examples above, input fastq.gz file should have file name as <ins>*SY1004_HuR_CLIP.fastq.gz*</ins>.

-x specifies path to index files for bowtie2 mapping and genome file for bedGraph generation:  
- Bowtie2 index files can be made by using bowtie2-build and supplying genome fasta files.
- Here are example Bowtie2 index files made using GRCh38.primary_assembly.genome.fa downloaded from GENCODE.
  - GRCh38.primary_assembly.genome.1.bt2
  - GRCh38.primary_assembly.genome.2.bt2
  - GRCh38.primary_assembly.genome.3.bt2
  - GRCh38.primary_assembly.genome.4.bt2
  - GRCh38.primary_assembly.genome.rev.1.bt2
  - GRCh38.primary_assembly.genome.rev.2.bt2
- For this example, -x option should be specified as:
  - /PATH/TO/GENOME/ANNOTATION/GRCh38.primary_assembly.genome

- Genome file for bedtools genomecov should also have same path and name scheme as bowtie2 index:
  - /PATH/TO/GENOME/ANNOTATION/GRCh38.primary_assembly.genome.fa.fai
  - This .fa.fai file can be easily made by:
    - samtools faidx /PATH/TO/GENOME/ANNOTATION/GRCh38.primary_assembly.genome.fa

-d option is used to turn on/off sample demultiplexing:
- If demultiplexing is required, make sure to set -d option to yes (default is no).
  - If -d option is set to 'yes', barcode file should be provided.
  - Barcode file format information: http://hannonlab.cshl.edu/fastx_toolkit/commandline.html#fastx_barcode_splitter_usage
  - Barcode file should have the same name convention as fastq.gz file, followed by '_BC.txt'.
   - For the example fastq.gz above, barcode file should have filename as <ins>*SY1004_HuR_CLIP_BC.txt*</ins>

-------------------------------------------------------------------------------------------------------------------
## MAPittyMap.sh 
### Options:

Options in **bold** are required.

| Options  |Descriptions                                                |Default|
| :-------:|------------------------------------------------------------|:-----:|
| h        | print usage information                                    |       |
| **i**    | **experiment ID**                                          |       |
| **y**    | **experiment type (e.g., Input, Enrich, or Fraction)**     |       |
| k        | keep intermediate files (yes/no)                           |no     |
| **x**    | **path to genome index**                                   |       |
| m        | allowed mismatch for mapping                               |0      |
| s        | seed length for mapping                                    |15     |
| t        | number of threads to be utilized by bowtie/samtools        |1      |
| p        | minimum distance between peaks for homer                   |50     |
| z        | size of peaks for homer                                    |20     |
| f        | fragment length for homer                                  |25     |

-------------------------------------------------------------------------------------------------------------------
### Usage:
Similar to CLIPittyClip.sh, but uses processed fasta files.

MAPittyMap.sh can be ran as following:
```
MAPittyMap.sh -i ID -y TYPE -x /PATH/TO/GENOME/ANNOTATION/ANNOTATION_FILE_NAME
```
Options **-i**, **-y**, and **-x** are required. 

See CLIPittyClip.sh usage section for more information.

- If using groomed.fasta output from CLIPittyClip.sh:
  - Copy *groomed.fasta to a separate location (not required, but good for organization.)
  - In the directory containing the *groomed.fasta, run MapittyMap.sh.
  - For option -i, use the same input as the original CLIPittyClip.sh run.
  - For option -y, add '_collapsed_rmBC_groomed' at the end of the original input.
  - For example: 
    - CLIPittyClip.sh -i JL100 -y Input -x path-to-genome-index
    - MapittyMap.sh -i JL100 -y Input_collapsed_rmBC_groomed -x path-to-diff-genome-index

-------------------------------------------------------------------------------------------------------------------
## PEAKittyPeak.sh 
### Options:

Options in **bold** are required.

| Options  |Descriptions                                                |Default|
| :-------:|------------------------------------------------------------|:-----:|
| h        | print usage information                                    |       |
| p        | minimum distance between peaks for homer                   |50     |
| z        | size of peaks for homer                                    |20     |
| f        | fragment length for homer                                  |25     |

-------------------------------------------------------------------------------------------------------------------
### Usage:
- Make a folder named \"BED\" that contains all bed files that you want to call peaks on.
  - Run this program inside the directory that contains \"BED\" folder (**not** inside the \"BED\" folder)
  - The program will then make a bed file that combines all the provided bed files.
  - Peak calling will be performed using the combined bed file.
-------------------------------------------------------------------------------------------------------------------
## Note on read architecture in CoCLIP Experiment:
The read architecture is as following:

5' - **UMI** | **Barcode** | **CCC** | **Reads** | **L32** - 3'

- 7nt long unique molecular identifier (UMI) used to collapse duplicate reads.
- 6nt long barcode used for library multiplexing.
- CCC Spacer
- Reads
- L32 RNA Linker

