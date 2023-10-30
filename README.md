# CLIPittyClip: Single-line CLIP data analysis pipeline
- Version 1.0
- Author: Soon Yi
- Last update: 2023-10-25
-------------------------------------------------------------------------------------------------------------------
CLIPittyClip is a single command line CLIP analysis pipeline, from fastq.gz to peaks.

The pipeline utilizes following programs: 
 - fastx_toolkit (fastq_quality_filter, collapser, barcode_splitter, trimmer, clipper)
 - ctk (stripBarcode, tag2collapse)
 - bowtie2
 - samtools (view, sort, index)
 - bedtools (bamtobed, genomecov, coverage)
 - Homer (makeTagDirectory, findPeaks)

Make sure conda environment has all the necessary programs installed.  
For more information, see CLIPittyClip_environment.yml and CLIPittyclip-CondaInfo.txt.

-------------------------------------------------------------------------------------------------------------------
CLIPittyClip has been tested in:
- Apple MacbOok Pro (2016) Intel Core i7 Montery 12.6
- Apple MacBook Pro M2 Montery 12.6
- Apple MacBook Air M2 Ventura 13.5
- Ubuntu 20.04.6 LTS 64bit
-------------------------------------------------------------------------------------------------------------------
## Installation:

Use the CLIPittyClip_environment.yml to install necessary packages:

    conda env create -n CLIPittyClip_env -f CLIpittyClip_environment.yml

For Apple silicon, you may need to add CONDA_SUBDIR=osx-64 in front of the conda commands:

    CONDA_SUBDIR=osx-64 conda env create -n CLIPittyClip_env -f CLIpittyClip_environment.yml

***Note***: Using mamba instead of conda will be much faster.

Activate the created environment and install Homer by following instructions on Homer website.  
CLIPittyClip_env environment created will have all the necessary prerequisites to install Homer.
> http://homer.ucsd.edu/homer/introduction/install.html

Make sure to add Homer to PATH variable.  

Finally, add directory in which CLIPittyClip.sh is located in to PATH variable by running install_zshrc.sh (or install_bashrc.sh) file.  
Use one that matches your Shell. Modern Apple products use zsh as default. install.sh file can be deleted afterwards.  

Restart terminal and activate the CLIPittyClip_env before using the CLIPittyClip pipeline.  
```
conda activate CLIPittyClip_env
```
 
-------------------------------------------------------------------------------------------------------------------
## CLIPittyClip.sh Options:
 - -h: print usage information
 - **-i: experiment ID**
 - **-y: experiment type (e.g., Input, Enrich, or Fraction)**
 - -q: quality score threshold                             (default: 30)
 - -d: barcode file for demultiplexing.                    (default: no)
 - -b: allowed barcode mismatch for demultiplexing         (default:  0)
 - -l: minimum read length after 5/3 end trimming/clipping (default: 16)
 - **-x: path to genome index for bowtie2 mapping**
 - -m: allowed mismatch for mapping                        (default:  0)
 - -s: seed length for mapping                             (default: 15)
 - -t: number of threads to utilize for mapping            (default:  1)
 - -p: minimum distance between peaks for homer            (default: 50)
 - -z: size of peaks for homer                             (default: 20)
 - -f: fragment length for homer                           (default: 25)

Options in **bold** are required.

-------------------------------------------------------------------------------------------------------------------
## Usage:
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
## Note on read architecture:
The read architecture is as following:

5' - **UMI** | **Barcode** | **CCC** | **Reads** | **L32** - 3'

- 7nt long unique molecular identifier (UMI) used to collapse duplicate reads.
- 6nt long barcode used for library multiplexing.
- CCC Spacer
- Reads
- L32 RNA Linker

