----------------------------------------------------------------------
CLIPittyClip: Single-line CLIP data analysis pipeline
Version 0.1
Author: Soon Yi
Last update: 2023-08-06
----------------------------------------------------------------------
This is a standard CLIP analysis pipeline, from fastq.gz to peaks.

The pipeline utilizes following programs: 
 - fastx_toolkit (collapser, barcode_splitter, trimmer, clipper)
 - ctk (stripBarcode, tag2collapse)
 - bowtie2
 - samtools (view, sort, index)
 - bedtools (bamtobed, genomecov)
 - Homer (makeTagDirectory, findPeaks)

Make sure conda environment has all the necessary programs installed.
For more information, see CLIPittyClip-env.txt and CLIPittyclip-CondaInfo.txt.

----------------------------------------------------------------------
Installation:

Use the CLIPittyClip-env.txt to install necessary packages:
conda env create -n CLIPittyClip_env -f CLIpittyClip_environment.yml
CONDA_SUBDIR=osx-64 conda env create -n CLIPittyCLIP_ENV -f CLIpittyClip_environment.yml

Also, install Homer by following instructions on Homer website.

Finally, add directory in which CLIPittyClip.sh and install.sh is located in to PATH variable by running install.sh file.

install.sh file can be deleted afterwards.

----------------------------------------------------------------------
Restart terminal and CLIPittyClip.sh can be ran as following:

CLIPittyClip.sh -i JL0361 -y Input -x /PATH/TO/GENOME/ANNOTATION/GRCh38.primary_assembly.genome -b 2 

Bowtie2 index files can be made by using bowtie2-build and supplying genome fasta files.
Bowtie2 index files consist of 6 files: 
GRCh38.primary_assembly.genome.1.bt2
GRCh38.primary_assembly.genome.2.bt2
GRCh38.primary_assembly.genome.3.bt2
GRCh38.primary_assembly.genome.4.bt2
GRCh38.primary_assembly.genome.rev.1.bt2
GRCh38.primary_assembly.genome.rev.2.bt2

Genome file for bedtools genomecov will have same path and name scheme as bowtie2 index: /PATH/TO/GENOME/ANNOTATION/GRCh38.primary_assembly.genome.fa.fai
This file can be easily made by: samtools faidx /PATH/TO/GENOME/ANNOTATION/GRCh38.primary_assembly.genome.fa

Indices are supplied in indices.tar.gz file, which can be unziped using gunzip indices.tar.gz