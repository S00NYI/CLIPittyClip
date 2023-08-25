-------------------------------------------------------------------------------------------------------------------
CLIPittyClip: Single-line CLIP data analysis pipeline
Version 0.1
Author: Soon Yi
Last update: 2023-08-06
-------------------------------------------------------------------------------------------------------------------
This is a standard CLIP analysis pipeline, from fastq.gz to peaks.

The pipeline utilizes following programs: 
 - fastx_toolkit (collapser, barcode_splitter, trimmer, clipper)
 - ctk (stripBarcode, tag2collapse)
 - bowtie2
 - samtools (view, sort, index)
 - bedtools (bamtobed, genomecov)
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
Installation:

Use the CLIpittyClip_environment.yml to install necessary packages:
    conda env create -n CLIPittyClip_env -f CLIpittyClip_environment.yml

For Apple silicon, you may need to add CONDA_SUBDIR=osx-64 in front of the conda commands:
    CONDA_SUBDIR=osx-64 conda env create -n CLIPittyClip_env -f CLIpittyClip_environment.yml

** Using mamba instead of conda will be much faster.

Activate the created environment and install Homer by following instructions on Homer website:
    http://homer.ucsd.edu/homer/introduction/install.html
    CLIPittyClip_env environment created will have all the necessary prerequisites to install Homer.

Finally, add directory in which CLIPittyClip.sh is located in to PATH variable by running install_zshrc.sh (or install_bashrc.sh) file.
    Use one that matches your Shell. Modern Apple products use zsh as default.
    install.sh file can be deleted afterwards.

-------------------------------------------------------------------------------------------------------------------
Restart terminal and CLIPittyClip.sh can be ran as following:

CLIPittyClip.sh -i PROJECT_ID -y PROJECT_TYPE -x /PATH/TO/GENOME/ANNOTATION/ANNOTATION_FILE_NAME

Options -i, -y, and -x are required.
    For more detailed usage information, type CLIPittyClip.sh in terminal after installation and press enter.
    Make sure to activate the conda environment.

-x specifies path to Bowtie2 index
    Bowtie2 index files can be made by using bowtie2-build and supplying genome fasta files.
        Here are example Bowtie2 index files made using GRCh38.primary_assembly.genome.fa downloaded from GENCODE.
    Bowtie2 index files consist of 6 files: 
        GRCh38.primary_assembly.genome.1.bt2
        GRCh38.primary_assembly.genome.2.bt2
        GRCh38.primary_assembly.genome.3.bt2
        GRCh38.primary_assembly.genome.4.bt2
        GRCh38.primary_assembly.genome.rev.1.bt2
        GRCh38.primary_assembly.genome.rev.2.bt2

    For this example, -x option should be specified as:
        /PATH/TO/GENOME/ANNOTATION/GRCh38.primary_assembly.genome

Genome file for bedtools genomecov should have same path and name scheme as bowtie2 index: 
    /PATH/TO/GENOME/ANNOTATION/GRCh38.primary_assembly.genome.fa.fai
    This .fa.fai file can be easily made by: 
        samtools faidx /PATH/TO/GENOME/ANNOTATION/GRCh38.primary_assembly.genome.fa

-------------------------------------------------------------------------------------------------------------------



