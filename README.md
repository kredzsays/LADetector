#DamID Sequencing Analysis and LADetector Algorithm for LAD Segmentation
For a detailed description of the DamID-sequencing protocol and the benchmarking of our program, please refer to:
*“LADetector: a rapid, large-domain calling algorithm that allows for accurate identification of edges and dips with both ChIP and DamID technologies“* Wong and Luperchio et al.

*If you use our program to analyze your data in a published work, please cite the above paper in your publication.*

For questions about usage, please email Xianrong Wong: xwong2@jhu.edu, Teresa Luperchio: trl@jhmi.edu, or Karen Reddy: kreddy4@jhmi.edu.

## Introduction

DNA Adenine Methyltransferase Identification, DamID, has been commonly used to probe lamina chromatin interactions in cells. The detailed description of the DamID protocol is found in *Vogel et al, Detection of in vivo protein-DNA interactions using DamID in mammalian cells. Nature Protocols 2007*. For a detailed description of the adaptation for sequencing, please refer to our paper. In brief, the end point of DamID amplification is a pool of DNA amplicons flanked by DamID adapters/primers which range typically between 200bp – 3kb. To obtain a library size of between 100-300bp for efficient flow cell clustering, the DamID amplicons are randomized by inter-fragment ligation followed by sonication, to prevent possible preferential loss of smaller DNA fragments by straight up sonication of the DamID amplicons. This, however, will generate DNA fragments with some probability of harboring DamID primer sequences which have to be bioinformatically removed prior to mapping to a reference index.

This program can be divided into three parts; **mapping**, **normalization** and **LAD detection**. The first part, pre-normalization mapping, deals with the removal of primer sequences and mapping of the reads to a reference genome. This is followed by counting the number of reads that have been mapped to each genomic bin, genome wide. The second part of the program normalizes the counts in each bin to the sequencing depth of both the experimental sample and the control sample. The normalized experimental scores are then further normalized to the control sample scores to yield genome wide log2-ratios. Part three of the program uses the LADetector version 2 algorithm to extract genomic intervals for LADs and DIPs from the log2-ratios (version 1 was published and provided in *Harr JC, Luperchio TR, Wong X, Cohen E, Wheelan SJ, Reddy KL “Directed targeting of chromatin to the nuclear lamina is mediated by chromatin state and A-type lamins”. JCB 2015 January 5*).

## Command line tool
### 1. Installation

  1.1 Usage of the LADetector requires downloading all files and scripts. The full set of scripts necessary can be found in [LADetector_8122016.tar.gz] (https://github.com/thereddylab/LADetector/raw/master/LADetector_8122016.tar.gz). Some necessary files required for including alignable and unalignable regions are included in the repo. 

  1.2 Installation of [bowtie] (http://bowtie-bio.sourceforge.net/tutorial.shtml), [bedtools] (http://bedtools.readthedocs.org/en/latest/content/installation.html), [samtools] (http://samtools.sourceforge.net/) and [R] (https://www.r-project.org/) (for LADetector) is required. These are executables, please be sure to include their locations in your PATH. For further help with your PATH, please see [this page as an example] (http://www.computerhope.com/issues/ch001647.htm).

 1.3 Create a folder to contain pre-built indexes for mapping (this workflow uses bowtie). Pre-built indexes for many genome builds that are specific for bowtie can be downloaded from the [bowtie website] (http://bowtie-bio.sourceforge.net/tutorial.shtml). You must have these for correct mapping. 

 1.4 Have a copy of the genomic bins in the unpacked LADetector folder (available for [download] (https://github.com/thereddylab/LADetector/raw/master/LADetector_8122016.tar.gz) as mentioned above). In our analysis, we used genomic bins flanked by GATC sites (DPNI and DPNII digestion sites) from the mm9 build. We have provided the bins for mm9 and hg38 downloadble from the repository. If you are using other genomes or builds you will need to provide them yourself.

 1.5 Define environment variables. Each script is defined in Usage section, and directions for setting environment variables are provided in the header of each script.

1.6 Samples to be analyzed should be stored on your system in uncompressed fastq format, one file for control (Dam) and one for experiment (LaminB1).

### 2. Usage

 2.1 Pre-normalization mapping

This part of the program is illustrated and described in figure 1.

![ladetector_figure1_1](https://cloud.githubusercontent.com/assets/17512466/16018097/e5e44ff2-3170-11e6-8e4f-d7d8aaf7300a.png)

You must set environment variables within *pre-normalization-scoring.sh* before running, directions are in the header:

    ######################################
    TO DO before you run this script
    SCRIPTS=PATH_TO_DamID-LADetector_folder
    #Enter the path of your scripts folder/1_fastq_quality_trimmer.pl.
    BUILD=PATH_TO_reference_index_folder
    #Enter the path of your reference index folder.
    #mm9 and hg38 are provided
    BINS=PATH_TO_genomic_bins
    #Enter the path to your genomic_bins.bed file.  

For example:

  
    SCRIPTS=/home-4/trl@jhmi.edu/scratch/ladetector/scripts/1_fastq_quality_trimmer.pl
    #Enter the path of your scripts folder/1_fastq_quality_trimmer.pl.
    BUILD=/home-4/trl@jhmi.edu/scratch/genomes/mm9/mm9
    #Enter the path of your reference index folder.
    #mm9 and hg38 are provided
    BINS=/home-4/trl@jhmi.edu/scratch/ladetector/scripts/Dpn.bed
    #Enter the path to your genomic_bins.bed file.  


To run the *pre-normalization-scoring.sh* script on a node in **slurm**, use:

    bash path_to_pre-normalization-scoring.sh path_to_input_fastq 1

The number '1' is necessary for running on slurm

To run on a **regular shell**, use:

    bash path-to-pre-normalization-scoring.sh path_to_input_fastq


Specifying a memory 10gb is sufficient to run the script on a 40-50 million read fastq input.

2.2 Normalization and obtaining genome wide log2-ratios

After *path-to-pre-normalization-scoring.sh* has been run, you will find several files in the folder with your input fastq file, 2 of which are needed for the *Normalization.pl* script.

To run this use:

    perl path-to-Normalization.pl path-to-Dam.preNormalization.score path-to-LmnB1. .preNormalization.score path-to-Dam.mappedReadCounts path-to-LmnB1.mappedReadCounts path-where-you-want-output.bed

The output of this program yields the normalized log2ratios bed file.

2.3 LADetector

The central part of the algorithm is illustrated and described in figure 2.

![ladetector_figure2_1](https://cloud.githubusercontent.com/assets/17512466/16018096/e5dd9662-3170-11e6-9d3a-794dc675b05c.png)

You must set environment variables within *LADs\_and\_DIPs.sh* 

    ######################################
    TO DO before you run this script
    SCRIPTS=PATH_TO_DamID-LADetector_folder
    #Enter the path of LADs_andDIPs.sh
    GENOME=PATH_TO_mm9.genome
    #Enter the path to mm9.genome file - format:   chr    chrSize   
    UNALIGNABLE=PATH_TO_mm9.unalignable.genome
    #Enter the path to mm9.unalignable.genome file

For example:

    SCRIPTS=/home-4/trl@jhmi.edu/scratch/ladetector/scripts/LADetector_scripts/LADs_and_DIPs.sh
    #Enter the path of LADs_andDIPs.sh
    GENOME=/home-4/trl@jhmi.edu/scratch/ladetector/scripts/LADetector_scripts/mouse.mm9.genome
    #Enter the path to mm9.genome file - format:   chr    chrSize   
    UNALIGNABLE=/home-4/trl@jhmi.edu/scratch/ladetector/scripts/LADetector_scripts/mm9.unalignable.txt

The GENOME file stores the information of the sizes of all the chromosomes in the mm9 build. The format of the file is:

    chr1 197195432
    chr2 181748087
    chr3 159599783
    chr4 155630120
    .
    .
    .

Additionally, you must have R running on your system

To run, use:

    bash path_to_LADs_and_DIPs.sh path_to_normalized.bed min_DIP_size max_DIP_size

For default LADetector settings, enter 0 and 0 for min\_DIP\_size and max\_DIP\_size.
This would set min\_DIP\_size to 2000 and max\_DIP\_size to 7000.

For example: 

    bash path_to_LADs_and_DIPs.sh path_to_normalized.bed 0 0

The min\_DIP\_size and max\_DIP\_size defines the minimum and maximum sizes for a negative region within a LAD to be recognized as a Dip. As shown in figure 2b, increasing the max\_DIP\_size increases the number of Dips reported. When a Dip is reported, the flanking positive LADetectorII outputs are stitched together indicating that THAT Dip is a feature of a LAD. If the negatively associated region exceeds the max\_DIP\_size, the flanking LADs do not get stitched together indicating that the negative region is defined as an interLAD.

The important output files for LADetector will be appended with **.LADs**, **.DIPs** and **.bedgraph** files. 

* **.LADs** stores intervals for **LADs** in bed format
* **.DIPs** stores intervals for **DIPs** in bed format
* **.bedgraph** stores the **log2-ratios** with unalignable regions removed.

Updated 12/1/2016 TRL
