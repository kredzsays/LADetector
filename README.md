­­­­

Vignette for DamID Sequencing Analysis

And LADetector algorithm for LAD Segmentation

For a detailed description of the DamID-sequencing protocol and the benchmarking of our program, please refer to:

“*\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_* “

If you use our program to analyze your data in a published work, please cite the above paper in the main text of your publication.

For questions about usage, please email Xianrong Wong: xwong2@jhu.edu, Teresa Luperchio: trl@jhmi.edu, or Karen Reddy: kreddy4@jhmi.edu.

**Introduction**

DNA Adenine Methyltransferase Identification, DamID, has been commonly used to probe lamina chromatin interactions in cells. The detailed description of the DamID protocol has been described in \[1\]. For a detailed description of the adaptation for sequencing, please refer to the paper on the cover page of this vignette. In brief, the end point of DamID amplification is a pool of DNA amplicons flanked by DamID adapters/primers which range typically between 200bp – 3kb. To obtain a library size of between 100-300bp for efficient flow cell clustering, the DamID amplicons are randomized by inter-fragment ligation followed by sonication, to prevent possible preferential loss of smaller DNA fragments by straight up sonication of the DamID amplicons. This, however, will generate DNA fragments with some probability of harboring DamID primer sequences which have to be bioinformatically removed prior to mapping to a reference index.

This program can be divided into three parts. The first part, pre-normalization mapping, deals with the removal of primer sequences and mapping of the reads to a reference genome. This is followed by counting the number of reads that have been mapped to each genomic bin, genome wide. The second part of the program normalizes the counts in each bin to the sequencing depth of both the experimental sample and the control sample. The normalized experimental scores are then further normalized to the control sample scores to yield genome wide log2-ratios. Part three of the program uses the LADetector version 2 algorithm to extract genomic intervals for LADs and DIPs from the log2-ratios (version 1 was published and used and provided in Harr JC, Luperchio TR, Wong X, Cohen E, Wheelan SJ, Reddy KL “Directed targeting of chromatin to the nuclear lamina is mediated by chromatin state and A-type lamins”. JCB 2015 January 5).

**2. Installation**

Usage of the LADetector requires downloading all files and scripts. The full set of files can be found in DamID-LADetector.tgz. Examples are provided.

Prerequisites include:

2.1 Installation of bowtie, bedtools, samtools and R (for LADetector). More information on this can be found at:

-   http://bowtie-bio.sourceforge.net/tutorial.shtml

-   <http://bedtools.readthedocs.org/en/latest/content/installation.html>

-   <http://samtools.sourceforge.net/>

-   https://www.r-project.org/

2.2 A folder with pre-built indexes for mapping (this workflow uses bowtie). Some pre-built indexes (for bowtie) can be found in <http://bowtie-bio.sourceforge.net/tutorial.shtml>

2.3 Have a copy of the genomic bins in the unpacked DamID-LADetector folder. In our analysis, we used genomic bins flanked by GATC sites (DPNI and DPNII digestion sites) from the mm9 build. We have provided the bins for mm9 and hg38 in the DamID-LADetector folder.

2.4 Define environment variables. Each script is defined in Usage section, and directions for setting environment variables are provided in header of each script.

**3. Usage**

**3.1 Pre-normalization mapping**

This part of the program is illustrated and described in figure 1.

To set environment variables for pre-normalization-scoring.sh

1.  Replace {PATH\_TO\_DamID-LADetector\_folder} to the pathlength of the DamID-LADetector folder.

2.  Replace {PATH\_TO\_reference\_index\_folder} to the pathlength of the folder where your reference index, eg mm9, is. If using another reference index other than mm9, replace \[mm9\] to the index you are using eg hg19.

3.  Replace {PATH\_TO\_genomic\_bins} with the pathlength of your genomic\_bin.bed file. For example, we used the Dpn.bed file and we would change the pathlength to where the file actually is, ~/scratch/DamID-LADetector/Dpn.bed

-   To run the pre-normalization-scoring.sh script on a node in slurm, type:

bash path\_to\_pre-normalization-scoring.sh path\_to\_input\_fastq 1

-   To run on a regular shell, type:

bash path-to-pre-normalization-scoring.sh path\_to\_input\_fastq

The number 1 at the end of the first command simply tells the program that you are running on slurm and certain modules need to be loaded.

-   To submit a job to slurm, type:

sbatch –t hh:mm:ss --mem=10G path\_to\_pre-normalization-scoring.sh path\_to\_input\_fastq 1

Specifying a memory 10gb is sufficient to run the script on a 40-50 million read fastq input.

<img src="./media/image1.jpg" width="576" height="842" />

**3.2 Normalization and obtaining genome wide log2-ratios**

After 3.1 has been run, you will find several files in the folder with your input fastq file, 2 of which are needed for the Normalization.pl script.

To run this use:

perl path-to-Normalization.pl path-to-Dam.preNormalization.score path-to-LmnB1. .preNormalization.score path-to-Dam.mappedReadCounts path-to-LmnB1.mappedReadCounts path-where-you-want-output.bed

The output of this program yields the normalized log2ratios bed file.

**3.3 LADetector**

The central part of the algorithm is illustrated and described in figure 2.

To set environment variable for LADs\_and\_DIPs.sh,

1.  Replace {PATH\_TO\_DamID-LADetector\_folder} to the pathlength of the DamID-LADetector folder.

2.  Replace {PATH\_TO\_mm9.genome} to the pathlength of the mm9.genome file.

> This file stores the information of the sizes of all the chromosomes in the mm9 build. The format of the file is shown below:
>
> chr1 197195432
>
> chr2 181748087
>
> chr3 159599783
>
> chr4 155630120

.

.

.

To run, type:

bash path\_to\_LADs\_and\_DIPs.sh path\_to\_normalized.bed min\_DIP\_size max\_DIP\_size

For default LADetector\_III settings, enter 0 and 0 for min\_DIP\_size and

<img src="./media/image2.jpg" width="564" height="783" />**max\_DIP\_size respectively. This would set min\_DIP\_size to 2000 and max\_DIP\_size** to 7000.

The min\_DIP\_size and max\_DIP\_size defines the minimum and maximum sizes for a negative region within a LAD to be recognized as a dip. As shown in figure 2b, increasing the max\_DIP\_size increases the number of dips reported. When a dip is reported, the flanking positive LADetectorII outputs are stitched together indicating that THAT dip is a feature of a LAD. If the negatively associated region exceeds the max\_DIP\_size, the flanking LADs do not get stitched together indicating that the negative region is defined as an interLAD.

The important output files for LADetector will be appended with .LADs, .DIPs and .bedgraph where .LADs stores intervals for LADs, .DIPs stores intervals for DIPs and .bedgraph stores the log2-ratios with unalignable regions removed.

**4. Galaxy Implementation**

The aforementioned analysis can be performed using pre-packaged and custom packaged tools in galaxy. We provide several workflows that will suit the needs of a range of end-users.

Where to find “working” image of the analysis?

**4.1 Comprehensive workflow**

This workflow strings together all 3 parts of the analysis (pre-normalization scoring, normalization and LAD detection) as shown in figure 3. Users need only upload the 2 input fastq files (experimental and normalizer) and the genomic bins of choice.

<img src="./media/image3.png" width="564" height="612" />

**4.2 Parsed Workflow**

**4.2.1 Scoring and normalization**

This workflow requires the same inputs as the comprehensive workflow (Figure 3). The endpoint of the workflow is the generation of a normalized log2ratio bedgraph file and an unalignable regions file.

**4.2.2 LADetection by LADetector**

<img src="./media/image5.png" width="576" height="420" />This workflow uses the custom built LADetector packages that binarizes log2ratio inputs into LAD intervals and DIP intervals. As shown in figure 4, the inputs to this workflow are the log2ratios generated from 4.2.1 and the unalignable regions bed file also generated in 4.2.1. This workflow exists as a separate entity to enhance its versatility so as to facilitate its use as a downstream tool that caters to a wide range sequencing analysis and normalization techniques.
