# Seismic Amplification
Programs and scripts written for the manuscript **"Chromothripsis followed by circular recombination drives oncogene amplification in human cancer"**. For detailed information on the methods, please refer to the supplementary material.

## (1) Seismic amplification detection:

The R-script **seismic_amplification_detection.R** detects regions with seismic amplification from given copy number and rearrangement data.

### Requirements:
R (written and tested with version 3.5.1)  
R packages:  
1. GenomicRanges
2. igraph
3. biovizBase (optional to load cytoband information, see below)

### Usage
1. Prepare your copy number segmentation data as a GRanges object (chromosome, start and end) with an additional column called "cn"
2. Prepare your rearrangement breakpoints as a data.frame with at least four columns named "chr1", "bp1", "chr2", "bp2"
3. Create a GRanges object with UCSC style cytoband information (requires the column "gieStain", with values like "gneg", "gpos25", "acen" etc.). A file for hg19 is provided
4. Run `detect_seismic_amplification(cnv=MY_CNVs, sv=MY_SVS, chrBands=MV_BANDS)`

#### Detailed list of (optional) parameters:

- cnv (required): GRanges object containing a complete copy number profile (all states, normal/gain/loss) and the integer copy number in an extra column called "cn"
- sv (required): data.frame with structural variation breakpoints; requires four columns named "chr1", "bp1", "chr2", "bp2"; can have any kind of additional columns
- chrBands (required): GRanges object, with UCSC style cytoband information (requires the column "gieStain", with values like "gneg", "gpos25", "acen" etc.) it is possible to download it directly into R using the package "biovizBase" and its function getIdeogram("myGenomeVersion", cytobands=TRUE)
- minInternalSVs (optional): numeric, minimum number of internal structural variations required by seismic amplifications (we found 14 to be a reasonable threshold)
- ploidy (optional): numeric, ploidy of the sample; used to distinguish amplifications by copy number
- cnvTol (optional): numeric, maximum distance between a sv breakpoint and a cnv (to be sensitive for inaccuracies in cnv and breakpoint detections)
- noXY (optional): TRUE/FALSE, to exclude chromosomes X and Y from analysis

#### Output:
- list with two items:
  - amplicon: GRanges object with potential seismic amplifications and various annotations
  - svs: data.frame with structural variations associated to called seismic amplifications (this is an annotated subset of the provided svs)

## (2) Rearrangement calling with SoReCa
The C++ source code for the SOmatic REarrangement CAlling (SoReCa) used in the manuscript can be downloaded from [here](http://www.uni-koeln.de/med-fak/soreca/soreca.tgz). All genomic annotation files for hg19,hg38 and mm10 and a database of rearrangement calls from normal genomes for filtering are provided in the download.

### Usage

1. Run `soreca enspan` separately for your tumor and normal alignment file to detect encompassing paired-end reads:
#### enspan parameters:
- -i: input bam-file (required)
- -o: output file (prefix only, required)
- -build: genome build (hg19,hg38,mm10) [hg19]
- -w: maximal insert size for correct pair [600]
- -all: show all results
- -maxr: maximal number of reads to realign [1000]
- -q: mapping quality filter [1] 

2. Run `soreca unsnarl` with the results from step 1 to generate the list of somatic rearrangement breakpoints:
#### unsnarl arameters:
- -inT: input text file for the tumor created by the enspan module (required)
- -inN: input text file for the matched normal created by the enspan module (required)
- -o: output file (prefix, required)
- -adb: add normal to the specified database (prefix: *.Ndb)
- -db: use database of normal genomes (prefix: *.Ndb)
- -seg: copy number seg-file (optional)
- -bam: bam-file of the matched normal (optional)
- -nc: minimal number of clipped bases [10]
- -build: genome build (hg19, hg19u, hg38, mm9, mm10) [hg19]
- -nrn: if this flag is set no read names are shown
- -w: maximal insert size [600]


## (3) Simulation of seismic amplification, chromothripsis and breakage-fusion-bridge amplification
R-scripts for the simulation of seismic amplification, chromothripsis and breakage-fusion-bridge and 
their combinations. Please refer to the header of the script for details on how to run the 
simulations. The script contains the complete set of variables and parameters used in the manuscript.

## (4) Simulation of the double minute (DM) evolutionary model
A C++ program for the simulation of the DM evolutionary model. Run "compile.sh" to compile the script 
and create the main executable "simulate-ring-trajectories". Initial chromothripsis setups for three 
chromosomal regions are provided.

## (5) Microhomology annotation for rearrangement breakpoints
A C++ program for the calling of microhomologies at rearrangement breakpoints. Run "make" from 
within the directory to create the main executable "mh-caller".
