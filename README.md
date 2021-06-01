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
4. Run **detect_seismic_amplification(cnv=MY_CNVs, sv=MY_SVS, chrBands=MV_BANDS)**

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


## (2) Simulation of seismic amplification, chromothripsis and breakage-fusion-bridge amplification
R-scripts for the simulation of seismic amplification, chromothripsis and breakage-fusion-bridge and 
their combinations. Please refer to the header of the script for details on how to run the 
simulations. The script contains the complete set of variables and parameters used in the manuscript.

## (3) Simulation of the double minute (DM) evolutionary model
A C++ program for the simulation of the DM evolutionary model. Run "compile.sh" to compile the script 
and create the main executable "simulate-ring-trajectories". Initial chromothripsis setups for three 
chromosomal regions are provided.

## (4) Microhomology annotation for rearrangement breakpoints
A C++ program for the calling of microhomologies at rearrangement breakpoints. Run "make" from 
within the directory to create the main executable "mh-caller".
