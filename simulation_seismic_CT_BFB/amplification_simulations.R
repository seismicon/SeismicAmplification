# Make sure you have the Bioconductor package "GenomicRanges" and its dependencies installed
library("GenomicRanges")


################################################################################################################################################################
################################################################################################################################################################
# Simulation setup and parameters used in the manuscript "Chromothripsis followed by circular recombination drives oncogene amplification in human cancer"
#
# Functions for the simulation of three mechanisms are provided further below:
# (a) seismic amplification / circular recombination; function seismic_amplification_simulation()
# (b) chromothripsis (CT); function chromothripsis_simulation()
# (c) breakage-fusion-bridge (BFB) amplification; function bfb_simulation()
#
# This is a good starting point to run your own simulations. To do so, modify the objects "initRegions", "preservedRegions", "initSVs", "outputDir", "nBreaks", 
# "nCycles", "pCrossEven" and "nRuns" below according to your specifications.
#
# Before starting the simulation, make sure to run all functions further below once to make them available in your R session.
# The mechamisms can be combined in any order by calling functions for chromothripis, seismic amplification, BFB on one
# regions and svs object and passing it on to the next function.
#
# The "amplicon" data structure that is part of the list returned by the functions is a (large) GRanges object containing all fragments within the
# rearranged/amplified region (or circle in case of seismic amplification); in case of CNVs, fragments appear multiple times (gains) or are removed (deletions);
# the order of fragments in corresponds to the order within the amplicon; the final copy-number profile is computed as a coverage across 
# the "amplicon" object; svs between fragments are kept in a separate data.frame


# initial genomic regions for simulation in three sizes (choose one of a,b,c)
# main effect of the regions is their size, no actual sequence or location info is used during simulation
# preservedRegions are regions that should be preserved during chromothripsis and circular recombination (e.g. oncogenes); they can be empty GRanges()
#
# Size (a): 6mb (1 region)
initRegions = GRanges(seqnames="chr12", IRanges(start=c(54000001), end=c(60000000)), inv=0)
preservedRegions = GRanges(seqnames="chr12", IRanges(start=c(58118980,58131796,58141510,58148881,58156117,58162254,58176372,58191159,58213710,58335324), end=c(58135940,58143994,58149796,58154190,58162769,58166576,58201854,58212487,58240522,58351052)), hgnc_symbol=c("AGAP2","TSPAN31","CDK4","MARCH9","CYP27B1","METTL1","TSFM","AVIL","CTDSP2","XRCC6BP1"))
initSVs = data.frame()
# Size (b): 13mb (2 regions)
initRegions = GRanges(seqnames="chr12", IRanges(start=c(54000001, 68000001), end=c(60000000, 75000000)), inv=c(0,0))
preservedRegions = GRanges(seqnames="chr12", IRanges(start=c(58118980,58131796,58141510,58148881,58156117,58162254,58176372,58191159,58213710,58335324,69201956,69235977,69633317,69742121,69753483,69864129,69979114,70002351,70037140,70132461,70219084), end=c(58135940,58143994,58149796,58154190,58162769,58166576,58201854,58212487,58240522,58351052,69239214,69365350,69668138,69748014,69784576,69973562,69995350,70004942,70093256,70216984,70352877)), hgnc_symbol=c("AGAP2","TSPAN31","CDK4","MARCH9","CYP27B1","METTL1","TSFM","AVIL","CTDSP2","XRCC6BP1","MDM2","CPM","CPSF6","LYZ","YEATS4","FRS2","CCT2","LRRC10","BEST3","RAB3IP","MYRFL"))
initSVs = data.frame("step"="init", "chr1"="chr12", "bp1"=60000000, "chr2"="chr12", "bp2"=68000001, "orientation"="t2h")
# Size (c): 25mb (3 regions)
initRegions = GRanges(seqnames="chr12", IRanges(start=c(54000001, 68000001, 114000001), end=c(60000000, 75000000, 126000000)), inv=c(0,0,0))
preservedRegions = GRanges(seqnames="chr12", IRanges(start=c(58118980,58131796,58141510,58148881,58156117,58162254,58176372,58191159,58213710,58335324,69201956,69235977,69633317,69742121,69753483,69864129,69979114,70002351,70037140,70132461,70219084,121837844,121866900), end=c(58135940,58143994,58149796,58154190,58162769,58166576,58201854,58212487,58240522,58351052,69239214,69365350,69668138,69748014,69784576,69973562,69995350,70004942,70093256,70216984,70352877,121868389,122018920)), hgnc_symbol=c("AGAP2","TSPAN31","CDK4","MARCH9","CYP27B1","METTL1","TSFM","AVIL","CTDSP2","XRCC6BP1","MDM2","CPM","CPSF6","LYZ","YEATS4","FRS2","CCT2","LRRC10","BEST3","RAB3IP","MYRFL","RNF34","KDM2B"))
initSVs = data.frame("step"="init", "chr1"="chr12", "bp1"=c(60000000, 75000000), "chr2"="chr12", "bp2"=c(68000001, 114000001), "orientation"=c("t2h", "t2h"))


# simulate three scenarios:
# (1) chromothripsis followed by BFB
# (2) BFB followed by chromothripsis
# (3) chromothripsis followed by seismic amplification

# chromothripsis parameters: number of random chromothripsis breakpoints
nBreaks = 10 # 20 # 30
# BFB parameters: number of BFB cycles
nCycles = 3 # 2 # 4 # 5 # 7 # 10

## scenario (1)
amp = chromothripsis_simulation(regions=initRegions, preservedRegions=preservedRegions, svs=initSVs, nBreaks=nBreaks, invP=0.5, delP=0.5)
amp = bfb_simulation(regions=amp$amplicon, svs=amp$svs, nCycles=nCycles)
          
## scenario (2)
amp = bfb_simulation(regions=initRegions, svs=initSVs, nCycles=nCycles)
amp = chromothripsis_simulation(regions=amp$amplicon, preservedRegions=preservedRegions, svs=amp$svs, nBreaks=nBreaks, invP=0.5, delP=0.5)

# Scenario (3)
# chromothripsis parameters: number of random chromothripsis breakpoints
nBreaks = 10 # 20 # 30
# number of circular recombination cycles
nCycles = 10 # 20 # 30
# probability for two crossing overs during circular recombination
pCrossEven = 0.9 # 0.25 # 0.5 # 0.75 #1

amp = chromothripsis_simulation(regions=initRegions, preservedRegions=preservedRegions, svs=initSVs, nBreaks=nBreaks, invP=0.5, delP=0.5)
amp = seismic_amplification_simulation(amp$amplicon, amp$svs, preservedRegions, nCycles=nCycles, maxCN=1e4, maxCircleSize=8e8, pCrossEven=pCrossEven, N=16, R=4, p0=0.3)


################################################################################################################################################################
################################################################################################################################################################

# Main function for the computational simulation of seismic amplification 
# as used in the manuscipt "Chromothripsis followed by circular recombination drives oncogene amplification in human cancer"
#
# In brief, the simulation consists of the repetition of 2 steps:
# (1) Circular recombination
# (2) DNA circle selection
# For more information on the model and methods, please refer to the supplementary material
#
# Parameters (the first four are required):
# regions: GRanges object with genomic coordinates that will be part of the seismic amplification 
# svs: data.frame with svs from previous simulation steps (usually chromothripsis), or an empty data.frame
# preservedRegions: GRanges object with genomic coordinates of regions (e.g. oncogenes) that provide a selective advantage and must be preserved during circular recombination
# nBreakds: numeric, number of chromothripsis breakpoints
# nCycles: numeric, maximum number of recombination cycles
# maxCN: numeric, maximum allowed copy number, recombination stops if amplification exceeds this copy number
# maxCircleSize: numeric, maximum size of circles during recombination; larger circles are rejected
# pCrossEven: numeric, probability for two crossing overs
# N: numeric, maximum copy number above which no further advantage can be gained; used during circle selection
# R: numeric, decay of additional selective advantage; used during circle selection
# p0: numeric, initial probability for selection model
#
# Returns:
# list with three items:
# regions: fragments of seismic amplicon as a GRanges object
# cnvs: copy number profile as GRanges object
# svs: sv breakpoints as data.frame with columns "step", "chr1", "bp1", "chr2", "bp2", "orientation"

seismic_amplification_simulation <- function(regions, svs, preservedRegions, nCycles, maxCN=150, maxCircleSize=1e9, pCrossEven=0.9, N=16, R=4, p0=0.3){
  
  message("Running seismic amplification simulation")
  
  amp = regions
  
  # (1) Circular recombination
  nCycle = 1
  prevScore = 1
  currentMaxCN = 1
  
  # recombination will stop at maximum number of cycles or if copy number exceeds given theshold
  while((nCycle <= nCycles) && (currentMaxCN <= maxCN)){
    
    message("Cycle ", nCycle)
    
    # randomly choose 1 or 2 crossing-overs with given probability
    nCross = sample(1:2, 1, prob=c(1-pCrossEven,pCrossEven))
    
    # replicate circles; modify both circles separately, with sv breakpoints for each circle
    circle1 = amp
    circle2 = amp
    circle1_sv1 = GRanges()
    circle1_sv2 = GRanges()
    circle2_sv1 = GRanges()
    circle2_sv2 = GRanges()
    
    # simulate crossing-overs
    for(i in 1:nCross){
      
      # distance between crossing-over bps is max. 10 percent of circle size
      maxSizeSV = round(sum(width(circle1))/10)
      # select crossing-over bp randomly within fragments of at least 3bp
      idx = which(width(circle1) > 3)
      r1 = ifelse(length(idx) == 1, idx, sample(idx, 1, prob=width(circle1[idx])))
      bp = sample(start(circle1)[r1]:(end(circle1)[r1]-2), 1)
      sizeSV = sample(0:maxSizeSV, 1)
      # in case of two crossing-overs, look for the same bp in the second circle; in case of multiple matches, due to previous duplications, choose randomly one of the matches
      # select which one of the circles get the deletion (whichDel), other one gets the duplication
      if(i %% 2 == 0){
        bp1 = GRanges(seqnames(circle1)[r1], IRanges(bp, bp))
        r2 = findOverlaps(bp1, circle2, select="arbitrary")
        whichDel = sample(1:2, 1)
      }else{
        r2 = r1
        whichDel = 1
      }
      
      if(whichDel == 1){
        # implement deletion into circle1; capture deleted fragment for insertion/duplication into other circle
        del = deletion(circle1, r1, bp, sizeSV)
        insert = del[[2]]
        circle1 = del[[1]]
        if(del[[3]] > r1){
          # deletion does not cross the beginning of the circle
          circle1_sv1 = c(circle1_sv1, addSV(circle1, r1, r1+1, 1))
          circle1_sv2 = c(circle1_sv2, addSV(circle1, r1, r1+1, 2))
        }else{
          # deletion crosses the beginning of the circle
          circle1_sv1 = c(circle1_sv1, addSV(circle1, length(circle1), 1, 1))
          circle1_sv2 = c(circle1_sv2, addSV(circle1, length(circle1), 1, 2))
        }
        # implement duplication into circle2 using the deleted fragment from above
        circle2 = insertion(circle2, r2, bp, sizeSV, insert)
        insert$inv = ifelse(insert$inv == 0, 1, 0)  # reverse insertion elements to mimic a h2t duplication rearrangement instead of t2h
        circle2_sv1 = c(circle2_sv1, addSV(insert, 1, length(insert), 1))
        circle2_sv2 = c(circle2_sv2, addSV(insert, 1, length(insert), 2))
      }else{
        # implement deletion into circle2; capture deleted fragment for insertion/duplication into other circle
        del = deletion(circle2, r2, bp, sizeSV)
        insert = del[[2]]
        circle2 = del[[1]]
        if(del[[3]] > r2){
          # deletion does not cross the beginning of the circle
          circle2_sv1 =c(circle2_sv1,  addSV(circle2, r2, r2+1, 1))
          circle2_sv2 =c(circle2_sv2,  addSV(circle2, r2, r2+1, 2))
        }else{
          # deletion crosses the beginning of the circle
          circle2_sv1 = c(circle2_sv1, addSV(circle2, length(circle2), 1, 1))
          circle2_sv2 = c(circle2_sv2, addSV(circle2, length(circle2), 1, 2))
        }
        # implement duplication into circle1 using the deleted fragment from above
        circle1 = insertion(circle1, r1, bp, sizeSV, insert)
        insert$inv = ifelse(insert$inv == 0, 1, 0)  # reverse insertion elements to mimic a h2t duplication rearrangement instead of t2h
        circle1_sv1 = c(circle1_sv1, addSV(insert, 1, length(insert), 1))
        circle1_sv2 = c(circle1_sv2, addSV(insert, 1, length(insert), 2))
      }
      # in case of just one crossing-over, combine both circles into one large circle
      if(nCross == 1){
        circle1 = c(circle1, circle2)
        circle1_sv1 = c(circle1_sv1, circle2_sv1)
        circle1_sv2 = c(circle1_sv2, circle2_sv2)
        circle2 = GRanges()
        circle2_sv1 = GRanges()
        circle2_sv2 = GRanges()
      }
      
    }
    
    # (2) DNA circle selection
    # a) if no parameter preservedRegions with regions of selective advantage is given, select the larger circle
    # b) otherwise compute a 'fitness score' based on the model of Garsed et al. "The architecture and evolution of cancer neochromosomes", Cancer Cell 26, 653-667 (2014) and select the circle having a higher score
    # please see the Garsed et al. paper, especially the supplementary material, for details on the score computation and whole selection model
    score1 = 1
    score2 = 1
    rejectCircle1 = FALSE
    rejectCircle2 = FALSE
    
    if(length(preservedRegions) > 0){
      
      # fitness score is a product of contributions from every region in preservedRegions depending on its (mean) copy number in a circle
      for(j in 1:length(preservedRegions)){
        # circle1
        cn1j = subsetByOverlaps(circle1, preservedRegions[j])
        if(length(cn1j) > 0){
          # compute mean copynumber of preservedRegion j
          start(cn1j)[which(start(cn1j) < start(preservedRegions[j]))] = start(preservedRegions[j])
          end(cn1j)[which(end(cn1j) > end(preservedRegions[j]))] = end(preservedRegions[j])
          cn1j = coverage(cn1j)[[as.character(seqnames(preservedRegions)[j])]]
          cn1j = cn1j[cn1j > 0]
          cn1j = round(sum(runValue(cn1j)) / nrun(cn1j))
        }else{
          # discard circle1 if it does not contain all given preservedRegions
          rejectCircle1 = TRUE
          cn1j = 0
        }
        # circle2
        cn2j = subsetByOverlaps(circle2, preservedRegions[j])
        if(length(cn2j) > 0){
          # compute mean copynumber of preservedRegion j
          start(cn2j)[which(start(cn2j) < start(preservedRegions[j]))] = start(preservedRegions[j])
          end(cn2j)[which(end(cn2j) > end(preservedRegions[j]))] = end(preservedRegions[j])
          cn2j = coverage(cn2j)[[as.character(seqnames(preservedRegions)[j])]]
          cn2j = cn2j[cn2j > 0]
          cn2j = round(sum(runValue(cn2j)) / nrun(cn2j))
        }else{
          # discard circle2 if it does not contain all given preservedRegions
          rejectCircle2 = TRUE
          cn2j = 0
        }
        # update score and proceed to next region in recRegion
        score1 = selectionScore(cn1j, R, N) * score1
        score2 = selectionScore(cn2j, R, N) * score2
      }
      # reject circles larger than maxCircleSize
      rejectCircle1 = ifelse(sum(as.numeric(width(circle1))) > maxCircleSize, TRUE, rejectCircle1)
      rejectCircle2 = ifelse(sum(as.numeric(width(circle2))) > maxCircleSize, TRUE, rejectCircle2)
      # compute a selection probability for each circle from fitness scores; probability uses the probability of the circles from the previous cycle
      # set selection probability to 0 for rejected circles
      prevProb = selectionProb("p", prevScore, score1, score2, p0)
      prob1 = ifelse(rejectCircle1, 0, selectionProb("d", prevScore, score1, score2, p0))
      prob2 = ifelse(rejectCircle2, 0, selectionProb("d", prevScore, score2, score1, p0))
      
      # if both selection probabilites are lower compared to the previous cycle, forget everything and skip this cycle (i.e. currentMaxCN, amp and prevScore stay the same)
      # otherwise, select circle with maximum selection probability
      maxProb = max(prevProb, prob1, prob2)
      if(prob1 == maxProb || prob2 == maxProb){
        # select circle1
        if(prob1 > prob2){
          prevScore = score1
          amp = circle1
          sv1 = circle1_sv1
          sv2 = circle1_sv2
        }
        # select circle2
        if(prob1 < prob2){
          prevScore = score2
          amp = circle2
          sv1 = circle2_sv1
          sv2 = circle2_sv2
        }
        # if both probabilites are equal, select the largest circle
        if(prob1 == prob2){
          if(sum(as.numeric(width(circle1))) > sum(as.numeric(width(circle2)))){
            prevScore = score1
            amp = circle1
            sv1 = circle1_sv1
            sv2 = circle1_sv2
          }else{
            prevScore = score2
            amp = circle2
            sv1 = circle2_sv1
            sv2 = circle2_sv2
          }
        }
      }
    }else{
      # just one crossing-over was chosen or no regions for selection model were given, so select the largest circle
      if(sum(as.numeric(width(circle1))) > sum(as.numeric(width(circle2)))){
        amp = circle1
        sv1 = circle1_sv1
        sv2 = circle1_sv2
      }else{
        amp = circle2
        sv1 = circle2_sv1
        sv2 = circle2_sv2
      }
    }
    if(!(rejectCircle1 & rejectCircle2)){
      # remove svs not overlapping the chosen circle
      idx = intersect(queryHits(findOverlaps(sv1, amp)), queryHits(findOverlaps(sv2, amp)))
      sv1 = sv1[idx]
      sv2 = sv2[idx]
      
      # format and save svs
      svs_seismic = data.frame("step"="seismic", "chr1"=seqnames(sv1), "bp1"=start(sv1), "chr2"=seqnames(sv2), "bp2"=start(sv2), "orientation"=sv1$Orientation)
      svs = rbind(svs, svs_seismic)
      
      nCycle = nCycle + 1
    }   
    
  }
  
  cov = coverage(amp)
  cnvs = rle2ranges(cov)

  message("Done")
  return(list("amplicon"=amp, "cnvs"=cnvs, "svs"=svs))
}

################################################################################################################################################################
################################################################################################################################################################

# Function for the computational simulation of chromothripsis
# as used in the manuscipt "Chromothripsis followed by circular recombination drives oncogene amplification in human cancer"
#
# In brief, the simulation splits a GRanges object at random breakpoints and deletes a fraction of ranges while reorganizing the rest
# For more information on the model and methods, please refer to the supplementary material
#
# Parameters (the first four are required):
# regions: GRanges object with genomic coordinates that will be affected by the chromothripsis simulation
# preservedRegions: GRanges object with genomic coordinates of regions (e.g. oncogenes) that provide a selective advantage and must be preserved during chromothripsis
# svs: data.frame with sv breakpoints; if chromothripsis is the first/only simulation, then this is an empty data.frame, otherwise it is the output of a previous simulation (seismic or BFB)
# nBreaks: numeric, number of random chromothripsis breakpoints
# invP: numeric, probability for a chromothripsis fragment to be inverted
# delP: numeric, probability for a chromothripsis fragment to be deleted
#
# Returns:
# list with three items:
# regions: reorganized fragments as a GRanges object including information if fragment is inverted or not
# cnvs: copy number profile as GRanges object
# svs: sv breakpoints as data.frame with columns "step", "chr1", "bp1", "chr2", "bp2", "orientation"

chromothripsis_simulation <- function(regions, preservedRegions, svs, nBreaks, invP=0.5, delP=0.2){
  
  message("Running chromothripsis simulation")
  
  # initialize inversion columns if not found
  if(is.null(mcols(regions)$inv)){
    mcols(regions) = NULL
    regions$inv = 0
  }
  amp_list = list(regions)
  
  # chromothripsis, breakage and reorganization of DNA fragments
  if(nBreaks > 0){
    
    for(n in 1:nBreaks){
      i = sample(1:length(amp_list), 1, prob=sapply(amp_list, function(x) sum(width(x))))
      amp = amp_list[[i]]
      idx = which(width(amp) > 3)
      # choose region with probability proportional to its size
      r = sample(idx, 1, prob=width(amp[idx]))
      # choose breakpoint randomly within region r
      bp = sample((start(amp)[tail(r,1)]+1):(end(amp)[tail(r,1)]-1),1)
      # break amplicon at breakpoint bp
      j = length(amp_list)
      # breakpoint before last fragment
      if(r < length(amp)){
        amp_list[[i]] = amp[1:r]
        end(amp_list[[i]])[r] = bp
        amp_list[[j+1]] = c(GRanges(seqnames=seqnames(amp)[r], ranges=IRanges(bp+1, end(amp)[r]), inv=amp$inv[r]), amp[(r+1):length(amp)])
      # breakpoint within last fragment
      }else{
        amp_list[[i]] = amp[1:r]
        end(amp_list[[i]])[r] = bp
        amp_list[[j+1]] = c(GRanges(seqnames=seqnames(amp)[r], ranges=IRanges(bp+1, end(amp)[r]), inv=amp$inv[r]))
      }
      
    }
    
  }
  
  # some framgents get lost with given probability delP
  # take care not to lose fragments overlapping with regions from preservedRegions
  idx = which(sapply(amp_list, function(x) all(countOverlaps(x, preservedRegions)==0)))
  if(length(idx) > 0){
    idx = idx[sample(c(FALSE,TRUE), length(idx), prob=c(1-delP, delP), replace=TRUE)]
  }
  if(length(idx) > 0){
    amp_list = amp_list[-(idx)]
  }
  # reorganization of fragments in random order
  idx = sample(1:length(amp_list),length(amp_list))
  amp_list = amp_list[idx]
  # invert segments with given probability invP
  amp_list = lapply(amp_list, function(x){
    if(sample(c(FALSE,TRUE), 1, prob=c(1-invP, invP))){
      x = rev(x)
      x$inv = abs(x$inv - 1)
    }
    return(x)})
  
  # add connections between chromothripsis fragments to svs
  sv1 = GRanges()
  sv2 = GRanges()
  for(i in 1:(length(amp_list)-1)){
    ct_bp1 = amp_list[[i]]
    ct_bp1 = ct_bp1[length(ct_bp1)]
    ct_bp2 = amp_list[[i+1]]
    ct_bp2 = ct_bp2[1]
    ct_bp = c(ct_bp1, ct_bp2)
    ct_bp1 = addSV(ct_bp, 1, 2, 1)
    ct_bp2 = addSV(ct_bp, 1, 2, 2)
    # sort breakpoints on the same chromosome
    if(as.character(seqnames(ct_bp1)) == as.character(seqnames(ct_bp2)) & start(ct_bp1) > start(ct_bp2)){
      s = start(ct_bp1)
      ranges(ct_bp1) = IRanges(start(ct_bp2), start(ct_bp2))
      ranges(ct_bp2) = IRanges(s, s)
      ct_bp1$Orientation = reverse(ct_bp1$Orientation)
      ct_bp2$Orientation = reverse(ct_bp2$Orientation)
    }
    
    sv1 = c(sv1, ct_bp1)
    sv2 = c(sv2, ct_bp2)
  }
  amp = do.call(c, amp_list)
  
  # calculate copy number profile
  cov = coverage(amp)
  cnvs = rle2ranges(cov)

  # format svs
  svs = rbind(svs, data.frame("step"="chromothripsis", "chr1"=seqnames(sv1), "bp1"=start(sv1), "chr2"=seqnames(sv2), "bp2"=start(sv2), "orientation"=sv1$Orientation))
  
  message("Done")
  
  return(list("amplicon"=amp, "cnvs"=cnvs, "svs"=svs))  
}


###############################################################################################################################################################
################################################################################################################################################################

# Function for the computational simulation of the breakage-fusion-bridge (BFB) amplification mechamism
# as used in the manuscipt "Chromothripsis followed by circular recombination drives oncogene amplification in human cancer"
#
# In brief, the simulation splits a GRanges object repeatedly at a random breakpoint, deletes the latter part and duplicates the former part in an inverted fashion (fold-back inversion)
# For more information on the model and methods, please refer to the supplementary material
#
# Parameters (all three are required):
# regions: GRanges object with genomic coordinates that will be part of the BFB simulation 
# svs: data.frame with sv breakpoints; if bfb is the first/only simulation, then this is an empty data.frame, otherwise it is the output of a previous simulation (seismic or chromothripsis)
# nCycles: numeric, number of BFB cycles
#
# Returns:
# list with three items:
# regions: fragments of the BFB amplicon as a GRanges object
# cnvs: copy number profile as GRanges object
# svs: sv breakpoints as data.frame with columns "step", "chr1", "bp1", "chr2", "bp2", "orientation"

bfb_simulation <- function(regions, svs, nCycles){
 
  message("Running breakage-fusion-bridge simulation")
  
  amp = regions
  
  # initialize inversion columns if not found
  if(is.null(mcols(amp)$inv)){
    mcols(amp) = NULL
    amp$inv = 0
  }
  
  cycle = 0
  
  while(cycle <= nCycles){
    
    if(cycle > 0){
      message("Cycle ", cycle)
    }
    
    # select the first fold-back bp (telomere loss) at the end of the given region, so that the whole region is affected by bfb
    if(cycle == 0){
      r = length(amp)
      bp = ifelse(amp[r]$inv == 0, end(amp[r]), start(amp[r]))
    # after the first fold-back bp, select the bp randomly
    }else{
      # select random fold-back breakpoint (position of loss of telomeric region or fusion sister chromatids)
      idx = which(width(amp) > 3)
      # set breakpoint in previously duplicated regions
      if(cycle > 0){
        idx = idx[(length(idx)/2+1):length(idx)]
      }
      r = ifelse(length(idx) == 1, idx, sample(idx, 1, prob=width(amp[idx])))
      bp = sample((start(amp)[r]+1):(end(amp)[r]-1), 1)
    }
    
    # apply breakpoint to regions:
    # split bp-region at bp and remove all regions up-/downstream of regions
    # bfb on q-arm: remove downstream regions
    if(amp[r]$inv == 0){
      end(amp[r]) = bp
      orientation = "t2t"
    }else{
      start(amp[r]) = bp
      orientation = "h2h"
    }
    if(r < length(amp)){
      amp = amp[-((r+1):length(amp))]
    }  
    
    # save current fold-back inversion
    sv1 = GRanges(seqnames(amp[r]), IRanges(bp-1, bp-1), "Orientation"=orientation)
    sv2 = GRanges(seqnames(amp[r]), IRanges(bp+1, bp+1), "Orientation"=orientation)
    svs = rbind(svs, data.frame("step"="BFB", "chr1"=seqnames(sv1), "bp1"=start(sv1), "chr2"=seqnames(sv2), "bp2"=start(sv2), "orientation"=sv1$Orientation))
    
    ## replication and fusion
    if(cycle < nCycles){
      # add inverted duplicated regions
      amp_fb = rev(amp)
      amp_fb$inv = abs(1-amp_fb$inv)
      amp = c(amp, amp_fb)
    }
    cycle = cycle + 1
    
  }
  
  # calculate copy number profile
  cov = coverage(amp)
  cnvs = rle2ranges(cov)

  message("Done")
  
  return(list("amplicon"=amp, "cnvs"=cnvs, "svs"=svs))
  
}


################################################################################################################################################################
################################################################################################################################################################
# Other smaller functions used by the simulations above


# Function used by seismic_amplification_simulation() to compute fitness score for a region
#
# Parameters:
# cn: numeric, copy number of region
# R: numeric, decay of additional selective advantage
# N: numeric, maximum copy number above which no further advantage can be gained
#
# Returns:
# fitness score, numeric
selectionScore <- function(cn, R, N){
  n = min(cn, N)
  sj = (R*N+1)/(N-1)
  if(n > 1){
    for(j in 2:n){
      sj = ((R*(N+1)-(R-1)*j) / (N-1)) * sj
    }
  }
  return(sj)
}

# Function used by seismic_amplification_simulation() to compute selection probability
#
# Parameters:
# nc: either "d" for daughter cycle, i.e. current cycle, or "p" for parent, i.e. previous cylce
# sp: numeric, parent/previous fitness score
# sd1: numeric, current daughter1/circle1 fitness score
# sd2: numeric, current daughter2/circle2 fitness score
# p0: numeric, initial/minimum probability
#
# Returns:
# selection probability, numeric
selectionProb <- function(nc="d", sp, sd1, sd2, p0=0.3){
  p = 0
  if(any(!is.na(c(sp,sd1,sd2)))){
    if(nc == "p"){
      p = p0 + 0.5*(1-3*p0)*(sp/(sd1+sd2+sp))
    }
    if(nc == "d"){
      p = p0 + 0.5*(1-3*p0)*(sd1/(sd1+sd2+sp)+sd1/(sd1+sd2))
    }
  }
  return(p)
}

# Function used by seismic_amplification_simulation() to delete a fragment from a circle
#
# Parameters:
# amp: GRanges object with circle regions
# r: numeric, row index for insertion into amp
# bp: numeric, position within row r
# sizeSV: numeric, deletion length
#
# Returns:
# GRanges object with modified circle regions
# GRanges with deleted region
# numeric row index in amp where deletion ends
deletion <- function(amp, r, bp, sizeSV){
  if(sizeSV > 0){
    # start of deletion in row r at position bp; split amp accordingly
    amp = amp[c(1:r,r:length(amp))]
    end(amp[r]) = bp
    start(amp[r+1]) = bp+1
    # go sizeSV positions further into amp to find the end row i of the deletion; remember rows that are completely deleted in vector idx
    i = r+1
    idx = c()
    while(width(amp)[i] < sizeSV){
      idx = c(idx, i)
      sizeSV = sizeSV - width(amp)[i]
      i = ifelse(i == length(amp), 1, i + 1)
    }
    # get deleted region
    insert = amp[c(idx, i)]
    end(insert)[length(insert)] = start(insert)[length(insert)]  + sizeSV - 2
    start(amp)[i] = start(amp)[i] + sizeSV - 1
    # delete region
    if(length(idx) > 0){
      amp = amp[-idx]
    }
  }else(
    insert = GRanges()
  )
  return(list(amp, insert, i))
}

# Function used by seismic_amplification_simulation() to insert a fragment into a circle
#
# Parameters:
# amp: GRanges object with circle regions
# r: numeric, row index for insertion into amp
# bp: numeric, position within row r
# insert: GRanges with region to be inserted into amp
#
# Returns:
# GRanges object with modified circle regions
insertion <- function(amp, r, bp, sizeSV, insert){
  if(sizeSV > 0){
    # split circle at row r and position bp and insert fragment
    amp = amp[c(1:r,r:length(amp))]
    end(amp[r]) = bp
    start(amp[r+1]) = bp+1
    amp = c(amp[1:r], insert, amp[(r+1):length(amp)])
  }
  return(amp)  
}

# Function used by seismic_amplification_simulation() and chromothripsis_simulation() to create an SV beakpoint object from given position
#
# Parameters:
# amp: GRanges object with circle regions
# i1: numeric, row index of amp where sv starts
# i2: numeric, row index of amp where sv ends
# bp: numeric, either 1 or 2, for breakpoint 1 or 2 of an SV
#
# Returns:
# GRanges object with one SV breakpoint
addSV <- function(amp, i1, i2, bp){
  
  # set orientation of SV to h2t, t2h, h2h, or t2t; h=head, t=tail
  orientation1 = ifelse(amp$inv[i1] == 0, "t", "h")
  orientation2 = ifelse(amp$inv[i2] == 0, "h", "t")
  orientation = paste(orientation1, orientation2, sep="2")
  # either first breakpoint
  if(bp == 1){
    if(amp$inv[i1] == 0){
      sv =GRanges(seqnames=seqnames(amp)[i1], IRanges(end(amp)[i1], end(amp)[i1]), "Orientation"=orientation)
    }else{
      sv =GRanges(seqnames=seqnames(amp)[i1], IRanges(start(amp)[i1], start(amp)[i1]), "Orientation"=orientation)
    }
    # or second breakpoint
  }else{
    if(amp$inv[i2] == 0){
      sv =GRanges(seqnames=seqnames(amp)[i2], IRanges(start(amp)[i2], start(amp)[i2]), "Orientation"=orientation)
    }else{
      sv =GRanges(seqnames=seqnames(amp)[i2], IRanges(end(amp)[i2], end(amp)[i2]), "Orientation"=orientation)
    }
  }
  return(sv)
}

# Function used by all simulations to convert a coverage rle object into a GRanges object
#
# Parameters:
# rle: rle object from GenomicRanges coverage() function
#
# Returns:
# GRanges object based on rle object
rle2ranges <- function(rle){
  
  ranges = GRanges()
  for(chr in names(rle)){
    r = rle[[chr]]
    for(i in 1:(nrun(r))){
      if(i == 1){
        ranges = suppressWarnings(c(ranges, GRanges(seqnames=chr, IRanges(1, runLength(r)[i]), cn = runValue(r)[i])))
      }else{
        ranges = suppressWarnings(c(ranges, GRanges(seqnames=chr, IRanges(tail(end(ranges),1)+1, tail(end(ranges),1)+runLength(r)[i]), cn = runValue(r)[i])))
      }
    }
  }
  ranges$cn = ranges$cn + 1
  ranges = ranges[-1]
  return(ranges)
}