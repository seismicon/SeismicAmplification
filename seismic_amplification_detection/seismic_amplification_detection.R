################################################################################################################################################################
################################################################################################################################################################
#
# Main function for the detection of seismic amplification 
# as used in the manuscipt "Chromothripsis followed by circular recombination drives oncogene amplification in human cancer"
#
# Before starting, make sure the following packages and their dependencies are installed:
# GenomicRanges
# igraph
# biovizBase (optional to load cytoband information, see below)
#
# In brief, the seismic amplification detection consists of 4 steps:
# (1) Find amplified regions from given copy number profile
# (2) Compute regions spanned by structural variations
# (3) Overlap amplified regions from (1) with  regions from (2) to retrieve potential seismic amplifications
# (4) Filter seismic amplifications by the number of internal structural variations
# For more detailed information on the detection method, please refer to the supplementary material
#
# Parameters (first three are required):
# cnv: GRanges object containing a complete copy number profile (all states, normal/gain/loss) and the integer copy number in an extra column called "cn"
# sv: data.frame with structural variation breakpoints; requires four columns named "chr1", "bp1", "chr2", "bp2"; can have any kind of additional columns
# chrBands: GRanges object, with UCSC style cytoband information (requires the column "gieStain", with values like "gneg", "gpos25", "acen" etc.)
#           it is possible to download it directly into R using the package "biovizBase" and its function getIdeogram("myGenomeVersion", cytobands=TRUE)
# minInternalSVs: numeric, minimum number of internal structural variations required by seismic amplifications (we found 14 to be a reasonable threshold)
# ploidy: numeric, ploidy of the sample; used to distinguish amplifications by copy number
# cnvTol: numeric, maximum distance between a sv breakpoint and a cnv (to be sensitive for inaccuracies in cnv and breakpoint detections)
# noXY: TRUE/FALSE, to exclude chromosomes X and Y from analysis
#
# Returns:
# list with two items:
# amplicon: GRanges object with potential seismic amplifications and various annotations
# svs: data.frame with structural variations associated to called seismic amplifications (this is an annotated subset of the provided svs)

detect_seismic_amplification <- function(cnv, sv, chrBands, minInternalSVs=14, ploidy=2, cnvTol=5000, noXY=TRUE){
  
  require(GenomicRanges)
  require(igraph)
  
  message("Detecting seismic amplification")
  
  # prepare genomic annotations for centromeres, telomeres and chromosome arms
  # centromeres and telomeres are used to filter false-positive cnvs (incl. surrounding regions up to 1mb)
  chrs = seqlevels(chrBands)
  centromeres = reduce(chrBands[chrBands$gieStain == "acen",])
  chrLengths = aggregate(width(chrBands), by=list("chr"=as.character(seqnames(chrBands))), FUN=sum)
  chrLengths = setNames(chrLengths$x, chrLengths$chr)[chrs]
  telomeres_p = GRanges(chrs, IRanges(rep(1, length(chrs)), rep(1e6, length(chrs))))
  telomeres_q = GRanges(chrs, IRanges(chrLengths-1e6, chrLengths))
  chrArms = setdiff(chrBands, centromeres)
  chrArms$arm = rep(c("p","q"), length(chrs))
  chrArms$name = paste0(seqnames(chrArms), chrArms$arm)
  regionFilter = c(telomeres_p, centromeres + 1e5, telomeres_q)
  
  # (1) Find amplified regions from given copy number profile
  # classify copy numbers regions into amplifications
  # amplification threshold (CN >= 5 or CN >= 9) depends on ploidy (2 or > 2)
  cnv$type = "non-amp"
  if(ploidy > 2){
    cnv$type[cnv$cn >= 9] = "amp"
  }else{
    cnv$type[cnv$cn >= 5] = "amp"
  }
  # normalize copy number to ploidy 2
  cnv$cn = (cnv$cn - (ploidy -2))
  
  # optional: do not detect seismic amplification on chromosomes X and Y
  if(noXY){
    cnv = cnv[!(seqnames(cnv) %in% c("chrX","chrY"))]
  }
  
  # remove cnvs overlapping with telomeres or centromeres
  cnv = filter_cnv(cnv, regionFilter)
  # add chromosome arm to chromosome names
  cnv = chr2chrArm(cnv, chrs, chrArms)
  cnv_segments = cnv
  cnv_segments_amp = cnv[cnv$type == "amp"]
  cnv_amp = reduce(cnv_segments_amp)
  
  # (2) Compute regions spanned by structural variants
  if(length(cnv_amp) > 0){
    
    # split svs into individual breakpoint objects (easier to compute sv overlaps with cnvs)
    bp1 = GRanges(sv$chr1, IRanges(sv$bp1, sv$bp1))
    bp2 = GRanges(sv$chr2, IRanges(sv$bp2, sv$bp2))
    
    # remove svs overlapping with telomeres or centromeres
    idx = ((countOverlaps(bp1, regionFilter) == 0) & (countOverlaps(bp2, regionFilter) == 0))
    bp1 = bp1[idx]
    bp2 = bp2[idx]
    sv = sv[idx,]
    
    # add chromosome arm to chromosome names
    bp1 = chr2chrArm(bp1, chrs, chrArms)
    bp2 = chr2chrArm(bp2, chrs, chrArms)
    
    # use only svs overlapping with amplified regions (+/- given tolerance "cnvTol")
    idx = (countOverlaps(bp1, cnv_amp+cnvTol) > 0) & (countOverlaps(bp2, cnv_amp+cnvTol) > 0)
    bp1 = bp1[idx]
    bp2 = bp2[idx]
    sv = sv[idx,]
    
    ## remove duplicated sv entries (both breakpoints identical or within 2 bp)
    if(length(bp1) >= 2){
      dup = rep(FALSE, length(bp1))
      for(i in 2:length(bp1)){
        if((abs(start(bp1[i]) - start(bp1[i-1])) <= 2) && abs(start(bp2[i]) - start(bp2[i-1])) <= 2){
          dup[i] = TRUE
        }
      }
      bp1 = bp1[!dup]
      bp2 = bp2[!dup]
      sv = sv[!dup,]
    }
    
    # seismic regions are computed from the union of sv spans (the interval between start and end of an sv)
    # (a) sv spans of intra-chromosomal svs
    idx = which(seqnames(bp1) == seqnames(bp2))
    svRangeIntra = GRanges(seqnames(bp1[idx]), IRanges(start(bp1)[idx], end(bp2[idx])))
    # (b) sv spans of inter-chromosomal svs
    # compute the span as the region between one breakpoint and the next neightbouring breakpoint of another sv on the same amplified region
    idx = which(seqnames(bp1) != seqnames(bp2))
    svNeighbour = c(bp1, bp2)
    rangeChr = rangeStart = rangeEnd = c()
    for(i in idx){
      amp_bp1 = subsetByOverlaps(cnv_amp+cnvTol, bp1[i])[1]
      amp_bp2 = subsetByOverlaps(cnv_amp+cnvTol, bp2[i])[1]
      bp1Neighbour = svNeighbour[start(svNeighbour) != start(bp1[i]) & (countOverlaps(svNeighbour, amp_bp1) > 0)]
      bp2Neighbour = svNeighbour[start(svNeighbour) != start(bp2[i]) & (countOverlaps(svNeighbour, amp_bp2) > 0)]
      if(length(bp1Neighbour) > 0 & length(bp2Neighbour) > 0){
        bp1Neighbour = bp1Neighbour[which.min(abs(start(bp1Neighbour)-start(bp1[i])))]
        rangeChr = c(rangeChr, as.character(seqnames(bp1Neighbour)))
        rangeStart = c(rangeStart, min(start(bp1[i]), start(bp1Neighbour)))
        rangeEnd = c(rangeEnd, max(start(bp1[i]), start(bp1Neighbour)))
        bp2Neighbour = bp2Neighbour[which.min(abs(start(bp2Neighbour)-start(bp2[i])))]
        rangeChr = c(rangeChr, as.character(seqnames(bp2Neighbour)))
        rangeStart = c(rangeStart, min(start(bp2[i]), start(bp2Neighbour)))
        rangeEnd = c(rangeEnd, max(start(bp2[i]), start(bp2Neighbour)))
      }
    }
    svRangeInter = GRanges(rangeChr, IRanges(rangeStart, rangeEnd))
    # combine intra- and inter-chromosomal sv-ranges
    svRange = c(svRangeIntra, svRangeInter)
    svRange = reduce(svRange)
    
    # (3) Overlap amplified regions from (1) with  regions from (2) to retrieve potential seismic amplifications
    # split amplified cnv segments at sv breakpoints (some cnvs may go way beyond breakpoints, leading to overestimations of amplicon sizes)
    amps = split_regions(cnv_segments_amp, svRange)
    # define amplicon regions as continuous sets of amplified segments within sv ranges
    amps = subsetByOverlaps(amps, svRange)
    amps = reduce(amps)
    # remove amplicons below the given cnv tolerance level
    amps = amps[width(amps) > (2*cnvTol)]
    
    if(length(amps) > 0){
      
      # amplified regions that are interconnected by svs get the same id (hence they belong to the same amplicon)
      ids = get_connected_regions(amps, bp1, bp2, cnvTol)
      amps$id = ids
      
      ## add annotations for each amplicon region
      for(i in 1:length(amps)){
        amp_i = amps[i]
        # number of segments within the region
        amp_i_seg = subsetByOverlaps(cnv_segments, amp_i)
        amps$nSegments[i] = length(amp_i_seg)
        # median copy number of segments within the region
        q = round(quantile(amp_i_seg$cn, c(0.5,1)))
        amps$medianCN[i] = q["50%"]
        # maximum copy number within the region
        amps$maxCN[i] = q["100%"]
      }
      # add annotations about the whole amplicon (values will be the same for each region within an amplicon)
      amps$size_amplicon = 0
      amps$nChrs_amplicon = 0
      amps$nRegions_amplicon = 0
      amps$medianCN_amplicon = 0 
      amps$cnSpan_amplicon = 0
      amps$cnStates_amplicon = 0
      amps$nSegments_amplicon = 0
      amps$nSVs_amplicon = 0
      amps$nSVsInternal_amplicon = 0
      # for each amplicon id (multiple regions may belong to one amplicon and have the same id)
      for(id in unique(amps$id)){
        
        i = which(amps$id == id)
        
        amp_i = amps[i]
        amp_i_seg = subsetByOverlaps(cnv_segments, amp_i)
        
        # svs belonging to current amplicon
        idx = ((countOverlaps(bp1, amp_i+cnvTol) > 0) & (countOverlaps(bp2, amp_i+cnvTol) > 0))
        bp1_i = bp1[idx]
        bp2_i = bp2[idx]
        
        # size of the amplicon
        amps$size_amplicon[i] = sum(width(amp_i))
        # number of different chromosomes involved in the amplicon structure
        amps$nChrs_amplicon[i] = length(unique(gsub("p|q","",seqnames(amp_i))))
        # number of regions the amplicon consists of
        amps$nRegions_amplicon[i] = length(amp_i)
        # median copy number of the whole amplicon
        q = round(quantile(amp_i_seg$cn, c(0.05,0.5,0.95,1)))
        amps$medianCN_amplicon[i] = q["50%"]
        # copy number span, i.e. range from 5% to 95% quantile (do not count extreme outliers, hence the quantiles here)
        amps$cnSpan_amplicon[i] = (q["95%"] - q["5%"] + 1)
        # number of different copy number states (same as above, remove extreme outliers below/above 5% and 95% quantile)
        amps$cnStates_amplicon[i] = length(unique(amp_i_seg$cn[amp_i_seg$cn >= q["5%"] & amp_i_seg$cn <= q["95%"]]))
        # number of copy number segments within the whole amplicon
        amps$nSegments_amplicon[i] = length(amp_i_seg)
        # total number of svs within the amplicon
        amps$nSVs_amplicon[i] = sum((countOverlaps(bp1_i, amp_i+cnvTol) > 0) & (countOverlaps(bp2_i, amp_i+cnvTol) > 0))
        # total number of internal svs within the amplicon; internal svs have one or both breakpoints within the amplicon, not at the border of regions
        amps$nSVsInternal_amplicon[i] = sum((countOverlaps(bp1_i, amp_i-cnvTol) > 0) | (countOverlaps(bp2_i, amp_i-cnvTol) > 0))
      }
      
      # (4) Filter seismic amplifications by the number of internal structural variations
      # remove low complexity amplicons below the given threshold "minInternalSVs"
      amps = amps[amps$nSVsInternal_amplicon >= minInternalSVs]
      
      if(length(amps) > 0){
        names(amps) = 1:length(amps)
        amps = chrArm2chr(amps, chrs, chrArms)
        
        # return only svs within seismic amplicons
        bp1 = chrArm2chr(bp1, chrs, chrArms)
        bp2 = chrArm2chr(bp2, chrs, chrArms)
        idx = (countOverlaps(bp1, amps+cnvTol) > 0) & (countOverlaps(bp2, amps+cnvTol) > 0)
        bp1 = bp1[idx]
        bp2 = bp2[idx]
        sv = sv[idx,]
        # classify svs as "internal" or "flanking"
        # "flanking": both breakpoints lie at the border of an amplicon region
        # "internal": at least one breakpoint lies within the amplicon region, not at the border
        idxInternal = ((countOverlaps(bp1, amps-cnvTol) > 0) | (countOverlaps(bp2, amps-cnvTol) > 0))
        idxFlanking = ((countOverlaps(bp1, amps-cnvTol) == 0) & (countOverlaps(bp2, amps-cnvTol) == 0))
        sv$type = ""
        sv$type[idxInternal] = "internal"
        sv$type[idxFlanking] = "flanking"
        # add amplicon id to svs
        o = findOverlaps(bp1, amps+cnvTol)
        sv$amplicon_id = 0
        sv$amplicon_id[queryHits(o)] = amps$id[subjectHits(o)]
        
        message("Done")
        message("Seismic amplification detected: ", length(unique(amps$id)), " amplicon(s)")
        return(list("amplicons"=amps, "svs"=sv))
      }
    }
    
  }

  # no seismic amplification found  
  message("Done")
  message("No seismic amplification detected")
  return(list("amplicons"=GRanges(), "svs"=data.frame()))

}

################################################################################################################################################################
################################################################################################################################################################
# Other smaller functions used by the detecion above

# Function used by detect_seismic_amplification() to filter a given GRanges object by another GRanges object
# It is a modification of GRanges setdiff() for filtering and trimming cnvs, but without reduce() and thus keeping the mcols
#
# Parameters:
# gr1: GRanges object
# gr2: GRanges object
#
# Returns:
# GRanges object gr1 without the regions from gr2, but with its initial annotations in mcols
filter_cnv <- function(gr1, gr2){
  gr3 = split(gr1, rep_len(c(1, 2), length(gr1)))
  gr3 = lapply(gr3, function(x){setdiff(x, gr2)})
  if(length(gr3) == 2){
    gr3 = c(gr3[[1]], gr3[[2]])
  }else{
    gr3 = c(gr3[[1]])
  }
  gr3 = gr3[order(gr3)]
  o = findOverlaps(gr3, gr1)
  mcols(gr3) = mcols(gr1[subjectHits(o)])
  return(gr3)
}

# Function used by detect_seismic_amplification() to add the p-/q-arm info to the chromosome name of a GRanges object (e.g. chr12 to chr12p)
#
# Parameters:
# gr: GRanges object, input
# chrs: character vector with all unique chromosome names (without p/q)
# chrArms: GRanges object with start and end of every chromosome arm
#
# Returns:
# GRanges object gr with modified seqnames and seqlevels
chr2chrArm <- function(gr, chrs, chrArms){
  o = findOverlaps(gr, chrArms)
  seqlevels(gr) = c(chrs, chrArms$name)
  seqnames(gr) = Rle(factor(paste0(seqnames(gr), chrArms[subjectHits(o)]$arm), levels=c(chrs, chrArms$name)))
  gr = dropSeqlevels(gr, chrs)
  return(gr)
}

# Function used by detect_seismic_amplification() to remove the p-/q-arm info to the chromosome name of a GRanges object (e.g. chr12 to chr12p)
#
# Parameters:
# gr: GRanges object, input
# chrs: character vector with all unique chromosome names (without p/q)
# chrArms: GRanges object with start and end of every chromosome arm
#
# Returns:
# GRanges object gr with modified seqnames and seqlevels
chrArm2chr <- function(gr, chrs, chrArms){
  seqlevels(gr) = c(chrs, chrArms$name)
  seqnames(gr) = Rle(factor(gsub("p|q", "", seqnames(gr)), levels=c(chrs, chrArms$name)))
  gr = dropSeqlevels(gr, chrArms$name)
  return(gr)
}

# Function used by detect_seismic_amplification() to split a GRanges object at given positions, it adds a new line with every split
#
# Parameters:
# gr: GRanges object, input
# pos: GRanges object with regions whose starts and ends are used to split gr
#
# Returns:
# GRanges object gr with regions split at given positions
split_regions <- function(gr, pos){
  if(length(pos) > 0){
    pos = c(GRanges(seqnames(pos), IRanges(start(pos),start(pos))), GRanges(seqnames(pos), IRanges(end(pos),end(pos))))
    for(i in 1:length(pos)){
      posi = pos[i]
      o = findOverlaps(gr, posi)
      for(j in queryHits(o)){
        gr = c(gr, gr[j])
        end(gr[j]) = start(posi)-1
        start(gr[length(gr)]) = start(posi)+1
      }
    }
    gr = gr[width(gr) > 0]
  }
  return(sort(gr))
}

# Function used by detect_seismic_amplification() to identify (amplified) regions which are connected by structural variants
#
# Parameters:
# gr: GRanges object with (amplified) regions
# bp1: GRanges object with the first breakpoint of structural variants
# bp2: GRanges object with the second breakpoint of structural variants
# tol: numeric, maximum distance between a breakpoint and the region in gr (to be sensitive for inaccuracies in cnv and breakpoint detections)
#
# Returns:
# numeric vector, a set of ids, one for each region; regions connected by structural variants share the same id
get_connected_regions <- function(gr, bp1, bp2, tol=5000){
  gr$id = 1:length(gr)
  if(length(gr) > 0){
    connectedAmp = list()
    gr = gr+tol
    for(i in 1:length(gr)){
      gri = gr[i]
      connectedAmp[[i]] = gri$id
      idx1 = queryHits(findOverlaps(bp2, gri))
      idx2 = queryHits(findOverlaps(bp1, gri))
      gr_connectedAmp = subsetByOverlaps(gr, c(bp1[idx1], bp2[idx2]))
      connectedAmp[[i]] = c(connectedAmp[[i]], gr_connectedAmp$id)
    }
    m = sapply(connectedAmp,function(x) sapply(connectedAmp,function(y) length(intersect(x,y))>0))
    amps = groups(components(graph_from_adjacency_matrix(m)))
    ampIDs = vector(mode="numeric", length=length(gr))
    for(i in 1:length(amps)){
      ampsi = amps[[i]]
      for(j in ampsi){
        ampIDs[j] = i
      }
    }
  }
  return(ampIDs)
}
