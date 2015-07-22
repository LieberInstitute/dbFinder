library('GenomicRanges')
library('devtools')

# Read data
info <- read.csv('journal.pgen.1003433.s006.csv')

# Construct GRanges object and save it
chr <- gsub(':.*', '', info$Peak.Coordinates)
pos <- gsub('chr[0-9|X|Y][0-9]*:', '', info$Peak.Coordinates)
start <- as.numeric(gsub('-.*', '', pos))
end <- as.numeric(gsub('.*-', '', pos))
neunPeaks <- GRanges(seqnames = chr, IRanges(start = start, end = end), strand = '*')
mcols(neunPeaks) <- info[, -1]
neunPeaks
save(neunPeaks, file = 'neunPeaks.Rdata')

# Reproducibility info
options(width = 120)
session_info()
Sys.time()

#> library('GenomicRanges')
#Loading required package: BiocGenerics
#Loading required package: parallel
#
#Attaching package: ‘BiocGenerics’
#
#The following objects are masked from ‘package:parallel’:
#
#    clusterApply, clusterApplyLB, clusterCall, clusterEvalQ, clusterExport, clusterMap, parApply, parCapply, parLapply, parLapplyLB,
#    parRapply, parSapply, parSapplyLB
#
#The following object is masked from ‘package:stats’:
#
#    xtabs
#
#The following objects are masked from ‘package:base’:
#
#    anyDuplicated, append, as.data.frame, as.vector, cbind, colnames, do.call, duplicated, eval, evalq, Filter, Find, get, intersect,
#    is.unsorted, lapply, Map, mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank, rbind, Reduce,
#    rep.int, rownames, sapply, setdiff, sort, table, tapply, union, unique, unlist, unsplit
#
#Loading required package: S4Vectors
#Loading required package: stats4
#Loading required package: IRanges
#Loading required package: GenomeInfoDb
#> library('devtools')
#> 
#> # Read data
#> info <- read.csv('journal.pgen.1003433.s006.csv')
#> 
#> # Construct GRanges object and save it
#> chr <- gsub(':.*', '', info$Peak.Coordinates)
#> pos <- gsub('chr[0-9|X|Y][0-9]*:', '', info$Peak.Coordinates)
#> start <- as.numeric(gsub('-.*', '', pos))
#> end <- as.numeric(gsub('.*-', '', pos))
#> neunPeaks <- GRanges(seqnames = chr, IRanges(start = start, end = end), strand = '*')
#> mcols(neunPeaks) <- info[, -1]
#> neunPeaks
#GRanges object with 1157 ranges and 11 metadata columns:
#         seqnames                 ranges strand   |    Length TSSGenomicPosition NearestTSS DistanceToTSS Proximal.Distal AveragePrenatal
#            <Rle>              <IRanges>  <Rle>   | <integer>          <integer>   <factor>     <integer>        <factor>       <numeric>
#     [1]     chr1 [110198190, 110199797]      *   |      1607          110198697      GSTM4             0        Proximal     0.017112328
#     [2]     chr1 [110452725, 110454687]      *   |      1962          110453232       CSF1             0        Proximal     0.005020967
#     [3]     chr1 [154377344, 154379306]      *   |      1962          154377668       IL6R             0        Proximal     0.006345178
#     [4]     chr1 [156610834, 156614446]      *   |      3612          156611739       BCAN             0        Proximal     0.006783004
#     [5]     chr1 [ 15735487,  15740224]      *   |      4737           15736390      EFHD2             0        Proximal     0.010312704
#     ...      ...                    ...    ... ...       ...                ...        ...           ...             ...             ...
#  [1153]    chr19 [ 23253131,  23254620]      *   |      1489           22952784      ZNF99        300347          Distal      0.01199730
#  [1154]     chr2 [134022446, 134025233]      *   |      2787          134326031     NCKAP5        300798          Distal      0.03507908
#  [1155]    chr19 [ 23257518,  23258919]      *   |      1401           22952784      ZNF99        304734          Distal      0.05137412
#  [1156]     chr4 [ 82964128,  82966354]      *   |      2226           83295149     HNRNPD        328795          Distal      0.03473416
#  [1157]     chr8 [ 49292251,  49293758]      *   |      1507           49647870     EFCAB1        354112          Distal      0.01418229
#         AverageInfant  AverageOld  LogRatio        Ttest Direction
#             <numeric>   <numeric> <numeric>    <numeric>  <factor>
#     [1]    0.02486927  0.03911535 -1.147182 4.870000e-06        Up
#     [2]    0.01094479  0.01177156 -1.084868 4.845153e-03        Up
#     [3]    0.00944907  0.01961488 -1.488816 4.050000e-14        Up
#     [4]    0.01055625  0.02001794 -1.433223 2.300000e-10        Up
#     [5]    0.01862396  0.02562585 -1.234884 4.300000e-17        Up
#     ...           ...         ...       ...          ...       ...
#  [1153]   0.007633866 0.002199031  2.022505  0.021097842      Down
#  [1154]   0.032532566 0.013445828  1.320510  0.004677639      Down
#  [1155]   0.037834675 0.012145966  1.994234  0.000148264      Down
#  [1156]   0.028601378 0.013827794  1.269000  0.000407679      Down
#  [1157]   0.012563380 0.005315922  1.265324  0.022638812      Down
#  -------
#  seqinfo: 23 sequences from an unspecified genome; no seqlengths
#> save(neunPeaks, file = 'neunPeaks.Rdata')
#> 
#> # Reproducibility info
#> options(width = 120)
#> session_info()
#Session info-----------------------------------------------------------------------------------------------------------
# setting  value                                             
# version  R Under development (unstable) (2014-11-01 r66923)
# system   x86_64, darwin10.8.0                              
# ui       AQUA                                              
# language (EN)                                              
# collate  en_US.UTF-8                                       
# tz       America/Detroit                                   
#
#Packages---------------------------------------------------------------------------------------------------------------
# package       * version date       source        
# BiocGenerics  * 0.13.11 2015-04-03 Bioconductor  
# devtools      * 1.6.1   2014-10-07 CRAN (R 3.2.0)
# GenomeInfoDb  * 1.3.16  2015-03-27 Bioconductor  
# GenomicRanges * 1.19.52 2015-04-04 Bioconductor  
# IRanges       * 2.1.43  2015-03-07 Bioconductor  
# rstudioapi      0.3.1   2015-04-07 CRAN (R 3.2.0)
# S4Vectors     * 0.5.22  2015-03-06 Bioconductor  
# XVector         0.7.4   2015-02-08 Bioconductor  
#> Sys.time()
#[1] "2015-05-12 10:49:06 EDT"
#> 
