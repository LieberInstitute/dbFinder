# Rscript ranges.R > logs/ranges_log.txt 2>&1

library('IRanges')
library('RColorBrewer')

plotRanges <- function(x, xlim = x, main = deparse(substitute(x)),
                       col = "black", sep = 0.5, ...) {
    height <- 1
    if (is(xlim, "Ranges"))
    xlim <- c(min(start(xlim)), max(end(xlim)))
    plot.new()
    bins <- disjointBins(IRanges(start(x), end(x) + 1))
    bins.tmp <- bins
    bins.tmp[bins == 1] <- 2
    bins.tmp[bins == 2] <- 1
    bins <- bins.tmp
    i <- seq_len(length(x))
    
    plot.window(xlim, c(0, max(bins)*(height + sep)))
    ybottom <- bins * (sep + height) - height
    rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col = col[i], border = col[i], ...)
    title(main)
    axis(1)
}

start <- 1e4 + c(500, 530, 540, 550, 730, 750, 770, 790, 1000, 1010, 1030, 1050)
end <- start + 300
ir <- IRanges(start, end)
ir2 <- IRanges(min(start), max(end))
ir3 <- c(ir, ir2)

bins <- disjointBins(IRanges(c(start, min(start)), c(end, max(end)) + 1))

pdf(file.path('plots', 'ranges.pdf'))
plotRanges(ir3, main = '', xlim = c(start(range(ir)), end(range(ir))) + c(-1, 1) * 150, col = rep(c(brewer.pal(4, 'Set1')[2:4], 'black'), c(4, 4, 4, 1)))
dev.off()

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()