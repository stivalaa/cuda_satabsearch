# R script for plotting SSE number histogram
# Alex Stivala
# August 2008
# $Id: plotssehistogram.r 3001 2009-12-01 05:41:12Z alexs $
#
# Uses commandArgs() R function to get trailing arguments from R command
# line ie after the --args option. The filename of the .ssenumslist file
# (which is just a list of numbers of SSEs in each domatin, one per line)
# is obtained from --args, and the output file is constructed form it
# eg astral-sel-gs-95-1.73.ssenumslist results in 
# astral-sel-gs-95-1.73-ssehistogram.eps with
# 
#    R --vanilla -f plotssehistogram.r --args astral-sel-gs-95-1.73.ssenumlist
#

filename <- commandArgs(trailingOnly=TRUE)

# EPS suitable for inserting into LaTeX
postscript(sub('[.]ssenumslist$','-ssehistogram.eps',filename),
           onefile=FALSE,paper="special",horizontal=FALSE, 
           width = 9, height = 6)
ssenums <- read.table(filename)

if (length(grep("1.73",filename)) > 0)  ymax <- 1000 else  ymax <- 1500

hist(ssenums$V1,freq=TRUE,density=20,main=NULL,xlab='Tableau size (number of SSEs in domain)',breaks=100, ylim=c(0,ymax), xlim=c(0,100))
dev.off()

