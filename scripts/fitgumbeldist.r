#############################################################################3
#
# File:    fitgumbeldist.r
# Author:  Alex Stivala
# Created: July 2010
#
# fitgumbeldist.r - compute Gumbest dist. parameters and build plots
#
# R script for fitting Gumbel distribution (finding parameters location
# and scale) and plotting histograms and fit over them.
# Also, compute the z-score and p-values using this, then
# compute coverage-vs-EPQ graph  as used in Ortiz et al 2002 (and earlier
# Levitt and Gerstein (1998) to compare the analytic and observed p-values)
#
#
# Uses commandArgs() R function to get trailing arguments from R command
# line ie after the --args option. The filename of the .slrtab file
# (score and label)
# is obtained from --args.
#
# label is used to separate same fold (1) and different fold (0) scores,
# different distributinos fitted to each
#
# outputs stuff on stdout (including the parameter values) and
# creates PostScript file with histograms and curves, with filename
# derived form input filename: basename_gumbel.eps
# 
#    R --vanilla -f plotssehistogram.r --args query200.slrtab
#
# creates query200_gumbel.eps
# and     query200_epq_gumbel.eps
#
# Uses the evir R package for gumbel fitting and the evd package for dbumbel()
#
# $Id: fitgumbeldist.r 3955 2010-07-21 05:52:47Z alexs $
#############################################################################3

library(evir)
library(evd)


#############################################################################3
#
# constants
#
#############################################################################3

colorvec=c('deepskyblue4','brown','red','turquoise','blue','purple','green','cyan','gray20','magenta','darkolivegreen2','midnightblue','magenta3','darkseagreen','violetred3','darkslategray3')
ltyvec=c(1,2,4,5,6,1,2,1,5,6,1,2,4,5,6,1,2)
pchvec=c(20,21,22,23,24,25,19,21,22)

eulergamma <- 0.5772156649015328606

#############################################################################3
#
# functions
#
#############################################################################3


# 
# z_gumbel() -  compute Z-score from Gumbel distribution
#
# Parameters:
#    x - score to compute Z-score for
#    gumbel - evd object for Gumbel distribution from gumbel()
#
# Return value:
#    Z-score computed for x according to gum distribution
#
z_gumbel <- function(x, gum )
{
      a<- gum0$par.ests[2]
      b <- gum0$par.ests[1]
      mu <- a + b * eulergamma
      sigma <- ( (pi  / sqrt(6)) * b )
      z <- (x - mu)/sigma
      return (z)
}

#
# pv_gumbel() - compute P-value for Z-score from Gumbel distribution
#
# Parameters: 
#    z - Z-score from z_gumbel()
#
# Return value:
#    P-value for the Z-score
#
pv_gumbel <- function(z)
{
    return ( 1 - exp(-exp(-( (pi/sqrt(6)*z + eulergamma) ) ) ) )
}

#############################################################################3
#
# main
#
#############################################################################3


filename <- commandArgs(trailingOnly=TRUE)

# EPS suitable for inserting into LaTeX
postscript(sub('[.]slrtab$','_gumbel.eps',filename),
           onefile=FALSE,paper="special",horizontal=FALSE, 
           width = 9, height = 6)

slrtab <- read.table(filename, header=T)
slrtab <- subset(slrtab, score >= 0) # one or two -ve scores sometimes

#
# Fit Gumbel distribution to same fold scores and different fold scores
#

gum0 <- gumbel(subset(slrtab, label == 0)$score)
gum1 <- gumbel(subset(slrtab, label == 1)$score)

cat('different folds: a = ', format(gum0$par.ests[2], digits=16), ' ' )
cat('b = ', format(gum0$par.ests[1], digits=16), '\n')
cat('same fold: a = ', format(gum1$par.ests[2], digits=16), ' ' )
cat('b = ', format(gum1$par.ests[1], digits=16), '\n')

#
# Plot histograms and Gumbel distributions fitted to them
#

h0 <- hist(subset(slrtab, label == 0)$score, breaks=20, plot=F)
h1 <- hist(subset(slrtab, label == 1)$score, breaks=20, plot=F)
xh0 <- c(min(h0$breaks), h0$breaks)
yh0 <- c(0, h0$density, 0)
xh1 <- c(min(h1$breaks), h1$breaks)
yh1 <- c(0, h1$density, 0)

plot(xh0, yh0, type='s', lty=ltyvec[1], col=colorvec[1], xlab='score',ylab='frequency', ylim=c(0,1.0)) # ylim=c(0, max(yh0, yh1)) )
lines(xh1, yh1, type='s', lty=ltyvec[2], col=colorvec[2])

curve(dgumbel(x, gum0$par.ests[2], gum0$par.ests[1]), col=colorvec[3], lty=ltyvec[3], add=T)
curve(dgumbel(x, gum1$par.ests[2], gum1$par.ests[1]), col=colorvec[4], lty=ltyvec[4], add=T)

legend("topright", lty=ltyvec, col=colorvec, bty='n',
       legend=c("Different folds, histogram",
         "Same fold, histogram",
         "Different folds, fitted Gumbel distribution",
         "Same fold, fitted Gumbel distribution") )
         

dev.off()

#
# compute P-values for all scores
#

pvslrtab <- slrtab
pvslrtab$pvalue <- pv_gumbel(z_gumbel(slrtab$score, gum0))

#
# compute coverage-vs-EPQ graph  as used in Ortiz et al 2002 (and earlier
# Levitt and Gerstein (1998) to compare the analytic and observed p-values)
#

# sort by ascending P-value so 'best' are at start of list
pvslrtab <- pvslrtab[sort.list(pvslrtab$pvalue),]

tptotal <- length(subset(pvslrtab, label == 1)$label)
tpcount <- 0
fpcount <- 0

coverage_vec <- c()
obs_pv_vec <- c()
analytic_pv_vec <- c()

for (i in 1:length(pvslrtab$label))
{
 if (pvslrtab$label[i] == 0)
 {
 fpcount <- fpcount + 1
 }
 else
 {
 tpcount <- tpcount + 1
 }
 observed_pvalue <- fpcount / (length(pvslrtab$label) - tptotal)
 coverage <- tpcount / tptotal

 coverage_vec <- c(coverage_vec, coverage)
 obs_pv_vec <- c(obs_pv_vec, observed_pvalue)
 analytic_pv_vec <- c(analytic_pv_vec, pvslrtab$pvalue[i])

 #cat (sprintf("%d\t%f\t%f\t%f\n", i, pvslrtab$pvalue[i], observed_pvalue, coverage) )

 if (observed_pvalue > 0.05)
   {
    break
  }
}

# EPS suitable for inserting into LaTeX
postscript(sub('[.]slrtab$','_epq_gumbel.eps',filename),
           onefile=FALSE,paper="special",horizontal=FALSE, 
           width = 9, height = 6)

plot(analytic_pv_vec*100, coverage_vec*100, type='n', xlab='Errors per query (%)',
     ylab='Coverage (%)')

lines(analytic_pv_vec*100, coverage_vec*100, col=colorvec[1], lty=ltyvec[1])
lines(obs_pv_vec*100, coverage_vec*100, col=colorvec[2], lty=ltyvec[2])

legend('topleft', col=colorvec, lty=ltyvec, legend=c('Analytic','Observed'),
       bty='n')

dev.off()


