# rocauc.r - plot ROC and compute AUC given score+label table with R using RROC
#
# Alex Stivala, October 2008
#
# Requires the ROCR package from CRAN (developed with version 1.0-2)
# (ROCR in turn requires gplots, gtools, gdata)
#
# Run this on the output of e.g. tsevalfn.py with the -l option,
# it is a table with one column of scores from classifier, and second
# column of true class label (0 or 1)
#
# Uses commandArgs() R function to get trailing arguments from R command
# line ie after the --args option. The filename of the .slrtab file
# is obtained from --args, and the output file is constructed form it
# eg foo.slrtab results in foo.eps file for ROC plot with
# 
#    R --vanilla -f rocauc.r --args foo.slrtab
#
#
# The citation for the ROCR package is
#   Sing et al 2005 "ROCR: visualizing classifier performance in R"
#   Bioinformatics 21(20):3940-3941
# 
# 
# Also calculate stdand eror (and 95% confidence inteval) for the AUC,
# according to Hanley-McNeil method:
#   Hanley & McNeil 1982 "The Meaning and Use of the Area under a Receiver
#   Operating Characteristic (ROC) Curve" Radiology 143(1):29-36
# ROCR has all sorts of featuers but this isn't one of them so had
# to implement it here.
#

# $Id: rocauc.r 3606 2010-05-04 06:03:55Z alexs $

library(ROCR)

#
# functions
#

#
# compute AUC and standard error and 95% CI by Hanley-McNeil method,
# using R wilcox.test for the Wilcoxon rank-sum (aka Mann-Whitney) test
# 
# Parameters:
#   tab : data frame with score and label columns
# 
# Return value:
#   list with auc and stderror members
#
compute_auc_error <- function(tab)
{
    # tab is a data frame with score and label columns
    x <- tab$score[tab$label == 1]
    y <- tab$score[tab$label == 0]
    wilcox <- wilcox.test(x,y)
    nA <- length(x)
    nN <- length(y)
    stopifnot(nA + nN == length(tab$label))
    nA <- as.double(nA)
    nN <- as.double(nN)
    U <- wilcox$statistic  # Mann-Whitney U statistic
    theta <- U / (nA * nN) # AUC
    theta2 <- theta*theta
    Q1 <- theta / (2 - theta)
    Q2 <- 2*theta2 / (1 + theta)
    SE2 <- (theta*(1-theta) + (nA - 1)*(Q1 - theta2) + (nN - 1)*(Q2 - theta2)) /
          (nA*nN)
    SE <- sqrt(SE2)

    retval <- list()
    retval$auc <- theta
    retval$stderror <- SE
    return(retval)
}

#
# plot the ROC curve and compute AUC using ROCR
#
plotroc <- function(filename)
{
    tab <- read.table(filename, header=TRUE)
    # tab is a data frame with score and label columns
    pred <- prediction(tab$score, tab$label)
    perfroc <- performance(pred, measure="tpr",x.measure="fpr")
    perfauc <- performance(pred, measure="auc")

    # EPS suitable for inserting into LaTeX
    postscript(sub('[.]slrtab$','.eps',filename),
               onefile=FALSE,paper="special",horizontal=FALSE, 
               width = 9, height = 6)
    plot(perfroc)
    auc <- perfauc@y.values
    legend('bottomright', legend=paste('AUC =', format(auc, digits=4)), bty='n')
    dev.off()

    aucerr <- compute_auc_error(tab)
    quantile <- 1.96 # 0.975 quantile of normal distrib
    low95 <- aucerr$auc - quantile * aucerr$stderr
    high95 <- aucerr$auc + quantile * aucerr$stderr

    cat(filename,':\n')
    cat('RROC          AUC =',format(perfauc@y.values,digits=4),'\n')
    cat('Hanley-McNeil AUC =',format(aucerr$auc,digits=4),'\n')
    cat('       std. error =',format(aucerr$stderror,digits=5),'\n')
    cat('           95% CI =',format(low95,digits=4),',',
                               format(high95,digits=4),'\n')
}


#
# main
#
filename <- commandArgs(trailingOnly=TRUE)
plotroc(filename)

