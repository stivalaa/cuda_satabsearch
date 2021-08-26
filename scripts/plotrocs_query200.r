# plotrocs_query200.r - plot ROCs for different methods on query200 data set
#
# Alex Stivala, March 2009
#
# Plot ROC curves for different methods on the query200 query set in
# ASTRAL SCOP 95% sequence nr data set
#
# Requires the ROCR package from CRAN (developed with version 1.0-2)
# (ROCR in turn requires gplots, gtools, gdata)
#
# Run this on the output of e.g. tsevalfn.py with the -l option,
# it is a table with one column of scores from classifier, and second
# column of true class label (0 or 1)
#
# The citation for the ROCR package is
#   Sing et al 2005 "ROCR: visualizing classifier performance in R"
#   Bioinformatics 21(20):3940-3941
# 
# 
# $Id: plotrocs_query200.r 2376 2009-05-14 01:40:32Z astivala $
 

library(ROCR)

#
# globals
#

colorvec=c('deepskyblue4','brown','red','turquoise','blue','purple','green','cyan','gray20','magenta','darkolivegreen2','midnightblue','magenta3','darkseagreen','violetred3','darkslategray3')
ltyvec=c(1,2,4,5,6,1,2,1,5,6,1,2,4,5,6,1,2)
namevec=c('QP tableau search (norm2)','VAST','SHEBA','TableauSearch (norm2)', 'TOPS')
slrtabs=c('query200/norm2/query200_roc.slrtab','../other_results/vast/vast_query200_res/vast_query200_roc.slrtab','../other_results/sheba/query200-pdbstyle-sel-gs-bib-95-1.73/sheba_query200_roc.slrtab','../other_results/TableauSearch/query200/norm2/tabsearch_query200_roc.slrtab','../other_results/tops/query200/tops_query200_roc.slrtab')

#
# functions
#

#
# Return the ROCR performance object for plotting ROC curve
#
# Parameters:
#   tab : data frame with score and label columns
# 
# Return value:
#   ROCR performance object with FPR and TPR for plotting ROC curve
#
compute_perf <- function(tab)
{
    # tab is a data frame with score and label columns
    pred <- prediction(tab$score, tab$label)
    perfroc <- performance(pred, measure="tpr",x.measure="fpr")
    return(perfroc)
}

#
# main
#



# EPS suitable for inserting into LaTeX

postscript('rocs_query200.eps',
           onefile=FALSE,paper="special",horizontal=FALSE, 
           width = 9, height = 6)
for (i in 1:length(slrtabs)) {
    tab <- read.table(slrtabs[i], header=TRUE)
    perfroc <- compute_perf(tab)
    plot(perfroc, lty=ltyvec[i], col=colorvec[i], add=(i>1),
         downsampling=0.5,
#         main='ROC for 200 queries in ASTRAL SCOP 95% sequence identity nonredundant data set' # no title since including in paper with caption
         )
}
legend('bottomright', col=colorvec, lty=ltyvec, legend=namevec)
#lines(c(0,1),c(0,1),type='l',lty=3)
dev.off()

