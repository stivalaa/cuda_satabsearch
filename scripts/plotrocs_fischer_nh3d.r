# plotrocfischer.r - plot ROCs for different method on Fischer and Nh3D data
#
# Alex Stivala, October 2008
#
# Plot ROC curves for different methods on the Fischer data set
# (Fischer et al 1996) and Nh3D data set (Thiruv et al 2005)
#  as used in Pelta et al 2008.
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
# $Id: plotrocs_fischer_nh3d.r 2376 2009-05-14 01:40:32Z astivala $
 

library(ROCR)

#
# globals
#

colorvec=c('deepskyblue4','brown','red','turquoise','blue','purple','green','cyan','gray20','magenta','darkolivegreen2','midnightblue','magenta3','darkseagreen','violetred3','darkslategray3')
ltyvec=c(1,2,4,5,6,1,2,1,5,6,1,2,4,5,6,1,2)
namevec=c('MSVNS3 norm1','MSVNS3 norm2', 'MSVNS3 norm3', 'QP tableau search norm1', 'QP tableau search norm2', 'QP tableau search norm3')

fischer_fold_files=c('../maxcmo_results/fischer/norm1/fold.slrtab', '../maxcmo_results/fischer/norm2/fold.slrtab','../maxcmo_results/fischer/norm3/fold.slrtab','fischer/norm1/fold.slrtab','fischer/norm2/fold.slrtab','fischer/norm3/fold.slrtab')
fischer_class_files=c('../maxcmo_results/fischer/norm1/class.slrtab','../maxcmo_results/fischer/norm2/class.slrtab','../maxcmo_results/fischer/norm3/class.slrtab','fischer/norm1/class.slrtab','fischer/norm2/class.slrtab','fischer/norm3/class.slrtab')
nh3d_arch_files=c('../maxcmo_results/nh3d/norm1/arch.slrtab','../maxcmo_results/nh3d/norm2/arch.slrtab','../maxcmo_results/nh3d/norm3/arch.slrtab','nh3d/norm1/arch.slrtab','nh3d/norm2/arch.slrtab','nh3d/norm3/arch.slrtab')
nh3d_class_files=c('../maxcmo_results/nh3d/norm1/class.slrtab','../maxcmo_results/nh3d/norm2/class.slrtab','../maxcmo_results/nh3d/norm3/class.slrtab','nh3d/norm1/class.slrtab','nh3d/norm2/class.slrtab','nh3d/norm3/class.slrtab')

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

postscript('rocs_fischer_fold.eps',
           onefile=FALSE,paper="special",horizontal=FALSE, 
           width = 9, height = 6)
for (i in 1:length(fischer_fold_files)) {
    tab <- read.table(fischer_fold_files[i], header=TRUE)
    perfroc <- compute_perf(tab)
    plot(perfroc, lty=ltyvec[i], col=colorvec[i], add=(i>1),
         #main='Fischer data set at fold level' # remove title for paper
         )
}
legend('bottomright', col=colorvec, lty=ltyvec, legend=namevec)
#lines(c(0,1),c(0,1),type='l',lty=3)
dev.off()


postscript('rocs_fischer_class.eps',
           onefile=FALSE,paper="special",horizontal=FALSE, 
           width = 9, height = 6)
for (i in 1:length(fischer_class_files)) {
    tab <- read.table(fischer_class_files[i], header=TRUE)
    perfroc <- compute_perf(tab)
    plot(perfroc, lty=ltyvec[i], col=colorvec[i], add=(i>1),
#         main='Fischer data set at class level'  # remove title for paper
         )
}
legend('bottomright', col=colorvec, lty=ltyvec, legend=namevec)
#lines(c(0,1),c(0,1),type='l',lty=3)
dev.off()

postscript('rocs_nh3d_arch.eps',
           onefile=FALSE,paper="special",horizontal=FALSE, 
           width = 9, height = 6)
for (i in 1:length(nh3d_arch_files)) {
    tab <- read.table(nh3d_arch_files[i], header=TRUE)
    perfroc <- compute_perf(tab)
    plot(perfroc, lty=ltyvec[i], col=colorvec[i], add=(i>1),
         # main='Nh3D data set at architecture level' # remove title for paper
         )
}
legend('bottomright', col=colorvec, lty=ltyvec, legend=namevec)
#lines(c(0,1),c(0,1),type='l',lty=3)
dev.off()


postscript('rocs_nh3d_class.eps',
           onefile=FALSE,paper="special",horizontal=FALSE, 
           width = 9, height = 6)
for (i in 1:length(nh3d_class_files)) {
    tab <- read.table(nh3d_class_files[i], header=TRUE)
    perfroc <- compute_perf(tab)
    plot(perfroc, lty=ltyvec[i], col=colorvec[i], add=(i>1),
         # main='Nh3D data set at class level' #  removed title for paper
         )
         
}
#lines(c(0,1),c(0,1),type='l',lty=3)
legend('bottomright', col=colorvec, lty=ltyvec, legend=namevec) 
dev.off()

