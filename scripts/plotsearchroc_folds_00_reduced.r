# R script for plotting ROC curve for different tableau search methods
# Alex Stivala
# July 2008
# Reduced version to make clearer in grayscale (showing only those that
# don't overlap too much).
# $Id: plotsearchroc_folds_00_reduced.r 1946 2008-10-04 05:14:43Z astivala $


#
# globals
#

colorvec=c('black','black','black','black','black','black')
ltyvec=c(1,2,5,4,5,6)
namevec=c('beta-grasp', 'Immunoglobulin', 'TIM-barrel', 'Plait (ferredoxin)', 'Key-barrel')


#
# main
#

betagrasp <- read.table('d1ubia_.tsrchn.tt.rtab', header=TRUE)
immunoglobulin <- read.table('d1ae6h1.tsrchn.orderpen0.typepen0.rtab', header=TRUE)
timbarrel <- read.table('d1tima_.tsrchn.orderpen0.typepen0.rtab', header=TRUE)
plait <- read.table('d1bhne_.tsrchn.orderpen0.typepen0.rtab', header=TRUE)
keybarrel <- read.table('d1tttb1.tsrchn.orderpen0.typepen0.rtab', header=TRUE)


# EPS suitable for inserting into LaTeX
postscript('searchroc_folds_reduced.eps',onefile=FALSE,paper="special",horizontal=FALSE, 
           width = 9, height = 6)

plot(c(0,1),c(0,1),type='l',xlim=c(0,1),ylim=c(0,1),xlab="False Positive Rate",ylab="True Positive Rate",lty=3,main='ROC curves for numeric version against ASTRAL 95% seq id SCOP as truth')

lines(betagrasp$fpr, betagrasp$tpr, col=colorvec[1], lty=ltyvec[1])
lines(immunoglobulin$fpr, immunoglobulin$tpr, col=colorvec[2], lty=ltyvec[2])
lines(timbarrel$fpr, timbarrel$tpr, col=colorvec[3], lty=ltyvec[3])
lines(plait$fpr, plait$tpr, col=colorvec[4], lty=ltyvec[4])
lines(keybarrel$fpr, keybarrel$tpr, col=colorvec[5], lty=ltyvec[5])

legend('bottomright', col=colorvec, lty=ltyvec, legend=namevec)


dev.off()

