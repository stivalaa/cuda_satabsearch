# R script for plotting ROC curve for different tableau search methods
# Alex Stivala
# July 2008
# Reduced version to make clearer in grayscale (showing only those that
# don't overlap too much).
# $Id: plotsearchroc_betagrasp_reduced.r 1946 2008-10-04 05:14:43Z astivala $


#
# globals
#

colorvec=c('black','black','black','black','black','black')
ltyvec=c(1,2,4,5,6)
namevec=c('QP numeric ordering SSE type',  'QP discrete SSE type', 'QP discrete ordering SSE type')


#
# main
#

qp_numeric <- read.table('betagrasp.tsrchn.ff.rtab',header=TRUE)
qp_numeric_ordering_penalizessetype <- read.table('betagrasp.tsrchn.tt.rtab',header=TRUE)

qp_discrete <- read.table('betagrasp.tsrchd.ff.rtab',header=TRUE)
qp_discrete_penalizessetype <- read.table('betagrasp.tsrchd.tf.rtab', header=TRUE)
qp_discrete_ordering_penalizessetype <- read.table('betagrasp.tsrchd.tt.rtab', header=TRUE)

# EPS suitable for inserting into LaTeX
postscript('searchroc_betagrasp_reduced.eps',onefile=FALSE,paper="special",horizontal=FALSE, 
           width = 9, height = 6)

plot(c(0,1),c(0,1),type='l',xlim=c(0,1),ylim=c(0,1),xlab="False Positive Rate",ylab="True Positive Rate",lty=3,main='ROC curves for query d1ubia_ against ASTRAL 95% seq id ProSMoS dataset as truth')

lines(qp_numeric_ordering_penalizessetype$fpr, qp_numeric_ordering_penalizessetype$tpr, col=colorvec[1], lty=ltyvec[1])
lines(qp_discrete_penalizessetype$fpr, qp_discrete_penalizessetype$tpr, col=colorvec[2], lty=ltyvec[2])
lines(qp_discrete_ordering_penalizessetype$fpr, qp_discrete_ordering_penalizessetype$tpr, col=colorvec[3], lty=ltyvec[3])

legend('bottomright', col=colorvec, lty=ltyvec, legend=namevec)


dev.off()

