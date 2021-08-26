# R script for plotting ROC curve for different tableau search methods
# Alex Stivala
# July 2008
# $Id: plotsearchroc.r 1946 2008-10-04 05:14:43Z astivala $


#
# globals
#

colorvec=c('deepskyblue4','brown','red','turquoise','blue','purple','green','cyan','gray20','magenta','darkolivegreen2','midnightblue','magenta3','darkseagreen','violetred3','darkslategray3')
ltyvec=c(1,2,4,5,6,1,2,1,5,6,1,2,4,5,6,1,2)
namevec=c('QP numeric', 'QP numeric SSE type', 'QP numeric ordering SSE type', 'QP discrete','QP discrete SSE type', 'QP discrete ordering SSE type')


#
# main
#

qp_numeric <- read.table('d1ubia_.tsrchn.ff.rtab',header=TRUE)
qp_numeric_penalizessetype<-read.table('d1ubia_.tsrchn.tf.rtab',header=TRUE)
qp_numeric_ordering_penalizessetype <- read.table('d1ubia_.tsrchn.tt.rtab', header=TRUE)

qp_discrete <- read.table('d1ubia_.tsrchd.ff.rtab',header=TRUE)
qp_discrete_penalizessetype<-read.table('d1ubia_.tsrchd.tf.rtab',header=TRUE)
qp_discrete_ordering_penalizessetype <- read.table('d1ubia_.tsrchd.tt.rtab', header=TRUE)


# EPS suitable for inserting into LaTeX
postscript('searchroc.eps',onefile=FALSE,paper="special",horizontal=FALSE, 
           width = 9, height = 6)

plot(c(0,1),c(0,1),type='l',xlim=c(0,1),ylim=c(0,1),xlab="False Positive Rate",ylab="True Positive Rate",lty=3,main='ROC curves for query d1ubia_ against ASTRAL 95% seq id SCOP as truth')

lines(qp_numeric$fpr, qp_numeric$tpr, col=colorvec[1], lty=ltyvec[1])
lines(qp_numeric_penalizessetype$fpr, qp_numeric_penalizessetype$tpr, col=colorvec[2], lty=ltyvec[2])
lines(qp_numeric_ordering_penalizessetype$fpr, qp_numeric_ordering_penalizessetype$tpr, col=colorvec[3], lty=ltyvec[3])
lines(qp_discrete$fpr, qp_discrete$tpr, col=colorvec[4], lty=ltyvec[4])
lines(qp_discrete_penalizessetype$fpr, qp_discrete_penalizessetype$tpr, col=colorvec[5], lty=ltyvec[5])
lines(qp_discrete_ordering_penalizessetype$fpr, qp_discrete_ordering_penalizessetype$tpr, col=colorvec[6], lty=ltyvec[6])


legend('bottomright', col=colorvec, lty=ltyvec, legend=namevec)


dev.off()

