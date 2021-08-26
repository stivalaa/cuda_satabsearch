# R script for plotting ROC curve for different tableau search methods
# Alex Stivala
# July 2008
# $Id: plotsearchroc_folds_nodistmatrix_discrete.r 2376 2009-05-14 01:40:32Z astivala $


#
# globals
#

colorvec=c('deepskyblue4','brown','red','turquoise','blue','purple','green','cyan','gray20','magenta','darkolivegreen2','midnightblue','magenta3','darkseagreen','violetred3','darkslategray3')
ltyvec=c(1,2,4,5,6,1,2,1,5,6,1,2,4,5,6,1,2)
namevec=c('beta-grasp', 'Immunoglobulin', 'TIM-barrel', 'Plait (ferredoxin)', 'GFP-like', 'Key-barrel', 'Jelly-roll', 'NAD-binding fold')


#
# main
#

betagrasp <- read.table('d1ubia_.nodistmatrix.tsrchd.tt.rtab', header=TRUE)
immunoglobulin <- read.table('d1ae6h1.nodistmatrix.tsrchd.tt.rtab', header=TRUE)
timbarrel <- read.table('d1tima_.nodistmatrix.tsrchd.tt.rtab', header=TRUE)
plait <- read.table('d1bhne_.nodistmatrix.tsrchd.tt.rtab', header=TRUE)
greekkey <- read.table('d1h6rb_.nodistmatrix.tsrchd.tt.rtab', header=TRUE)
keybarrel <- read.table('d1tttb1.nodistmatrix.tsrchd.tt.rtab', header=TRUE)
jellyroll <- read.table('d2phlb1.nodistmatrix.tsrchd.tt.rtab', header=TRUE)
nadbinding <- read.table('d1f6dc_.nodistmatrix.tsrchd.tt.rtab', header=TRUE)


# EPS suitable for inserting into LaTeX
postscript('searchroc_folds_nodistmatrix_discrete.eps',onefile=FALSE,paper="special",horizontal=FALSE, 
           width = 9, height = 6)

#plot(c(0,1),c(0,1),type='l',xlim=c(0,1),ylim=c(0,1),xlab="False Positive Rate",ylab="True Positive Rate",lty=3,main='ROC curves (no distance information) against ASTRAL 95% seq id SCOP as truth')
# no figure title for BMC Bioinformatics
#plot(c(0,1),c(0,1),type='l',xlim=c(0,1),ylim=c(0,1),xlab="False Positive Rate",ylab="True Positive Rate",lty=3)
plot(c(0,1),c(0,1),type='n',xlim=c(0,1),ylim=c(0,1),xlab="False Positive Rate",ylab="True Positive Rate",lty=3)

lines(betagrasp$fpr, betagrasp$tpr, col=colorvec[1], lty=ltyvec[1])
lines(immunoglobulin$fpr, immunoglobulin$tpr, col=colorvec[2], lty=ltyvec[2])
lines(timbarrel$fpr, timbarrel$tpr, col=colorvec[3], lty=ltyvec[3])
lines(plait$fpr, plait$tpr, col=colorvec[4], lty=ltyvec[4])
lines(greekkey$fpr, greekkey$tpr, col=colorvec[5], lty=ltyvec[5])
lines(keybarrel$fpr, keybarrel$tpr, col=colorvec[6], lty=ltyvec[6])
lines(jellyroll$fpr, jellyroll$tpr, col=colorvec[7], lty=ltyvec[7])
lines(nadbinding$fpr, nadbinding$tpr, col=colorvec[8], lty=ltyvec[8])

legend('bottomright', col=colorvec, lty=ltyvec, legend=namevec)


dev.off()

