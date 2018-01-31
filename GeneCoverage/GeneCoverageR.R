#!/usr/bin/Rscript
##########
library(fields)

#Sample selection #####
args <- commandArgs(trailingOnly=TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned
Sample <- args[1]
Sample
beds_folder <- args[2]
beds_folder
plotType <- args[3]
plotType
outDir <- args[4]
outDir
sPath <- paste0(beds_folder, "/", Sample, "/")
outFolder <- paste0(outDir, "/", Sample)
dir.create(outFolder, showWarnings = F, recursive = T)
#####

#read in plus file
plus.input=read.delim(paste0(sPath, Sample, ".plus.dist.1k.bed"),head=F)
#develop strand directional positioning
real.dist=matrix(ifelse(plus.input[,10]=='+',-1*plus.input[,17],plus.input[,17]),ncol=1)
plus.input=cbind(plus.input,real.dist)

#read in minus file
minus.input=read.delim(paste0(sPath, Sample, ".minus.dist.1k.bed"),head=F)
#develop strand directional positioning
real.dist=matrix(ifelse(minus.input[,10]=='+',-1*minus.input[,17],minus.input[,17]),ncol=1)
minus.input=cbind(minus.input,real.dist)

if(plotType=="TTS"){
  #create relative distance measure from coverage data start and entire gene (not just exon or intron) start stop.
plus.input <- subset(plus.input, plus.input$real.dist >= 0)
rel.dist=matrix(ifelse(plus.input[,10]=="+", (plus.input[,2] - (plus.input[,15])), (plus.input[,14] - (plus.input[,2]))), ncol=1)
plus.input=cbind(plus.input,rel.dist)
plus.1k.input <- subset(plus.input, plus.input$rel.dist >= -999)
  
#subset to exons only
plus.exon=subset(plus.1k.input,plus.1k.input$V12=='exon')

#look to subset coverage based on primary strand and input coverage information. Primary = matched, offstrand should be non genic strand
plus.primary=subset(plus.exon,plus.exon$V10=='+')
plus.offstrand=subset(plus.exon,plus.exon$V10=='-')

#stats bin in 300 bins for both promary and offstrand
plus.primary.bin=stats.bin(plus.primary$rel.dist,log(abs(plus.primary[,4])+1),N=200)
ppb=cbind(matrix(plus.primary.bin$centers,ncol=1),plus.primary.bin$stats["mean",])

plus.offstrand.bin=stats.bin(plus.offstrand$rel.dist,log(abs(plus.offstrand[,4])+1),N=200)
pob=cbind(matrix(plus.offstrand.bin$centers,ncol=1),plus.offstrand.bin$stats["mean",])

#create relative distance measure from coverage data start and entire gene (not just exon or intron) start stop.
minus.input <- subset(minus.input, minus.input$real.dist >= 0)
rel.dist=matrix(ifelse(plus.input[,10]=="+", (minus.input[,2] - (minus.input[,15])), (minus.input[,14] - (minus.input[,2]))), ncol=1)
minus.input=cbind(minus.input,rel.dist)
minus.1k.input <- subset(minus.input, minus.input$rel.dist >= -999)
  
#subset to exons only
minus.exon=subset(minus.1k.input,minus.1k.input$V12=='exon')

#look to subset coverage based on primary strand and input coverage information. Primary = matched, offstrand should be non genic strand
minus.primary=subset(minus.exon,minus.exon$V10=='-')
minus.offstrand=subset(minus.exon,minus.exon$V10=='+')

#stats bin in 300 bins for both promary and offstrand using log transformed values
minus.primary.bin=stats.bin(minus.primary$rel.dist,log(abs(minus.primary[,4])+1),N=200)
mpb=cbind(matrix(minus.primary.bin$centers,ncol=1),minus.primary.bin$stats["mean",])

minus.offstrand.bin=stats.bin(minus.offstrand$rel.dist,log(abs(minus.offstrand[,4])+1),N=200)
mob=cbind(matrix(minus.offstrand.bin$centers,ncol=1),minus.offstrand.bin$stats["mean",])

pob[,2]=-pob[,2]
mob[,2]=-mob[,2]

sense=ppb
sense[,2]=sense[,2]+mpb[,2]

antisense=pob
antisense[,2]=antisense[,2]+mob[,2]

out.table <- data.frame(cbind(sense,antisense[,2]))
colnames(out.table) <- c("Position", "Sense", "Antisense")
out.table$SampleName <- Sample

write.csv(out.table, paste0(outFolder, "/",Sample, '_average_coverage.csv'))

pdf(paste0(outFolder, "/",Sample, '_gene_coverge_plot_WC.pdf'),h=10,w=12)
#pdf("test_strands.pdf",h=10,w=12)
plot(x=NULL,y=NULL,xlim=c(-1000,1000),ylim=c(-5,5), main=Sample)
lines(ppb,col=1,lwd=2)
lines(pob,col=2,lwd=2)
lines(mpb,col=3,lwd=2)
lines(mob,col=4,lwd=2)

abline(v=0,lty=2)
abline(h=0,lty=1,col='grey')
legend('topright',c('plus_primary','plus_offstrand','minus_primary','minus_offstrand'),lty=1,col=c(1,2,3,4))
dev.off()

pdf(paste0(outFolder, "/",Sample, '_gene_coverge_plot.pdf'),h=10,w=12)
#pdf("test_sense.pdf",h=10,w=12)
plot(x=NULL,y=NULL,xlim=c(-1000,1000),ylim=c(-5,5),main=Sample)
lines(sense,col=1,lwd=2)
lines(antisense,col=2,lwd=2)

abline(v=0,lty=2)
abline(h=0,lty=1,col='grey')
legend('topright',c('sense','antisense'),lty=1,col=c(1,2))
dev.off()
  
  
  } else {
   #create relative distance measure from coverage data start and entire gene (not just exon or intron) start stop.
rel.dist=matrix(ifelse(plus.input$real.dist==0,
ifelse(plus.input[,10]=="-",
((plus.input[,15] - (plus.input[,2]))/(plus.input[,15] - plus.input[,14]))*1000,
(((plus.input[,2]) - plus.input[,14])/(plus.input[,15] - plus.input[,14]))*1000),
ifelse(plus.input$real.dist>0,
plus.input$real.dist + 1000,
plus.input$real.dist)),
ncol=1)
plus.input=cbind(plus.input,rel.dist)
  
#subset to exons only
plus.exon=subset(plus.input,plus.input$V12=='exon')  

#look to subset coverage based on primary strand and input coverage information. Primary = matched, offstrand should be non genic strand
plus.primary=subset(plus.exon,plus.exon$V10=='+')
plus.offstrand=subset(plus.exon,plus.exon$V10=='-')

#stats bin in 300 bins for both promary and offstrand
plus.primary.bin=stats.bin(plus.primary$rel.dist,log(abs(plus.primary[,4])+1),N=300)
ppb=cbind(matrix(plus.primary.bin$centers,ncol=1),plus.primary.bin$stats["mean",])

plus.offstrand.bin=stats.bin(plus.offstrand$rel.dist,log(abs(plus.offstrand[,4])+1),N=300)
pob=cbind(matrix(plus.offstrand.bin$centers,ncol=1),plus.offstrand.bin$stats["mean",])

#create relative distance measure from coverage data start and entire gene (not just exon or intron) start stop.
rel.dist=matrix(ifelse(minus.input$real.dist==0,
ifelse(minus.input[,10]=="-",
((minus.input[,15] - (minus.input[,2]))/(minus.input[,15] - minus.input[,14]))*1000,
(((minus.input[,2]) - minus.input[,14])/(minus.input[,15] - minus.input[,14]))*1000),
ifelse(minus.input$real.dist>0,
minus.input$real.dist + 1000,
minus.input$real.dist)),
ncol=1)
minus.input=cbind(minus.input,rel.dist)
  
#subset to exons only
minus.exon=subset(minus.input,minus.input$V12=='exon')  

#look to subset coverage based on primary strand and input coverage information. Primary = matched, offstrand should be non genic strand
minus.primary=subset(minus.exon,minus.exon$V10=='-')
minus.offstrand=subset(minus.exon,minus.exon$V10=='+')

#stats bin in 300 bins for both promary and offstrand using log transformed values
minus.primary.bin=stats.bin(minus.primary$rel.dist,log(abs(minus.primary[,4])+1),N=300)
mpb=cbind(matrix(minus.primary.bin$centers,ncol=1),minus.primary.bin$stats["mean",])

minus.offstrand.bin=stats.bin(minus.offstrand$rel.dist,log(abs(minus.offstrand[,4])+1),N=300)
mob=cbind(matrix(minus.offstrand.bin$centers,ncol=1),minus.offstrand.bin$stats["mean",])

pob[,2]=-pob[,2]
mob[,2]=-mob[,2]

sense=pepb
sense[,2]=sense[,2]+mpb[,2]

antisense=pob
antisense[,2]=antisense[,2]+mob[,2]

out.table <- data.frame(cbind(sense,antisense[,2]))
colnames(out.table) <- c("Position", "Sense", "Antisense")
out.table$SampleName <- Sample

#write.csv(out.table, 'average_coverage.csv')
write.csv(out.table, paste0(outFolder, "/",Sample, '_average_coverage.csv'))

pdf(paste0(outFolder, "/",Sample, '_gene_coverge_plot_WC.pdf'),h=10,w=12)
#pdf("test_strands.pdf",h=10,w=12)
plot(x=NULL,y=NULL,xlim=c(-1000,2000),ylim=c(-5,5), main=Sample)
lines(ppb,col=1,lwd=2)
lines(pob,col=2,lwd=2)
lines(mpb,col=3,lwd=2)
lines(mob,col=4,lwd=2)

abline(v=0,lty=2)
abline(v=1000,lty=2)
abline(h=0,lty=1,col='grey')
legend('topright',c('plus_primary','plus_offstrand','minus_primary','minus_offstrand'),lty=1,col=c(1,2,3,4))
dev.off()

pdf(paste0(outFolder, "/",Sample, '_gene_coverge_plot.pdf'),h=10,w=12)
#pdf("test_sense.pdf",h=10,w=12)
plot(x=NULL,y=NULL,xlim=c(-1000,2000),ylim=c(-5,5),main=Sample)
lines(sense,col=1,lwd=2)
lines(antisense,col=2,lwd=2)

abline(v=0,lty=2)
abline(v=1000,lty=2)
abline(h=0,lty=1,col='grey')
legend('topright',c('sense','antisense'),lty=1,col=c(1,2))
dev.off()
  }
