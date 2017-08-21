
####Yan Xia -Arale Hsia

setwd("/zs32/home/yxia/01sexmethy/01rush/20161009WGCNAexpression")
#save(datat,dataTraits,file="00rushdatat_dataTraits.RData")

load(file="00rushdatat_dataTraits.RData")

datat[1:5,1:5]
dataTraits[1:5,]
summary(as.factor(dataTraits$sex))
#  0   1 
#403 230 


#########
net<- blockwiseModules(datat,
 	power = 4,
 	networkType = "unsigned",
 	TOMType ="signed",
    corType="bicor",
 	deepSplit=4,
	maxBlockSize=24660, 
 	minModuleSize =50,
	reassignThreshold = 0, 
	mergeCutHeight = 0.15,
 	numericLabels = TRUE,
 	pamStage=FALSE,
 	pamRespectsDendro = FALSE,
	saveTOMs = TRUE,
	saveTOMFileBase ="rushexpall",
	verbose = 3);

c("######return rush table net$colors #########")
table(net$colors)

pdf(file = "07rushexpr633all_clusters.pdf", width = 12, height = 9)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(net, file = "07Net-rushmethy.RData")


allinfor<-data.frame(ID=colnames(datat),module=net$colors,color=labels2colors(net$colors))
dim(allinfor)

c("#######write.csv rush methy probes infor")
write.csv(allinfor,file="07rushexpr_probemodule.csv")

######### Step 3  gene significance and module membership
c("#######plot matrix of module and trait")
nGenes<- ncol(datat)
nSamples<- nrow(datat)
MEs0 <- moduleEigengenes(datat, moduleColors)$eigengenes
MEs=orderMEs(MEs0)
moduleTraitCor=cor(MEs,dataTraits,use="p")
moduleTraitPvalue=corPvalueStudent(moduleTraitCor,nSamples)

pdf("008rushall_module_traitrelationship.pdf",width=12,height=10)
textMatrix=paste(signif(moduleTraitCor,2),"\n(",signif(moduleTraitPvalue,1),")",sep="");
dim(textmatrix)=dim(moduleTraitCor)
par(mar=c(6,12,3,3))

labeledHeatmap(Matrix=moduleTraitCor,
			   xLabels=names(dataTraits),
			   yLabels=names(MEs),
			   ySymbols=names(MEs),
			   colorLabels=FALSE,
			   colors=greenWhiteRed(50),
			   textMatrix=textMatrix,
			   setStdMargins=FALSE,
			   cex.text=0.5,
			   zlim=c(-1,1),
			   main=paste("Rush Module-trait relationships"))
dev.off()












