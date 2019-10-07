#====================================
# utilities.R
#====================================

### function

topGOterms = function( fg.genes = NULL,
                       bg.genes = NULL,
                       organism = "Mouse", 
                       ontology.use = "BP",
                       stats.use = "fisher",
                       algorithm.use = "weight01",
                       topnodes.print=20,
                       num.char=100){
  
  if (is.null(fg.genes) | is.null(bg.genes)){
    stop("Error : Both gene lists are empty")
  }
  
  require(topGO)
  if (organism == "Mouse"){
    mapping.use = "org.Mm.eg.db"
    library(org.Mm.eg.db)
  } else if (organism == "Human"){
    mapping.use = "org.Hs.eg.db"
    library(org.Hs.eg.db)
  } else {
    stop("Error : Organisms other than mouse not supported currently")
  }
  
  n = length(bg.genes)
  geneList = integer(n)
  names(geneList) = bg.genes
  geneList[intersect(names(geneList), fg.genes)]=1
  print(paste0("Total ", length(geneList), " genes. ", sum(geneList), " genes in the foreground"))
  geneList = factor(geneList)
  
  if (ontology.use %in% c("BP", "CC", "MF")){
    print(paste0("Using Ontology : ", ontology.use))
  } else {
    stop("Error: Ontology not available. Should be one of BP, CC or MF")
  }
  # Make GO object
  GOdata <- new("topGOdata",
                description = "GOanalysis",
                ontology = ontology.use,
                allGenes = geneList,
                annot = annFUN.org,
                mapping = mapping.use,
                ID = "SYMBOL",
                nodeSize = 10)
  print(paste0("Using the ", stats.use, " statistic with the ", algorithm.use, " algorithm"))
  res.result <- runTest(GOdata, statistic = stats.use, algorithm = algorithm.use)
  to.return = list()
  to.return$GOdata = GOdata
  to.return$res.table <- GenTable(GOdata, pval = res.result, topNodes = topnodes.print, numChar = num.char)
  return(to.return)
}


if(F){ #=== Begin
# Plot confusion matrix
plotConfusionMatrix = function(X,row.scale=TRUE, col.scale=FALSE, col.low="blue", col.high="red", max.size=5, ylab.use="Known", xlab.use="Predicted", order=NULL, x.lab.rot=FALSE, plot.return=TRUE){
  
  if (!col.scale & row.scale){ X = t(scale(t(X), center=FALSE, scale=rowSums(X)));  X=X*100 }
  if (col.scale & !row.scale){ X = scale(X, center=FALSE, scale=colSums(X)); X = X*100 }
  if(col.scale & row.scale){
    print("Only one of row.scale or col.scale should be true. performing row scaling by default")
    X = t(scale(t(X), center=FALSE, scale=rowSums(X)))
    X=X*100
  }
  X[is.na(X)] = 0
  if (max(X) > 100){
    X=X/100
  }
  
  orig.rownames = rownames(X)
  orig.colnames = colnames(X)
  
  if (!is.null(order)){
    if (order == "Row"){  
      factor.levels = c()
      for (i1 in colnames(X)){
        if (max(X[,i1]) < 50) next
        ind.sort = rownames(X)[order(X[,i1], decreasing=TRUE)]
        ind.sort = ind.sort[!(ind.sort %in% factor.levels)]
        factor.levels = c(factor.levels, ind.sort[1])
      }
      factor.levels = c(factor.levels, setdiff(rownames(X), factor.levels))
      factor.levels = factor.levels[!is.na(factor.levels)]
    } 
    
    if (order == "Col") {
      factor.levels = c()
      for (i1 in rownames(X)){
        if (max(X[i1,]) < 50) next
        ind.sort = rownames(X)[order(X[i1,], decreasing=TRUE)]
        ind.sort = ind.sort[!(ind.sort %in% factor.levels)]
        factor.levels = c(factor.levels, ind.sort[1])
      }
      factor.levels = c(factor.levels, setdiff(rownames(t(X)), factor.levels))
      factor.levels = factor.levels[!is.na(factor.levels)]
    } 
  } else {
    factor.levels = rownames(t(X))
  }
  
  factor.levels = c(factor.levels, setdiff(rownames(X), factor.levels))
  X = melt(X)
  colnames(X) = c("Known", "Predicted", "Percentage")
  #X$Known = factor(X$Known, levels=rev(unique(X$Known)));
  #X$Predicted = factor(X$Predicted, levels = rev(factor.levels))
  
  if (!is.null(order)){
    if (order == "Row"){ 
      X$Known = factor(X$Known, levels=rev(factor.levels));
      X$Predicted = factor(X$Predicted, levels = orig.colnames)
      
    }
    if (order == "Col"){
      X$Predicted = factor(X$Predicted, levels = factor.levels);
      X$Known = factor(X$Known, levels=rev(orig.rownames));
    }
  } else {
    X$Known = factor(X$Known, levels=rev(unique(X$Known)));
    X$Predicted = factor(X$Predicted, levels=unique(X$Predicted));
  }
  
  
  
  
  p = ggplot(X, aes(y = Known,  x = Predicted)) + geom_point(aes(colour = Percentage,  size =Percentage)) + 
    scale_color_gradient(low =col.low,   high = col.high, limits=c(0, 100 ))+scale_size(range = c(1, max.size))+   theme_bw() #+nogrid
  p = p + xlab(xlab.use) + ylab(ylab.use) + theme(axis.text.x=element_text(size=12, face="italic", hjust=1)) + 
    theme(axis.text.y=element_text(size=12, face="italic"))  
  
  if (x.lab.rot) {
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  }
  print(p)
  
  if (plot.return) return(p)
}


plot_sample_dist = function(object, id.use = "orig.ident",max.val.perc=100,
                            max.size=10,row.scale=FALSE,top.mar=1, left.mar = 4.5,
                            right.mar=1.5, bottom.mar = -1,do.transpose=FALSE,
                            ident.order=NULL,x.lab.rot=0,to.return=TRUE,alpha.use=1) {
  
  if (is.null(ident.order)){
    ident.order = levels(object@ident)
  }
  
  
  # Get sample distribution
  ExpMat = table(object@meta.data[,id.use],object@ident)
  ExpMat = ExpMat[,ident.order]
  if(row.scale){
    ExpMat = t(scale(t(ExpMat), center=FALSE, scale=rowSums(ExpMat)))*100
  } else {
    ExpMat = scale(ExpMat, center=FALSE, scale=colSums(ExpMat))*100
  }
  ExpVal = melt(ExpMat)
  colnames(ExpVal) = c("sample","cluster","perc")
  ExpVal$sample = factor(ExpVal$sample, levels=rev(levels(object@meta.data[,id.use])))
  ExpVal$cluster = factor(ExpVal$cluster, levels= ident.order)
  p=ggplot(ExpVal, aes(y = factor(sample),  x = factor(cluster))) + geom_point(aes(colour = perc,  size =perc),alpha=alpha.use) + 
    scale_color_gradient(low ="white",   high = "darkblue", limits=c(min(ExpVal$perc), max(ExpVal$perc) ))+scale_size(range = c(1, max.size))+   theme_bw() +nogrid
  p = p + xlab("Cluster") + ylab("Sample") + theme(axis.text.x=element_text(size=12, face="italic", hjust=1)) + 
    theme(axis.text.y=element_text(size=12, face="italic")) + theme(axis.text.x = element_text(angle = x.lab.rot))
  
  if (do.transpose){
    ExpVal$sample = factor(ExpVal$sample, levels=levels(object@meta.data[,id.use]))
    p=ggplot(ExpVal, aes(y = factor(cluster),  x = factor(sample))) + geom_point(aes(colour = perc,  size =perc), alpha=alpha.use) +
      scale_color_gradient(low ="white",   high = "darkblue", limits=c(min(ExpVal$perc), max(ExpVal$perc) )) + 
      scale_size(range = c(1, max.size))+   theme_bw() +nogrid
    p = p + xlab("Sample") + ylab("Cluster") + theme(axis.text.x=element_text(size=12, face="italic", hjust=1)) + 
      theme(axis.text.x=element_text(size=12, face="italic", angle=45, hjust=1))
  }
  
  print(p)
  
  if (to.return){
    return(p)
  }
  
}

nogrid=theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
sort.column=function(x, col) {
  return(x[order(x[,col]),])
}


RF_train = function(train_object0, var.genes, do.scale=FALSE){
  
  predictor_Data = as.matrix(train_object0@data[var.genes,])
  if (do.scale) predictor_Data = t(scale(t(predictor_Data)))
  max.cells.per.ident = 700; train.frac = 0.6
  training.set = c(); validation.set=c()
  training.label = c(); validation.label=c();
  print(paste0("Using mininum of ", 0.5*100, " percent cells or ", max.cells.per.ident, " cells per cluster for training"))
  for (i in as.character(levels(train_object0@ident))){
    cells.in.clust = WhichCells(train_object0,i);
    n = min(max.cells.per.ident, round(length(cells.in.clust)*train.frac))
    train.temp = cells.in.clust[sample(length(cells.in.clust))][1:n]
    validation.temp = setdiff(cells.in.clust, train.temp)
    training.set = c(training.set,train.temp); validation.set=c(validation.set,validation.temp)
    training.label = c(training.label, rep(as.numeric(i)-1,length(train.temp))); validation.label = c(validation.label, rep(as.numeric(i)-1, length(validation.temp)));
  }
  
  print(length(training.set))
  tmp = as.vector(table(training.label))
  sampsizes = rep(min(tmp),length(tmp))
  rf = randomForest(x=t(predictor_Data[,training.set]), y=factor(training.label), importance = TRUE, ntree = 301, proximity=TRUE, sampsize=sampsizes, keep.inbag=TRUE, replace=FALSE) 
  validation_pred_rf = predict(rf,t(predictor_Data[,validation.set]))
  #plotConfusionMatrix(table(validation.label, validation_pred_rf))
  
  return(rf)
}


} 
#==END