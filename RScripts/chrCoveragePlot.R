dr <- function(dirPath, extension){
  strand <- "-"
  
  dataFiles <- list.files(dirPath,pattern = extension,full.names = T , recursive = F)
  
  
  png(filename = "E:/Chris_UM/Analysis/CoreData/5_MiaoKaiSC/coverage/run1/BRCA1_run1_coverage.png", width=5000, height=500*length(dataFiles)[1], res = 100)
  #par(mfrow = c())
  col = "blue";
  
  #This parameter decides the number of rows and columns in a given plot PNG file
  par(mfrow = c(length(dataFiles),1))
  
  
  for(covFile in  dataFiles){
    
    print (covFile)
    
    if(file.info(covFile)$size!=0){
      cov <- read.table(covFile ,colClasses=c("character","numeric","numeric","numeric"))
    }
    
    start <- min(cov$V2)
    end <- max(cov$V3)
    
    plot(c(0, end-start+1000), c(0, 60), xlab = "Chr Length", ylab = "Number of Reads",axes = FALSE)
    
    title(c(unlist(strsplit(basename(covFile),"_"))[1]),font.main= 10, col.main= "blue",cex.main = 5)
    
    #axis(1, at=seq(start, end, by=100),rev(seq(-2000, end-start-2000, by=100)) ,tick = T, hadj = 'vertical')
    axis(side = 1, at = seq(0, end-start+1000, by=100), labels = seq(-2000, end-start-1000, by=100) ,tick = T, hadj = 'vertical')
    axis(side = 1, at = c(1000,2000,end-start-500,end-start+500), labels = c("Promoter","Start","End","Downstream"), cex.axis=1, col.axis="red", padj = 2, tick = F)
    
    axis(2, at=seq(0,100,5),labels=T,tick = T )
    
    #draw promoter
    rect(xleft = 0 ,ybottom = -2, ytop = -0.5, xright = 2000,col = "red", border = col)
    
    #Shown GeneDownStream
    rect(xleft = end-start-500 ,ybottom = -2, ytop = -0.5, xright = end-start ,col = "red", border = col)
    
    lines(c(2000,2000),c(-10,0),col="red",lwd = 2)
    lines(c(end-start-500, end-start-500),c(0,-5),col="red",lwd = 2)
    
    # Draw peaks
    i<-1
    
    if(strand == "-"){
      while(i <= nrow(cov)){
        rect(xleft = end-cov[i,3] ,ybottom = 0, ytop = cov[i,4], xright = end-cov[i,2],col = col, border = col)
        i <- i+1;
      }
    }
    
    # Draw exons
    exons <- read.table("E:/Chris_UM/Analysis/CoreData/5_MiaoKaiSC/coverage/run1/exons.txt", colClasses=c("numeric","numeric"))
    i<-1
    
    while(i <= nrow(exons)){
      rect(xleft = end-exons[i,2] ,ybottom = -2, ytop = 0, xright = end-exons[i,1] ,col = "green", border = col)
      i <- i+1;
    }
    
  }
  
  # print(warnings())
  dev.off()
}

