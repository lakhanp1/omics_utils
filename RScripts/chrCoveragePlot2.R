covPlot <- function(dirPath, extension){
 
  dataFiles <- list.files(dirPath,pattern = extension,full.names = T , recursive = F)
  png(filename = "E:/Chris_UM/Analysis/3_cAlbicanceVirulance/coverageFiles/BRAC1_coverage.png", width=5000,height=500*length(dataFiles)[1])
  col = "blue";
  
  par(mfrow = c(length(dataFiles),1))
  
  for(covFile in  dataFiles){
    print (covFile)
    
    if(file.info(covFile)$size!=0){
      cov <- read.table(covFile ,colClasses=c("character","numeric","numeric","numeric"))
    }
    
    start <- min(cov$V2)
    end <- max(cov$V3)
    
    plot(c(0, end-start), c(0, 500), xlab = "Chr Length", ylab = "Number of Reads",axes = FALSE)
    
    title(c(unlist(strsplit(basename(covFile),"_"))[1]),font.main= 10, col.main= "blue",cex.main = 5)
    axis(1, at=seq(0, end-start, by=100),seq(start, end, by=100) ,tick = T, hadj = 'vertical')
    axis(2, at=seq(0,500,50),labels=T,tick = T )
    
    # Draw peaks
    i<-1
    
    #if(strand == "-"){
      while(i <= nrow(cov)){
        rect(xleft = end-cov[i,3] ,ybottom = 0, ytop = cov[i,4], xright = end-cov[i,2],col = col, border = col)
        i <- i+1;
      }
    #}
    
    
     
  }
  
}