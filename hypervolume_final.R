####1. Load packages
library(hypervolume) #For hypervolume operations
library(cluster) #for Calculation of Gower distances

####2. Functions 
#Input prep (not run, saved it just to have a record of the workflow)
{
  lst <- list.files("D:\\PhD\\ModellingViviparity\\FullDist_Sal_Lsal\\MV_analyses\\hypervolume", pattern = "\\.csv$", full.names = T)
  input <- read.csv(lst[1], sep = ";")
  for (i in 2:length(lst))
  {
    input2 <- read.csv(lst[i], sep = ";")
    colnames(input2) <- colnames(input)
    input <- rbind(input, input2)
  } 
  rm(input2, lst, i)
  #eliminate unused variables and rows with NAs
  #input <- input[,c(1:3,5,9,10,19,20,22,23:25)]
  input$Kbin_1[input$Kbin_1 == -9999] <- 0
  input$Kbin_1 <- as.factor(input$Kbin_1)
  input[input==-9999] <- NA
  input <- na.omit(input)
  #Create new column with broad Larv/Viv categories
  
  for(i in 1: length(input[, 25]))
  {
    if (input[i,25] == "Larv")
    {
      input[i, 26] <- "Larv"
    }
    else
      input[i, 26] <- "Viv"
  }
  input[,26] <- as.factor(input[,26])
  #Assign viviparous S. sal and S. alg to distinct groups (we were combining them initially)
  input$mode <- as.character(input$mode)
  for(i in 1: length(input[, 25]))
  {
    if (input[i,25] == "Viv_salalg" & input[i,24] == "Salg")
    {
      input[i, 25] <- "Viv_Salg"
    }
    else
      if (input[i,25] == "Viv_salalg" & input[i,24] == "Ssal")
      {
        input[i, 25] <- "Viv_Ssal"
      }
  }
  input$mode <- as.factor(input$mode)
  #Assign subgroups for sp comparisons
  
  for(i in 1: length(input[, 24]))
  {
    if (input[i,24] == "Lycia")
    {
      input[i, 27] <- "Lycia"
    }
    else
      if (input[i,24] == "Sinf")
      {
        input[i, 27] <- "Sinf"
      }
    else
      if (input[i,24] == "Satra" | input[i,24] == "Slan")
      {
        input[i,27] <- "Alps"
      }
    else
    {
      input[i,27] <- "Selse"
    }
  }
  input$V27 <- as.factor(input$V27)
  #Check if everything is ok. You should have no NAs, and categorical predictors should be classed as factorial with all relevant levels present
  summary(input)
  colnames(input) <- c(colnames(input[1:25]), "rep", "subgroup")
  write.table(input, file = ".\\input.csv", quote = F, sep = ";", dec = ".")
}
outputs <- function(nm, reps)
{
  dist <- vector(mode = "numeric", length = length(nm))
  names(dist) <- nm
  volumes <- matrix(NA, nrow = length(nm), ncol = 6, byrow = T)
  rownames(volumes) <- nm
  colnames(volumes) <- c("hv1", "hv2", "Intersection", "Union", "Unique_hv1", "unique_hv2")
  overlap <- matrix(NA, nrow = length(nm), ncol = 5, byrow = T)
  rownames(overlap) <- nm
  colnames(overlap) <- c("SorObs", "pVal_SorObs", "Smax", "OI", "pVal_OI")
  perm.list <- list()
  for (i in 1:length(nm))
  {
    perm.list[[i]] <- matrix(NA, nrow = reps+1, ncol = 3, byrow = T)
  }
  list(volumes = volumes, dist = dist, overlap = overlap, perm.list = perm.list)
}
#Sorenson
sorenson.coef <- function(X, Y, I) {
  ## Calculates the Soerenson coefficient
  ## X - number of species in sample X
  ## Y - number of species in sample Y
  ## I - number of species common to both samples
  ## NOTE: the formula is equivalent to (2a)/(2a+b+c) but avoids to
  ## include the intersection at the denominator and it is simpler
  ## to calculate the maximum Sorenson.
  (2*I)/(X+Y)
}#X, Y and I are hv1, hv2 and Intersection, respectively
#OI 
OI <- function(sp1, sp2, nm, reps = rep,
               repsperpoint = 2000, verbose = TRUE, ...) 
{
  ### Test hv siginificance with Soerenson and overlap index
  
  
  ## Build the observed hypervolumes for each species
  hv.sp1 <- hypervolume(sp1, method = "svm",
                        samples.per.point = repsperpoint)
                        
  
  hv.sp2<- hypervolume(sp2, method = "svm",
                       samples.per.point = repsperpoint)
  
  hv.set <- hypervolume_set(hv.sp1, hv.sp2, check.memory=FALSE,
                            verbose = verbose > 1)
  #Save hypervolume plots to disk
  pdf(file = paste("./plots/", nm[i], ".pdf", sep = ""))
  plot(hv.set)
  dev.off()
  #Compute volumes and distances
  volumes <- get_volume(hv.set)
  dist <- hypervolume_distance(hv.sp1, hv.sp2, check.memory = F)
  ## merge both species for a common niche/null distribution
  sp <- rbind(sp1, sp2)
  
  sp1.n <- nrow(sp1)
  sp2.n <- nrow(sp2)
  sp.n <- nrow(sp)
  
  null.dist <- matrix(NA, reps+1,3) #create matrix to store null dist. # rows = # permutations +1; # columns = 3 (names below)
  colnames(null.dist) <- c("Soerensen", "MaxSorensen", "OverlapIndex")
  for (rep in 1:reps) {
    if (verbose) cat(sprintf("%.1f", rep/reps*100), "%\r")
    
    ## Sample the common niche 
    samples <- sample(1:sp.n, sp1.n)
    sp1.rep <- sp[samples,]
    sp2.rep <- sp[-samples,]

    sp1.hv.rep <- hypervolume(sp1.rep, method = "svm",
                              samples.per.point = repsperpoint,
                              verbose = verbose > 1, ...)
    
    sp2.hv.rep <- hypervolume(sp2.rep, method = "svm",
                              samples.per.point = repsperpoint,
                              verbose = verbose > 1, ...)
    
    
    msg <- capture.output(
      hv.rep.set <- hypervolume_set(sp1.hv.rep, sp2.hv.rep,
                                    verbose = FALSE,
                                    check.memory=FALSE))
    
    v.rep <- get_volume(hv.rep.set)
    
    null.dist[rep, 1] <- sorenson.coef(v.rep[1], v.rep[2], v.rep[3]) #Sorensen for each rep
    null.dist[rep, 2] <- sorenson.coef(v.rep[1], v.rep[2],
                                       min(v.rep[1:2])) #Max Sorensen for each rep
  }
  
  
  S.obs <- sorenson.coef(volumes[1], volumes[2], volumes[3])
  null.dist[reps+1, 1] <- S.obs #Observed Sorensen
  
  Smax.obs <- sorenson.coef(volumes[1], volumes[2],
                            min(volumes[1:2]))
  null.dist[reps+1, 2] <- Smax.obs #Observed Max Sorensen
  
  null.dist[,3] <- null.dist[,1] / null.dist[,2] #OI for all permutations and observed (calculated as Sorensen/Max. Sorensen)
  
  OI.obs <- null.dist[reps+1, 3]
  
  p.value.sor <- sum(null.dist[,1] <= S.obs) / (reps+1)
  p.value.oi <- sum(null.dist[,3] <= OI.obs) / (reps+1)

  overlap <- c(S.obs,p.value.sor,Smax.obs,OI.obs, p.value.oi)
  
  ### plot Sorensen and null distribution
  pdf(file = paste("./plots/null_sor_", nm[i], ".pdf", sep = ""))
  hist(null.dist[,1], col='gray', breaks=50,
       xlab="Sorensen index", main = "")
  arrows(null.dist[reps+1,1], 100, null.dist[reps+1,1], 0, col='red')
  dev.off()
  ### plot Overlap Index and null distribution
  pdf(file = paste("./plots/null_OI_", nm[i], ".pdf", sep = ""))
  hist(null.dist[,3], col='gray', breaks=50,
       xlab="Overlap index", main="")
  arrows(null.dist[reps+1,3], 100, null.dist[reps+1,3], 0, col='red')
  dev.off()
  
  list(volumes = volumes, dist = dist, overlap = overlap, null.dist = null.dist)
}
#Number of permutations for Sorensen/OI
rep = 999 

####3. Load input
input <- read.csv(".\\input.csv", sep = ";")
input <- input[,c(5,9,10,19,20,22,24:27)]
pcoa <- cmdscale(daisy(scale(input[,1:6])),k =3, eig = T, list. = T)

###4. Among sp
#Create empty objects to receive outputs
nm.sp <- c("sp_Lyc_v_salall", "sp_Sinf_v_Selse", "sp_Scor_v_Alps", "sp_Satra_v_Slan", "sp_Ssal_v_Salg")
sp <- outputs(nm.sp, rep)
#Create subsets of data for pairwise analyses
sp1.list <- list(subset(pcoa$points, input$sp == "Lycia"), subset(pcoa$points, input$sp == "Sinf"), subset(pcoa$points, input$sp == "Scor"), subset(pcoa$points, input$sp == "Satra"), subset(pcoa$points, input$sp == "Ssal"))
sp2.list <- list(subset(pcoa$points, input$sp != "Lycia"), subset(pcoa$points, input$subgroup == "Selse"), subset(pcoa$points, input$subgroup == "Alps"), subset(pcoa$points, input$sp == "Slan"), subset(pcoa$points, input$sp == "Salg"))

for (i in 1:length(nm.sp))
{
  sp1 <- sp1.list[[i]]
  sp2 <- sp2.list[[i]]
  temp <- OI(sp1, sp2, nm.sp, rep)
  sp$volumes[i,] <- temp$volumes
  sp$dist[i] <- temp$dist
  sp$overlap[i,] <- temp$overlap
  sp$perm.list[[i]] <- temp$null.dist
  rm(temp)
}
frac.hv1 <-(sp$volumes[,5]/sp$volumes[,1])*100
frac.hv2 <-(sp$volumes[,6]/sp$volumes[,2])*100
sp$volumes <- cbind(sp$volumes, frac.hv1, frac.hv2)
###5. Larviv (larviparous vs. pueriparous)
#Create empty objects to receive outputs (only the first was kept for publication)
nm.larviv <- c("larviv_larv_v_viv", "larviv_larv_v_alps", "larviv_larv_v_lyc", "larviv_larv_v_ssal", "larviv_larv_v_salg")
larviv <- outputs(nm.larviv, rep)
#Create subsets of data for pairwise analyses
sp1.list <- list(subset(pcoa$points, input$mode == "Larv"), subset(pcoa$points, input$mode == "Larv"), subset(pcoa$points, input$mode == "Larv"), subset(pcoa$points, input$mode == "Larv"), subset(pcoa$points, input$mode == "Larv"))
sp2.list <- list(subset(pcoa$points, input$rep == "Viv"), subset(pcoa$points, input$mode == "Viv_Alp"), subset(pcoa$points, input$mode == "Viv_Lycia"), subset(pcoa$points, input$mode == "Viv_Ssal"), subset(pcoa$points, input$mode == "Viv_Salg"))

for (i in 1:length(nm.larviv))
{
  sp1 <- sp1.list[[i]]
  sp2 <- sp2.list[[i]]
  temp <- OI(sp1, sp2, nm.larviv, rep)
  larviv$volumes[i,] <- temp$volumes
  larviv$dist[i] <- temp$dist
  larviv$overlap[i,] <- temp$overlap
  larviv$perm.list[[i]] <- temp$null.dist
  rm(temp)
}

frac.hv1 <-(larviv$volumes[,5]/larviv$volumes[,1])*100
frac.hv2 <-(larviv$volumes[,6]/larviv$volumes[,2])*100
larviv$volumes <- cbind(larviv$volumes, frac.hv1, frac.hv2)
###6. Viviv (among pueriparous groups)
#Create empty objects to receive outputs
nm.viviv <- c("viviv_Ssal_v_Salg", "viviv_Ssal_v_Alps", "viviv_Ssal_v_lyc", "viviv_Salg_v_Alps", "viviv_Salg_v_Lyc","viviv_Alps_v_Lyc")
viviv <- outputs(nm.viviv, rep)
#Create subsets of data for pairwise analyses
sp1.list <- list(subset(pcoa$points, input$mode == "Viv_Ssal"), subset(pcoa$points, input$mode == "Viv_Ssal"), subset(pcoa$points, input$mode == "Viv_Ssal"), subset(pcoa$points, input$mode == "Viv_Salg"), subset(pcoa$points, input$mode == "Viv_Salg"), subset(pcoa$points, input$mode == "Viv_Alp"))
sp2.list <- list(subset(pcoa$points, input$mode == "Viv_Salg"), subset(pcoa$points, input$mode == "Viv_Alp"), subset(pcoa$points, input$mode == "Viv_Lycia"), subset(pcoa$points, input$mode == "Viv_Alp"), subset(pcoa$points, input$mode == "Viv_Lycia"), subset(pcoa$points, input$mode == "Viv_Lycia"))

for (i in 1:length(nm.viviv))
{
  sp1 <- sp1.list[[i]]
  sp2 <- sp2.list[[i]]
  temp <- OI(sp1, sp2, nm.viviv, rep)
  viviv$volumes[i,] <- temp$volumes
  viviv$dist[i] <- temp$dist
  viviv$overlap[i,] <- temp$overlap
  viviv$perm.list[[i]] <- temp$null.dist
  rm(temp)
}
frac.hv1 <-(viviv$volumes[,5]/viviv$volumes[,1])*100
frac.hv2 <-(viviv$volumes[,6]/viviv$volumes[,2])*100
viviv$volumes <- cbind(viviv$volumes, frac.hv1, frac.hv2)

rm(frac.hv1, frac.hv2)


#########PCOA##########
#This section deals with plotting an approximation of the visual representation of hypervolumes by plotting the density function around each group along the axes of the PCoA used as input for hypervolume analysis

pts <- pcoa$points
df <- as.data.frame(cbind(pts,input$mode))
colnames(df) <- c("x", "y", "z", "mode")
df$mode <- as.factor(df$mode)
levels(df$mode) <- levels(input$mode)
gg <- merge(df,aggregate(cbind(mean.x=x,mean.y=y, mean.z=z)~mode,df,mean),by="mode")

#
#axes 1-2
centroids <- aggregate(cbind(x,y)~mode,df,mean)
f <- function(z) qt(0.025,df=length(z)-1, lower.tail=F)* sd(z)/sqrt(length(z))  # function to calculate 95% CI
se        <- aggregate(cbind(se.x=x,se.y=y)~mode,df,f)
centroids <- merge(centroids,se, by="mode")    # add CI column to centroids
ggplot(gg, aes(x,y,color=factor(mode)))+
  geom_point(data=centroids, size=5)+
  geom_errorbar(data=centroids,aes(ymin=y-se.y,ymax=y+se.y),width=0.1)+
  geom_errorbarh(data=centroids,aes(xmin=x-se.x,xmax=x+se.x),height=0.1)

#axes 1-3
centroids <- aggregate(cbind(x,z)~mode,df,mean)
f <- function(z) qt(0.025,df=length(z)-1, lower.tail=F)* sd(z)/sqrt(length(z))  # function to calculate 95% CI
se        <- aggregate(cbind(se.x=x,se.y=z)~mode,df,f)
centroids <- merge(centroids,se, by="mode")    # add CI column to centroids
ggplot(gg, aes(x,z,color=factor(mode)))+
  geom_point(data=centroids, size=5)+
  geom_errorbar(data=centroids,aes(ymin=z-se.y,ymax=z+se.y),width=0.1)+
  geom_errorbarh(data=centroids,aes(xmin=x-se.x,xmax=x+se.x),height=0.1)

#axes 2-3
centroids <- aggregate(cbind(y,z)~mode,df,mean)
f <- function(z) qt(0.025,df=length(z)-1, lower.tail=F)* sd(z)/sqrt(length(z))  # function to calculate 95% CI
se        <- aggregate(cbind(se.x=y,se.y=z)~mode,df,f)
centroids <- merge(centroids,se, by="mode")    # add CIcolumn to centroids
ggplot(gg, aes(y,z,color=factor(mode)))+
  geom_point(data=centroids, size=5)+
  geom_errorbar(data=centroids,aes(ymin=z-se.y,ymax=z+se.y),width=0.1)+
  geom_errorbarh(data=centroids,aes(xmin=y-se.x,xmax=y+se.x),height=0.1)


##3d PCoA
#library(pca3d)
#library(rgl)
#df.sub <- as.matrix(subset(df[,1:3], df$mode != "Larv"))
#as.factor(subset(df$mode, df$mode != "Larv"))

#colViviv <- c("#000000","#F50ABE", "#0A74F5", "#F50A0A") #codes derived from filling in function rgb() with rgb values from ArcGIS (for color homogeneity across figures). The order corresponds to the alphabetical order of the groups 

#pca3d(df.sub, group=mode.sub, palette = colViviv, axes.color = "black", show.plane = FALSE, bg = "white", radius = 0.5, show.ellipses = TRUE, ellipse.ci = 0.95, shape=1, show.shapes=TRUE, show.centroids = F, new = F)
#rgl.snapshot(".\\pcoa_viviv3.png") #Saving to file. 

#pca3d(pcoa$points, group=input$rep, axes.color = "black", show.plane = FALSE, bg = "white", radius = 0.5, show.ellipses = TRUE, ellipse.ci = 0.95, shape=1, show.shapes=TRUE, show.centroids = F)
#rgl.snapshot(".\\pcoa_larviv.png") #Saving to file. 


##3dPCoA part 2
library(scatterplot3d)
df.sub <- as.matrix(subset(df[,1:4], df$mode != "Larv"))
as.factor(subset(df$mode, df$mode != "Larv"))


colors <- c(NA, "#000000","#0AF5F5", "#0A74F5", "#F50A0A")
colors <- colors[as.numeric(mode.sub)]
viviv.s3d <- scatterplot3d(df.sub[,1:3], pch = 16, xlim = c(-6,6), ylim = c(-6,3), zlim = c(-3, 4), color=colors, angle = 120, box = F, main="PCoA pueriparous groups",
              xlab = "Dim 1",
              ylab = "Dim 2",
              zlab = "Dim 3")
legend(viviv.s3d$xyz.convert(-9.3,-1,5), legend = levels(mode.sub),
       col =  c(NA, "#000000","#0AF5F5", "#0A74F5", "#F50A0A"), pch = 16)

colarviv <- c("#E3E38D","#5331c7")
colarviv <- colarviv[as.numeric(input$rep)]
larviv.s3d <- scatterplot3d(pcoa$points, pch = 20, xlim = c(-6,6), ylim = c(-6,3), zlim = c(-3, 4), color=colarviv, angle = 120, box = F, main="PCoA reproductive modes",
              xlab = "Dim 1",
              ylab = "Dim 2",
              zlab = "Dim 3")
legend(larviv.s3d$xyz.convert(-9.3,-1,5), legend = levels(input$rep),
       col =  c("#E3E38D","#5331c7"), pch = 16)

##Aux stuff (not run)
#plot(1:length(pcoa$eig[1:10]),pcoa$eig[1:10],'bo-'); #aux code for running a scree plot (to visualize where the eigenvalue curve smooths over and therefore how many axes you need to properly represent your dataset.In this dataset it would eb the first 6 or 7)

