rm(list=ls())
setwd("D:/JanHeng_Rept1")

filesIwant <- dir(pattern="_cell_seg_data.txt") #13 images containing 2 sections of tumors each, in field of views (multiple fields -> 1 section, 2 sections -> 1 slide) 

# these are the unique 13 images ID
filenames <- unique(sapply(strsplit(filesIwant, split='[', fixed=TRUE), function(x) (x[[1]])))
filenames_all <- sapply(strsplit(filesIwant, split='[', fixed=TRUE), function(x) (x[[1]]))

# "Lin_Ctrl1_Scan1_" has 801 fields of view.
idx <- which(filenames_all %in% filenames[1])

#now to break the fields of view into the 2 sections
a <- filesIwant[idx]
a1 <- sapply(strsplit(a, split=']', fixed=TRUE), function(x) (x[[1]]))
a2 <- as.numeric(sapply(strsplit(a1, split=',', fixed=TRUE), function(x) (x[[2]])))
plot(a2)
mean(a2)
idx1 <- which(a2 < mean(a2)) #these are fields of view from the first section
idx2 <- which(a2 > mean(a2)) #these are fields of view from the second section

topsection <- filesIwant[idx1]
open <- read.csv(topsection[1], sep="\t")
open1 <- open[which(open$Phenotype == "CD4+" | open$Phenotype == "CD8+"),] #only want CD4 or CD8
colnames(open1)
open2 <- open1[,c(1:4,63,68,73,78)] #cytoplasmic intensity is very similar to entire cell intensity, just pick cytoplasmic

for (j in 2:length(topsection)){
  open <- read.csv(topsection[j], sep="\t")
  open1 <- open[which(open$Phenotype == "CD4+" | open$Phenotype == "CD8+"),] #only want CD4 or CD8
  colnames(open1)
  open3 <- open1[,c(1:4,63,68,73,78)] #cytoplasmic intensity is very similar to entire cell intensity, just pick cytoplasmic
  open2 <- rbind(open2, open3)
}

bottomsection <- filesIwant[idx2]
openB <- read.csv(bottomsection[1], sep="\t")
openB1 <- openB[which(openB$Phenotype == "CD4+" | openB$Phenotype == "CD8+"),] #only want CD4 or CD8
openB2 <- openB1[,c(1:4,63,68,73,78)] #cytoplasmic intensity is very similar to entire cell intensity, just pick cytoplasmic

for (j in 2:length(bottomsection)){
  openB <- read.csv(bottomsection[j], sep="\t")
  openB1 <- openB[which(openB$Phenotype == "CD4+" | openB$Phenotype == "CD8+"),] #only want CD4 or CD8
  openB3 <- openB1[,c(1:4,63,68,73,78)] #cytoplasmic intensity is very similar to entire cell intensity, just pick cytoplasmic
  openB2 <- rbind(openB2, openB3)
}

open2$section <- "Top"
openB2$section <- "Bottom"

colnames(open2)
colnames(openB2)

comb <- rbind(open2, openB2) #merged into 1 sheet for 1 image with 2 sections
combfilename <- paste0(filenames[1], "CD4CD8only.csv")
write.csv(comb, file=combfilename, row.names = T, col.names = T)


#########################################
## repeat in loop to clean up the remaining 11 images

for (i in 2:12){
# "Lin_Ctrl1_Scan1_" has 801 fields of view.
idx <- which(filenames_all %in% filenames[i])

#now to break the fields of view into the 2 sections
a <- filesIwant[idx]
a1 <- sapply(strsplit(a, split=']', fixed=TRUE), function(x) (x[[1]]))
a2 <- as.numeric(sapply(strsplit(a1, split=',', fixed=TRUE), function(x) (x[[2]])))
plot(a2)
mean(a2)
idx1 <- which(a2 < mean(a2)) #these are fields of view from the first section
idx2 <- which(a2 > mean(a2)) #these are fields of view from the second section

topsection <- filesIwant[idx1]
open <- read.csv(topsection[1], sep="\t")
open1 <- open[which(open$Phenotype == "CD4+" | open$Phenotype == "CD8+"),] #only want CD4 or CD8
colnames(open1)
open2 <- open1[,c(1:4,63,68,73,78)] #cytoplasmic intensity is very similar to entire cell intensity, just pick cytoplasmic

for (j in 2:length(topsection)){
  open <- read.csv(topsection[j], sep="\t")
  open1 <- open[which(open$Phenotype == "CD4+" | open$Phenotype == "CD8+"),] #only want CD4 or CD8
  colnames(open1)
  open3 <- open1[,c(1:4,63,68,73,78)] #cytoplasmic intensity is very similar to entire cell intensity, just pick cytoplasmic
  open2 <- rbind(open2, open3)
}

bottomsection <- filesIwant[idx2]
openB <- read.csv(bottomsection[1], sep="\t")
openB1 <- openB[which(openB$Phenotype == "CD4+" | openB$Phenotype == "CD8+"),] #only want CD4 or CD8
openB2 <- openB1[,c(1:4,63,68,73,78)] #cytoplasmic intensity is very similar to entire cell intensity, just pick cytoplasmic

for (j in 2:length(bottomsection)){
  openB <- read.csv(bottomsection[j], sep="\t")
  openB1 <- openB[which(openB$Phenotype == "CD4+" | openB$Phenotype == "CD8+"),] #only want CD4 or CD8
  openB3 <- openB1[,c(1:4,63,68,73,78)] #cytoplasmic intensity is very similar to entire cell intensity, just pick cytoplasmic
  openB2 <- rbind(openB2, openB3)
}

open2$section <- "Top"
openB2$section <- "Bottom"

colnames(open2)
colnames(openB2)

comb <- rbind(open2, openB2) #merged into 1 sheet for 1 image with 2 sections
combfilename <- paste0(filenames[i], "CD4CD8only.csv")
write.csv(comb, file=combfilename, row.names = T, col.names = T)
}

## for PosCtrl2 -> only has 1 section
filenames[13]
idx <- which(filenames_all %in% filenames[13])

#now to break the fields of view into the 2 sections
topsection <- filesIwant[idx]

open <- read.csv(topsection[1], sep="\t")
open1 <- open[which(open$Phenotype == "CD4+" | open$Phenotype == "CD8+"),] #only want CD4 or CD8
colnames(open1)
open2 <- open1[,c(1:4,63,68,73,78)] #cytoplasmic intensity is very similar to entire cell intensity, just pick cytoplasmic

for (j in 2:length(topsection)){
  open <- read.csv(topsection[j], sep="\t")
  open1 <- open[which(open$Phenotype == "CD4+" | open$Phenotype == "CD8+"),] #only want CD4 or CD8
  colnames(open1)
  open3 <- open1[,c(1:4,63,68,73,78)] #cytoplasmic intensity is very similar to entire cell intensity, just pick cytoplasmic
  open2 <- rbind(open2, open3)
}

open2$section <- "Top"

comb <- open2 #only has 1 section
combfilename <- paste0(filenames[13], "CD4CD8only.csv")
write.csv(comb, file=combfilename, row.names = T, col.names = T)


#########################################
## Data analysis

rm(list=ls())
setwd("D:/JanHeng_Rept1")

filesIwant <- dir(pattern="CD4CD8only.csv") #13 images, each image has "top and bottom".

mtx <- matrix(NA, nrow=13, ncol=7)
mtx[,1] <- filesIwant
colnames(mtx) <- c("Image", "meanCD4counts", "meanCD8counts", "meaniFNCD4", "meaniFNCD8", "meanGxBCD4", "meanGxBCD8")

for (i in 1:12){
open <- read.csv(filesIwant[i])
#table(open$Phenotype, open$section) #different section has different number of CD4 and CD8 detected

#mean number of detected CD4+ between 2 sections
mtx[i,2] <- round((length(which(open$Phenotype == "CD4+" & open$section == "Top")) + length(which(open$Phenotype == "CD4+" & open$section == "Bottom")))/2,0)
mtx[i,3] <- round((length(which(open$Phenotype == "CD8+" & open$section == "Top")) + length(which(open$Phenotype == "CD8+" & open$section == "Bottom")))/2,0)

mtx[i,4] <- (mean(open$Cytoplasm.iFN..Opal.520..Mean..Normalized.Counts..Total.Weighting.[which(open$Phenotype == "CD4+" & open$section == "Top")]) + mean(open$Cytoplasm.iFN..Opal.520..Mean..Normalized.Counts..Total.Weighting.[which(open$Phenotype == "CD4+" & open$section == "Bottom")]) /2 )
mtx[i,5] <- (mean(open$Cytoplasm.iFN..Opal.520..Mean..Normalized.Counts..Total.Weighting.[which(open$Phenotype == "CD8+" & open$section == "Top")]) + mean(open$Cytoplasm.iFN..Opal.520..Mean..Normalized.Counts..Total.Weighting.[which(open$Phenotype == "CD8+" & open$section == "Bottom")]) /2 )

mtx[i,6] <- (mean(open$Cytoplasm.GzB..Opal.570..Mean..Normalized.Counts..Total.Weighting.[which(open$Phenotype == "CD4+" & open$section == "Top")]) + mean(open$Cytoplasm.GzB..Opal.570..Mean..Normalized.Counts..Total.Weighting.[which(open$Phenotype == "CD4+" & open$section == "Bottom")]) /2 )
mtx[i,7] <- (mean(open$Cytoplasm.GzB..Opal.570..Mean..Normalized.Counts..Total.Weighting.[which(open$Phenotype == "CD8+" & open$section == "Top")]) + mean(open$Cytoplasm.GzB..Opal.570..Mean..Normalized.Counts..Total.Weighting.[which(open$Phenotype == "CD8+" & open$section == "Bottom")]) /2 )
}


for (i in 13:13){
  open <- read.csv(filesIwant[i])
  #table(open$Phenotype, open$section) #different section has different number of CD4 and CD8 detected
  
  #mean number of detected CD4+ between 2 sections
  mtx[i,2] <- length(which(open$Phenotype == "CD4+" & open$section == "Top"))
  mtx[i,3] <- length(which(open$Phenotype == "CD8+" & open$section == "Top"))
  
  mtx[i,4] <- (mean(open$Cytoplasm.iFN..Opal.520..Mean..Normalized.Counts..Total.Weighting.[which(open$Phenotype == "CD4+" & open$section == "Top")]))
  mtx[i,5] <- (mean(open$Cytoplasm.iFN..Opal.520..Mean..Normalized.Counts..Total.Weighting.[which(open$Phenotype == "CD8+" & open$section == "Top")]))
  
  mtx[i,6] <- (mean(open$Cytoplasm.GzB..Opal.570..Mean..Normalized.Counts..Total.Weighting.[which(open$Phenotype == "CD4+" & open$section == "Top")]))
  mtx[i,7] <- (mean(open$Cytoplasm.GzB..Opal.570..Mean..Normalized.Counts..Total.Weighting.[which(open$Phenotype == "CD8+" & open$section == "Top")]))
}

write.csv(mtx, "Summary_CD4CD8.csv")