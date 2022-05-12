###################################################
###  GBLUP/population differentiation	analysis  ###
###                   Ben Harms                 ###
###################################################
### randomization for prediction sets - modified from Arthur Bernardeli's code


require(AGHmatrix)
require(asreml)


###############################################
#######  Phenotype data prep for GBLUP  #######
###############################################

# Phenotype file: 3 columns, yield stability response in col 3
# Phenotypic data for training population(YS values), and VP - NA's for phenotype
#TPpheno <- read.csv("/Users/BenHarms/Documents/GBLUP/final932/tp_YSpheno.csv") #TP phenotypic data yield stability
#VPpheno <- read.csv("/Users/BenHarms/Documents/GBLUP/final932/vp_NA.csv")  #VP phenotype data - NA
TPpheno <- read.csv("data/tp2_YSpheno.csv")
VPpheno <- read.csv("data/vp2_NA.csv")


# set up TP and VP phenotype files
colnames(TPpheno) <- c("ENTRY_ID","Genotype","YieldStability")
colnames(VPpheno) <- c("ENTRY_ID","Genotype","YieldStability")
head(TPpheno)
dim(TPpheno)
head(VPpheno)
dim(VPpheno)

TPpheno$ENTRY_ID<-factor(TPpheno$ENTRY_ID)
VPpheno$ENTRY_ID<-factor(VPpheno$ENTRY_ID)
TPpheno$Genotype<-factor(TPpheno$Genotype)
VPpheno$Genotype<-factor(VPpheno$Genotype)
# set up G-matrix later when TP individuals are assigned





################################################
####### prepare GS sets w/ randomized TP #######
################################################

# randomly assign fold 1:2 for 10 reps - each TP line
# TP lines in fold 1 will be used for prediction/Fst calculation 92 TP - 118 VP
set.seed(1)
f <- data.frame(f01=sample(rep(1:2,92)[1:185], 185, replace=F),f02=sample(rep(1:2,92)[1:185], 185, replace=F),
                f03=sample(rep(1:2,92)[1:185], 185, replace=F),f04=sample(rep(1:2,92)[1:185], 185, replace=F),
                f05=sample(rep(1:2,92)[1:185], 185, replace=F),f06=sample(rep(1:2,92)[1:185], 185, replace=F),
                f07=sample(rep(1:2,92)[1:185], 185, replace=F),f08=sample(rep(1:2,92)[1:185], 185, replace=F),
                f09=sample(rep(1:2,92)[1:185], 185, replace=F),f10=sample(rep(1:2,92)[1:185], 185, replace=F))


INF <- data.frame(H=levels(factor(TPpheno$Genotype)),f)
for(i in 1:nrow(INF)){
  TPpheno$f01[TPpheno$Genotype==paste0(INF[i,1])] <- INF[i,2] 
  TPpheno$f02[TPpheno$Genotype==paste0(INF[i,1])] <- INF[i,3] 
  TPpheno$f03[TPpheno$Genotype==paste0(INF[i,1])] <- INF[i,4] 
  TPpheno$f04[TPpheno$Genotype==paste0(INF[i,1])] <- INF[i,5] 
  TPpheno$f05[TPpheno$Genotype==paste0(INF[i,1])] <- INF[i,6] 
  TPpheno$f06[TPpheno$Genotype==paste0(INF[i,1])] <- INF[i,7]  
  TPpheno$f07[TPpheno$Genotype==paste0(INF[i,1])] <- INF[i,8]  
  TPpheno$f08[TPpheno$Genotype==paste0(INF[i,1])] <- INF[i,9]  
  TPpheno$f09[TPpheno$Genotype==paste0(INF[i,1])] <- INF[i,10]  
  TPpheno$f10[TPpheno$Genotype==paste0(INF[i,1])] <- INF[i,11]  
}

head(TPpheno)

# take fold 1 from each rep (92 lines)
tp1 <- TPpheno[TPpheno$f01 == 1,]
tp1 <- tp1[,1:3]
tp2 <- TPpheno[TPpheno$f02 == 1,]
tp2 <- tp2[,1:3]
tp3 <- TPpheno[TPpheno$f03 == 1,]
tp3 <- tp3[,1:3]
tp4 <- TPpheno[TPpheno$f04 == 1,]
tp4 <- tp4[,1:3]
tp5 <- TPpheno[TPpheno$f05 == 1,]
tp5 <- tp5[,1:3]
tp6 <- TPpheno[TPpheno$f06 == 1,]
tp6 <- tp6[,1:3]
tp7 <- TPpheno[TPpheno$f07 == 1,]
tp7 <- tp7[,1:3]
tp8 <- TPpheno[TPpheno$f08 == 1,]
tp8 <- tp8[,1:3]
tp9 <- TPpheno[TPpheno$f09 == 1,]
tp9 <- tp9[,1:3]
tp10 <- TPpheno[TPpheno$f10 == 1,]
tp10 <- tp10[,1:3]

# compile phenotype file for each GS rep with random TP, same VP w/ NA
gs1 <- rbind(tp1,VPpheno)
gs2 <- rbind(tp2,VPpheno)
gs3 <- rbind(tp3,VPpheno)
gs4 <- rbind(tp4,VPpheno)
gs5 <- rbind(tp5,VPpheno)
gs6 <- rbind(tp6,VPpheno)
gs7 <- rbind(tp7,VPpheno)
gs8 <- rbind(tp8,VPpheno)
gs9 <- rbind(tp9,VPpheno)
gs10 <- rbind(tp10,VPpheno)

#drop any rows with all NA
gs1 <- gs1[!is.na(gs1$Genotype),]
gs2 <- gs2[!is.na(gs2$Genotype),]
gs3 <- gs3[!is.na(gs3$Genotype),]
gs4 <- gs4[!is.na(gs4$Genotype),]
gs5 <- gs5[!is.na(gs5$Genotype),]
gs6 <- gs6[!is.na(gs6$Genotype),]
gs7 <- gs7[!is.na(gs7$Genotype),]
gs8 <- gs8[!is.na(gs8$Genotype),]
gs9 <- gs9[!is.na(gs9$Genotype),]
gs10 <- gs10[!is.na(gs10$Genotype),]

#set row names as GenotypeID column - read in SNP data to datag
rownames(gs1) <- gs1[,2]
rownames(gs2) <- gs2[,2]
rownames(gs3) <- gs3[,2]
rownames(gs4) <- gs4[,2]
rownames(gs5) <- gs5[,2]
rownames(gs6) <- gs6[,2]
rownames(gs7) <- gs7[,2]
rownames(gs8) <- gs8[,2]
rownames(gs9) <- gs9[,2]
rownames(gs10) <- gs10[,2]





#############################################
######## genomic relationship matrix ########
#############################################
# calculate separate matrix for each GS run - probably not necessary

library(AGHmatrix)

#read SNP file
datag <- read.csv(file = "data/geno.csv")
datag <- as.data.frame(datag)

#set row names as GenotypeID column - read in SNP data to datag
rownames(datag) <- datag[,1]

#rename first column to genotype in snp matrix
names(datag)[names(datag) == 'X'] <- 'Genotype'

dropdatap <- subset(gs1,!(Genotype%in%datag$Genotype))
dropdatapp <- as.character(dropdatap[,1])
gs3 <- gs3[!(row.names(gs1) %in% dropdatapp),]

#do the same with genotype matrix (in case genotypes in datag are not present in datap)
dropdatag <- subset(datag,!(Genotype%in%gs9$Genotype))
dropdatagg <- as.character(dropdatag[,1])
datag <- datag[!(row.names(datag) %in% dropdatagg),]

#format matrix for gmatrix calculation
rnames <- datag[,1]
mat_data <- data.matrix(datag[,2:ncol(datag)])
rownames(mat_data) <- rnames
rownames(datag) <- rnames

# check SNP matrix
mat_data[1:5,1:5]
str(mat_data)

# calculate G-matrix using AGHmatrix package
gmatrix <- Gmatrix(SNPmatrix=mat_data, maf=0.05, missingValue=-9,
                   method="VanRaden", ploidy=2)
det(gmatrix) #if 0, matrix has no inverse

# prep G-matrix for ASReml
GinvVanRaden3 <- solve(gmatrix)
Ginv3 <- formatmatrix(GinvVanRaden3, return = TRUE)

colnames(Ginv3) <- c("Row", "Column", "Value")
head(Ginv3)
attr(Ginv3, "rowNames") <- as.character(gs3$ENTRY_ID)
attr(Ginv3, "colNames") <- as.character(gs3$ENTRY_ID)
attr(Ginv3, "INVERSE") <- TRUE





####################################
#############  GBLUP   #############
####################################

# run modelfor each GS pheno set (gs_) and matrix (Ginv_)
modelGBLUP1 <- asreml(fixed=YieldStability ~1, random=~vm(ENTRY_ID,Ginv1), 
                     workspace=7e08, na.action=na.method(x="include"), data=gs1)

# plot residuals
  pdf("modelgblupr1.pdf")
  plot(modelGBLUP1)
  dev.off()

# summary of variance components
summary(modelGBLUP1)$varcomp

# get BLUP values for prediction accuracy, TP accuracy
BLUP<-summary(modelGBLUP1,coef=TRUE)$coef.random
head(BLUP)

# TP -> TP accuracy
cor(gs1$YieldStability,BLUP,method="pearson",use="complete.obs")
tp1 <- 0.65
tp2 <- 0.86
tp3 <- 0.47
tp4 <- 0.87
tp5 <- 0.8881035
tp6 <- 0.6646904
tp7 <- 0.731858
tp8 <- 0.6713137
tp9 <- 0.9055388
tp10 <- 0.7571073

# prediction accuracy for validation population (pheno = NA)
vpys <- read.csv("data/vp2_YSpheno.csv") # make from phenovalid_total.csv file
VPBLUP <- BLUP
rnames <- gs1[,1]
rownames(VPBLUP) <- rnames
VPBLUP <- VPBLUP[93:nrow(VPBLUP),]
#bind observed and predicted to same dataframe for correlation/accuracy
df1 <- vpys[,2:3] #strainID and YieldStability from observed VP
blupVP <- VPBLUP[,1] #predicted BLUP from GS model
merged <- cbind(df1,blupVP)
r1 <- cor(merged$YieldStability,merged$blupVP,method="pearson",use="complete.obs")   


# final results table
# TP -> VP accuracy
VP_r <- c(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10)
# TP -> TP accuracy
TP_r <- c(tp1,tp2,tp3,tp4,tp5,tp6,tp7,tp8,tp9,tp10)
rep <- c(1:10)
results <- data.frame(rep,TP_r,VP_r)
write.csv(results, "pred_results.csv")




############################################
############# Fst calculation ##############
############################################

##### read in SNP data, run separately with prediction iterations (gs1,gs2,gs3...)
#geno1 <- read.csv("/Users/BenHarms/Documents/Classes/AGRO 932 - Biometrical breeding and genetics/midterm/data/TPandVP_MIPs.csv")
geno1 <- read.csv("data/TPandVP_MIPs.csv")

geno1 <- geno1[,5:ncol(geno1)]

############ start here for each Fst rep ############
# rearrange SNP data to line/genotype-row, marker-column
genotr <- t(geno1)
genotr <- as.data.frame(genotr)
colnames(genotr) <- genotr[1,]

# formatting SNP data
genotr[,1] <- rownames(genotr)
colnames(genotr) <- genotr[1,]
names(genotr)[names(genotr) == 'ID'] <- 'Genotype'
genotr$Genotype <- gsub("\\.", "-", genotr$Genotype)
rownames(genotr) <- genotr[,1]

# remove all lines not in specific prediction TP/VP set
# run this with gs1,gs2,gs3...gs10 for correct combination of TP/VP lines
dropdatag <- subset(genotr,!(Genotype%in%gs10$Genotype))
dropdatagg <- as.character(dropdatag[,1])
genotr <- genotr[!(row.names(genotr) %in% dropdatagg),]

# transpose back to marker-row, line/genotype-column
geno <- t(genotr)
geno <- as.data.frame(geno)
geno <- geno[-1,]

##### Format VP and TP in one dataset, split alleles into separate columns
for(i in 1:ncol(geno)){
  # replace slash and everything after it as nothing
  geno$newcol <- gsub("/.*", "", geno[,i] )
  # extract the line name
  nm <- names(geno)[i]
  # assign name for this allele
  names(geno)[ncol(geno)] <- paste0(nm, sep="_a1")
  geno$newcol <- gsub(".*/", "", geno[,i] )
  names(geno)[ncol(geno)] <- paste0(nm, sep="_a2")
}
dim(geno)

############### calculate Fst for 92 TP lines - 118 VP lines
# replace . with NA in allele columns
for(i in 211:630){
  geno[,i] <- gsub("\\.","NA",geno[,i]) 
  # need \\ before . to identify as special character
}

######## Count alleles for TP-1 and VP-2
## for two populations
geno$pop10 <- rowSums(geno[, 211:394] == "0")
geno$pop11 <- rowSums(geno[, 211:394] == "1")
geno$pop20 <- rowSums(geno[, 395:630] == "0")
geno$pop21 <- rowSums(geno[, 395:630] == "1")
geno$npop  <- rowSums(geno[, 631:634])
geno$npop1 <- rowSums(geno[, 631:632])
geno$npop2 <- rowSums(geno[, 633:634])

# Compute p, p1, p2 from allele columns(20 for 20 individuals) add 1's together and divide by # of individuals in population
# pop1 is l1 - l10
# pop2 is l11 - l20
geno$p <- apply(geno[, 211:630], 1, function(x) {sum(as.numeric(as.character(x)),na.rm = TRUE)})
geno$p <- geno$p/geno$npop
geno$p1 <- apply(geno[, 211:394], 1, function(x) {sum(as.numeric(as.character(x)),na.rm = TRUE)})
geno$p1 <- geno$p1/geno$npop1
geno$p2 <- apply(geno[, 395:630], 1, function(x) {sum(as.numeric(as.character(x)), na.rm = TRUE)})
geno$p2 <- geno$p2/geno$npop2


# Calculate Fst using p, p1, p2
geno$fst <- with(geno, ((p1-p)^2 + (p2-p)^2)/(2*p*(1-p)))

# remove NaN values in geno$Fst
geno <- subset(geno, fst != "NaN")

#calculate average fst: sum Fst/#loci  -- change object name each separate calculation
totalfst1 <- sum(geno$fst, na.rm = TRUE)
fst3 <- totalfst1/nrow(geno)
################### stop here for each Fst rep


#compile fst values
fst <- c(fst1,fst2,fst3,fst4,fst5,fst6,fst7,fst8,fst9,fst10)
rep <- c("r1","r2","r3","r4","r5","r6","r7","r8","r9","r10")
fst_results <- data.frame(rep,fst)
write.csv(fst_results, "fst_results.csv")



#################################################
############# theta pi calculation ##############
#################################################

######## calculate nuc. diversity for both populations at each position then average across all loci
### theta for VP and TP across same loci
pi <- function(n=20, p=geno){
  return(n/(n-1)*(1-p^2-(1-p)^2))
}
## Create new column for population 1 pi at each position
geno$pi1 <- pi(n=20, p=geno[, 639])
## Create new column for population 2 pi at each position
geno$pi2 <- pi(n=20, p=geno[, 640])

### Population 1 theta: add pi values at all loci and divide by number of loci
totalpi1 <- sum(geno$pi1, na.rm = TRUE) #for some reason there are missing values
loci <- 809 #number of loci
thetatp10 <- totalpi1/loci
# output: 0.3498

### Population 2 theta: add pi values at all loci and divide by number of loci
totalpi2 <- sum(geno$pi2, na.rm = TRUE)
loci <- 809
thetavp5 <- totalpi2/loci
# output: 0.302

pi <- c(thetatp1,thetatp2,thetatp3,thetatp4,thetatp5,thetatp6,thetatp7,thetatp8,thetatp9,thetatp10)
rep <- c(1:10)
pi_results <- data.frame(rep,pi)
write.csv(pi_results, "pi_results.csv")



###########################################
#############  Plot results   #############
###########################################

library(ggplot2)
library(ggpubr)

# visualize Fst results
# fst <- read.csv("cache/fst.csv")
pdf("TPandVP_Fst.pdf")
plot(geno$POS, geno$fst, xlab="Physical position", ylab="Fst value", main="TP-VP Fst")
dev.off()

############ visualize genomic prediction results accuracy vs. Fst
fst <- read.csv("/Users/BenHarms/Documents/GBLUP/final932/gblupfinal2/fst_results.csv")
predacc <- read.csv("/Users/BenHarms/Documents/GBLUP/final932/gblupfinal2/pred_results.csv")
rfst <- data.frame(fst[,3],predacc)
rfst <- rfst[,-2]
colnames(rfst) <- c("fst","rep","TP_r","VP_r")
#ggplot scatter with regression line VP prediction accuracy
t <- ggplot(rfst, aes(x=fst, y=VP_r)) + 
      geom_point() +
      geom_smooth(method=lm, se=FALSE) +
      stat_regline_equation(label.y = -0.08, aes(label = ..eq.label..)) +
      stat_regline_equation(label.y = -0.1, aes(label = ..rr.label..))
t + xlab("population differentiation (Fst)") + ylab("TP-VP prediction accuracy")
#ggplot scatter with regression line TP prediction accuracy
t1 <- ggplot(rfst, aes(x=fst, y=TP_r)) + 
        geom_point()+
        geom_smooth(method=lm, se=FALSE)+
        stat_regline_equation(label.y = 0.52, aes(label = ..eq.label..)) +
        stat_regline_equation(label.y = 0.48, aes(label = ..rr.label..))
t1 + xlab("population differentiation (Fst)") + ylab("TP-TP prediction accuracy")


############ visualize genomic prediction results accuracy vs. theta pi
pi <- read.csv("/Users/BenHarms/Documents/GBLUP/final932/gblupfinal2/pi_results.csv")
predacc <- read.csv("/Users/BenHarms/Documents/GBLUP/final932/gblupfinal2/pred_results.csv")
rpi <- data.frame(pi[,3],predacc)
rpi <- rfst[,-2]
colnames(rfst) <- c("pi","rep","TP_r","VP_r")
#ggplot scatter with regression line VP prediction accuracy
p <- ggplot(rpi, aes(x=pi, y=VP_r)) + 
      geom_point() +
      geom_smooth(method=lm, se=FALSE) +
      stat_regline_equation(label.y = -0.08, aes(label = ..eq.label..)) +
      stat_regline_equation(label.y = -0.1, aes(label = ..rr.label..))
p + xlab("TP diversity") + ylab("TP-VP prediction accuracy")
#ggplot scatter with regression line TP prediction accuracy
p2 <- ggplot(rpi, aes(x=pi, y=TP_r)) + 
        geom_point()+
        geom_smooth(method=lm, se=FALSE)+
        stat_regline_equation(label.y = 0.52, aes(label = ..eq.label..)) +
        stat_regline_equation(label.y = 0.48, aes(label = ..rr.label..))
p2 + xlab("TP diversity") + ylab("TP-TP prediction accuracy")



####  Yield stability phenotype distributions
ys <- read.csv(file = "/Users/BenHarms/Documents/GBLUP/final932/figures/yspheno_tpvp.csv")
# violin plot
p <- ggplot(ys, aes(x=pop_ID , y=YieldStability)) + 
  geom_violin()
p


