library(tidyverse)
install.packages("writexl")
install.packages("arsenal")
library(arsenal)
library("writexl")

current1Data <- read.table(file = "clipboard", sep = "\t", header=TRUE)
View(current1Data)

raw1Data <- read.table(file = "clipboard", sep = "\t")
View(raw1Data)
temp = raw1Data
temp[temp == ""] <- NA 
raw1Data = temp
raw1Data['new_col'] <- NA


tempDF = data.frame(matrix(nrow = 0, ncol = ncol(raw1Data)))
moms = data.frame(matrix(nrow = 0, ncol = ncol(raw1Data)))
for (i in seq(1, nrow(raw1Data), by = 2)) {
  mother = raw1Data[c(i),]
  #child = raw1Data[c(i+1), c(4:ncol(raw1Data))]
  moms = rbind(moms, mother)
}

View(moms)

kids = data.frame(matrix(nrow = 0, ncol = ncol(raw1Data)))
for (i in seq(2, nrow(raw1Data), by = 2)) {
  #mother = raw1Data[c(i),]
  child = raw1Data[c(i),]
  kids = rbind(kids, child)
}

View(kids)
warnings()
m = moms[-c(1:3)]
k = kids[-c(1:3)]

identical(m,k)  

summary(comparedf(m,k))


`for (i in seq(2, nrow(raw1Data), by = 2)) {
  tempDF = setdiff(moms[i],child[i])
}


datacols = c("ID", "smt","Sample","Variant")
convertedDF = data.frame(matrix(nrow = 0, ncol = 4))
rowsRawData = nrow(raw1Data)
for(i in 1:rowsRawData){
  
  colNum = 4
  while(is.na(raw1Data[i,colNum]) == FALSE){
    newRow = c(raw1Data[i,1],raw1Data[i,2],raw1Data[i,3],raw1Data[i,colNum])
    convertedDF[nrow(convertedDF) + 1,] <- newRow
    colNum = colNum +1
  }
}



View(convertedDF)
write_xlsx(convertedDF,"study1data.xlsx")

noM = current1Data$X0116.009M
noM = substr(noM,1,8)
noM

noM = current1Data$X0116.010C
noM = substr(noM,1,8)
noM

current1Data['X0116.010C'] <- noM

noM = convertedDF$X4
noM<-as.numeric(gsub("([0-9]+).*$", "\\1", noM))
noM
convertedDF['pos'] <- noM
testDF = convertedDF
View(testDF)
testDF <- testDF[-c(1), ]




num = 1
for (i in 1:nrow(convertedDF)) {
  if(convertedDF$X1[i] %in% current1Data$X0116.009M){
    if(convertedDF$pos[i] %in% current1Data$X14000){
      testDF <- testDF[-c(i), ]
      print(num)
      num = num +1
    }
  }
}

dim(testDF)
dim(convertedDF)

for (i in 1:nrow(convertedDF)) {
  if(convertedDF$X1[i] %in% current1Data$X0116.010C){
    if(convertedDF$pos[i] %in% current1Data$X14000){
      testDF <- testDF[-c(i), ]
      print(num)
      num = num +1
    }
  }
}
dim(testDF)
dim(convertedDF)

write_xlsx(testDF,"study1data.xlsx")


my_data <- read.table(file = "clipboard", sep = "\t", header=TRUE)
View(my_data)




(0.0002093501 /20)*1000000


10.4675/2.6




len = length(my_data$Proband.ID)
len


momFreq = c()

for(i in 1426:1439){
  proband = my_data[i,8]
  mother = my_data[i,7]
  pos = my_data[i,3]
  sample = my_data[i,11]
  print(filter(my_data, Proband.ID == mother & Position.rCRS. == pos & Sample == sample)[[6]])
  
  
  
  if(any(my_data$Proband.ID==mother)){
    momFreq = c(momFreq,filter(my_data, Proband.ID == mother & Position.rCRS. == pos & Sample == sample)[[6]])
  }
  else{
    momFreq = c(momFreq,NA)
  }
  print(i)  
}

df = data.frame(momFreq)
write_xlsx(df,"newDATA.xlsx")


any(my_data$Proband.ID=="F113g1")


which(my_data[8]=="F113g1", arr.ind=TRUE)

print(my_data[8],max = 2000)
momFreq
View(momFreq)
length((momFreq))

newData = my_data$mom.freq <-momFreq
View(my_data)
filter(my_data, Proband.ID == "11002m" & Position.rCRS. == "9029")[[6]]



#####################################################################################################################


library("readxl")
vcf <- read_excel("New Format Data 02.07.xlsx",sheet = 2)
View(vcf)

newvcf = head(vcf,-1)
View(newvcf)

NumOfMuts = nrow(newvcf)
NumOfMuts

bloodVCF <- newvcf %>% filter(Sample == "bl")
View(bloodVCF)


cheekVCF <- newvcf %>% filter(Sample == "ch")
View(cheekVCF)

probandsBl = unique(bloodVCF$`Proband ID`)
probandsCh = unique(cheekVCF$`Proband ID`)
probandsTot = unique(newvcf$`Proband ID`)

length(probandsCh)

mtlength = 16569

bloodMutRate = (nrow(bloodVCF)/length(probandsBl))/mtlength
cheekMutRate = (nrow(cheekVCF)/length(probandsCh))/mtlength
totMutRate = (nrow(newvcf)/length(probandsTot))/mtlength
bloodMutRate
cheekMutRate
totMutRate


b30 = (bloodMutRate/30)*1000000
c30 = (cheekMutRate/30)*1000000
t30 = (totMutRate/30)*1000000


lowend = 0.43
highend = 2.6

b30/lowend
b30/highend

c30/lowend
c30/highend

t30/lowend
t30/highend

plot(bloodVCF$`Family ID`, bloodVCF$`Mom Freq`)


study1VCF <- newvcf %>% filter(`Study ID` == 1)
probandsS1 = unique(study1VCF$`Proband ID`)
length(probandsS1)
s1MutRate = (nrow(study1VCF)/length(probandsS1))/mtlength
s1MutRate

study2VCF <- newvcf %>% filter(`Study ID` == 2)
probandsS2 = unique(study2VCF$`Proband ID`)
length(probandsS2)
s2MutRate = (nrow(study2VCF)/length(probandsS2))/mtlength
s2MutRate

study3VCF <- newvcf %>% filter(`Study ID` == 3)
probandsS3 = unique(study3VCF$`Proband ID`)
length(probandsS3)
s3MutRate = (nrow(study3VCF)/length(probandsS3))/mtlength
s3MutRate
View(study3VCF)

study4VCF <- newvcf %>% filter(`Study ID` == 4)
probandsS4 = unique(study4VCF$`Proband ID`)
length(probandsS4)
s4MutRate = (nrow(study4VCF)/length(probandsS4))/mtlength
s4MutRate
studyRates = c(s1MutRate,s2MutRate,s3MutRate,s4MutRate)
StudyNums = c(1,2,3,4)
df = data.frame(StudyNums,studyRates)
View(df)
plot(df)




#filtering 

filteredBloodVCF <- bloodVCF %>% filter(`Minor Allele Freq?` <= 0.002)
filterNum = length(bloodVCF) - length(filteredBloodVCF)
#View(filteredBloodVCF)
filteredProbandsBl = unique(filteredBloodVCF$`Proband ID`)
filteredBloodMutRate = (nrow(filteredBloodVCF)/length(filteredProbandsBl))/(mtlength-filterNum)
filteredBloodMutRate

filteredCheekVCF <- cheekVCF %>% filter(`Minor Allele Freq?` <= 0.001)
filterNum = length(cheekVCF) - length(filteredCheekVCF)
filteredProbandsCh = unique(filteredCheekVCF$`Proband ID`)
filteredCheekMutRate = (nrow(filteredCheekVCF)/length(filteredProbandsCh))/(mtlength-filterNum)
filteredCheekMutRate





#plotting
#blood
datacols = c("Filter", "Rate Estimates (mutations/site/generation)")
filterDFBl = data.frame(matrix(nrow = 0, ncol = 2))
colnames(filterDFBl) = datacols
View(filterDFBl)
for (i in seq(0, 0.6, by = 0.0001)) {
  filter = i
  filteredBloodVCF <- bloodVCF %>% filter(`Minor Allele Freq?` <= filter)
  #filteredBloodVCF <- bloodVCF %>% filter(`Mom Freq` <= filter)
  filterNum = length(bloodVCF) - length(filteredBloodVCF)
  #View(filteredBloodVCF)
  filteredProbandsBl = unique(filteredBloodVCF$`Proband ID`)
  filteredBloodMutRate = (nrow(filteredBloodVCF)/length(filteredProbandsBl))/(mtlength-filterNum)
  #filteredBloodMutRate
  newFilterData<-data.frame(filter, filteredBloodMutRate)
  names(newFilterData)<-c("Filter", "Rate Estimates")
  
  filterDFBl <- rbind(filterDFBl, newFilterData)
  
}

View(filterDFBl)
par(mfrow=c(1,2))
plot(filterDFBl, ylim = c(0,0.00014),xaxt='n', main='A', ylab = "Rate Estimates (mutations/site/generation)")
axis(1,at=seq(0,0.6,by=0.05))

#cheek
filterDFCh = data.frame(matrix(nrow = 0, ncol = 2))
colnames(filterDFCh) = datacols
View(filterDFCh)
for (i in seq(0, 0.6, by = 0.0001)) {
  filter = i
  filteredCheekVCF <- cheekVCF %>% filter(`Minor Allele Freq?` <= filter)
  #filteredCheekVCF <- cheekVCF %>% filter(`Mom Freq` <= filter)
  filterNum = length(cheekVCF) - length(filteredCheekVCF)
  filteredProbandsCh = unique(filteredCheekVCF$`Proband ID`)
  filteredCheekMutRate = (nrow(filteredCheekVCF)/length(filteredProbandsCh))/(mtlength-filterNum)
  newFilterData<-data.frame(filter, filteredCheekMutRate)
  names(newFilterData)<-c("Filter", "Rate Estimates")
  
  filterDFCh <- rbind(filterDFCh, newFilterData)
  
}

View(filterDFCh)
plot(filterDFCh, ylim = c(0,0.00025),xaxt='n',main="B", ylab = "Rate Estimates (mutations/site/generation)")
axis(1,at=seq(0,0.6,by=0.05))

#MOM
datacols = c("Filter", "Rate Estimates (mutations/site/generation)")
filterDFBl = data.frame(matrix(nrow = 0, ncol = 2))
colnames(filterDFBl) = datacols
View(filterDFBl)
for (i in seq(0, 0.6, by = 0.0001)) {
  filter = i
  filteredBloodVCF <- bloodVCF %>% filter(`Mom Freq` <= filter)
  filterNum = length(bloodVCF) - length(filteredBloodVCF)
  #View(filteredBloodVCF)
  filteredProbandsBl = unique(filteredBloodVCF$`Proband ID`)
  filteredBloodMutRate = (nrow(filteredBloodVCF)/length(filteredProbandsBl))/(mtlength-filterNum)
  #filteredBloodMutRate
  newFilterData<-data.frame(filter, filteredBloodMutRate)
  names(newFilterData)<-c("Filter", "Rate Estimates")
  
  filterDFBl <- rbind(filterDFBl, newFilterData)
  
}

View(filterDFBl)
par(mfrow=c(1,2))
plot(filterDFBl, ylim = c(0,0.00014),xaxt='n', main="A",ylab = "Rate Estimates (mutations/site/generation)")
axis(1,at=seq(0,0.6,by=0.05))

#cheek
filterDFCh = data.frame(matrix(nrow = 0, ncol = 2))
colnames(filterDFCh) = datacols
View(filterDFCh)
for (i in seq(0, 0.6, by = 0.0001)) {
  filter = i
  filteredCheekVCF <- cheekVCF %>% filter(`Mom Freq` <= filter)
  filterNum = length(cheekVCF) - length(filteredCheekVCF)
  filteredProbandsCh = unique(filteredCheekVCF$`Proband ID`)
  filteredCheekMutRate = (nrow(filteredCheekVCF)/length(filteredProbandsCh))/(mtlength-filterNum)
  newFilterData<-data.frame(filter, filteredCheekMutRate)
  names(newFilterData)<-c("Filter", "Rate Estimates")
  
  filterDFCh <- rbind(filterDFCh, newFilterData)
  
}

View(filterDFCh)
plot(filterDFCh, ylim = c(0,0.00025),xaxt='n',main="B",ylab = "Rate Estimates (mutations/site/generation)")
axis(1,at=seq(0,0.6,by=0.05))














layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
hist(newvcf$`Position(rCRS)`,xlim = c(0,18000),ylim = c(0,300),xaxt='n', yaxt='n', breaks=100, 
     main = "A",
     xlab = "Variant Position")
axis(1,at=seq(0,18000,by=2000))
axis(2,at=seq(0,300,by=50))

hist(bloodVCF$`Position(rCRS)`,xlim = c(0,18000),ylim = c(0,200),xaxt='n', yaxt='n', breaks=100, 
     main = "B",
     xlab = "Variant Position")
axis(1,at=seq(0,18000,by=2000))
axis(2,at=seq(0,300,by=50))
hist(cheekVCF$`Position(rCRS)`,xlim = c(0,18000),ylim = c(0,100),xaxt='n', yaxt='n', breaks=100, 
     main = "C",
     xlab = "Variant Position")
axis(1,at=seq(0,18000,by=2000))
axis(2,at=seq(0,100,by=10))
#MAF


h1 = hist(cheekVCF$`Minor Allele Freq?`,xlim = c(0,0.4),ylim = c(0,1000),xaxt='n', yaxt='n', breaks=95, 
     main = "B",
     xlab = "Proband Minor Allele Frequencies")

h1$density=h1$counts/sum(h1$counts)
plot(h1,freq=FALSE,xlim = c(0,0.4),xaxt='n',
     main = "B",
     xlab = "Child Minor Allele Frequencies",ylab = "Normalized Frequency")


axis(1,at=seq(0,1.0,by=0.01))
axis(2,at=seq(0,1000,by=250))



h2 = hist(bloodVCF$`Minor Allele Freq?`,xlim = c(0,0.4),ylim = c(0,1000),xaxt='n', yaxt='n', breaks=100, 
          main = "A",
          xlab = "Proband Minor Allele Frequencies")

h2$density=h2$counts/sum(h2$counts)
plot(h2,freq=FALSE,xlim = c(0,0.4),xaxt='n',yaxt = 'n',ylim = c(0,0.7),
     main = "A",
     xlab = "Child Minor Allele Frequencies",ylab = "Normalized Frequency")


axis(1,at=seq(0,1.0,by=0.01))
axis(2,at=seq(0,0.7,by=0.1))


axis(1,at=seq(0,1.0,by=0.01))
axis(2,at=seq(0,0.7,by=0.1))
options(scipen = 0)

install.packages("plotly")
library(plotly)


marker_style <- list(line = list(width = 1,color = 'rgb(0, 10, 0)'));

p3 <- plot_ly(x = bloodVCF$`Minor Allele Freq?`,
              type = "histogram",
              histnorm = "probability",
              marker = list(color = 'white',
                            line = list(color = ("black"),width = 1.0))) %>% 
  layout(title = "<b>A<b>",
         xaxis = list(title = "<b>Child Minor Allele Frequency<b>",
                      zeroline = FALSE,
                      ticks="outside",
                      nticks=20,
                      dtick = 0.03),
         yaxis = list(title = "<b>Normalized Frequency<b>",
                      zeroline = TRUE,
                      ticks="outside",
                      nticks=20), bargap = 0.1)
p3

p4 <- plot_ly(x = cheekVCF$`Minor Allele Freq?`,
              type = "histogram",
              histnorm = "probability",
              marker = list(color = 'white',
              line = list(color = ("black"),width = 1.0))) %>% 
  layout(title = "<b>B<b>",
         xaxis = list(title = "<b>Child Minor Allele Frequency<b>",
                      zeroline = FALSE,
                      ticks="outside",
                      nticks=20,dtick = 0.03),
         yaxis = list(title = "<b>Normalized Frequency<b>",
                      zeroline = TRUE,
                      ticks="outside",
                      nticks=20), bargap = 0.1)
p4

subplot(p3, p4,shareX = TRUE,shareY =TRUE)%>%layout(
  annotations = list( 
    list(x = -0.01 , y = 1.04, font = list(size = 16),text = "<b>A<b>",showarrow = F,xref='paper', yref='paper')),
                margin = list(b = 40, l = 100))
                




p5 <- plot_ly(x = bloodVCF$`Mom Freq`,
              type = "histogram",
              histnorm = "probability",
              marker = list(color = 'white',
                            line = list(color = ("black"),width = 1.0))) %>% 
  layout(title = "<b>A<b>",
         xaxis = list(title = "<b>Mother Minor Allele Frequency<b>",
                      zeroline = FALSE,
                      ticks="outside"),
         yaxis = list(title = "<b>Normalized Frequency<b>",
                      zeroline = TRUE,
                      ticks="outside"))
p5

p6 <- plot_ly(x = cheekVCF$`Mom Freq`,
              type = "histogram",
              histnorm = "probability",
              marker = list(color = 'white',
                            line = list(color = ("black"),width = 1.0))) %>% 
  layout(title = "<b>B<b>",
         xaxis = list(title = "<b>Mother Minor Allele Frequency<b>",
                      zeroline = FALSE,
                      ticks="outside",
                      nticks=20,dtick = 0.03),
         yaxis = list(title = "<b>Normalized Frequency<b>",
                      zeroline = TRUE,
                      ticks="outside",
                      nticks=20), bargap = 0.1)
p6

subplot(p5, p6,shareX = TRUE,shareY =TRUE)%>%layout(
  annotations = list( 
    list(x = -0.01 , y = 1.04, font = list(size = 16),text = "<b>A<b>",showarrow = F,xref='paper', yref='paper')),
  margin = list(b = 40, l = 100))



options(scipen = 999)
cmf = cheekVCF[!grepl("NA", cheekVCF$`Mom Freq`),]

View(cmf)
cmf$`Mom Freq` = as.numeric(cmf$`Mom Freq`) 
h3 = hist(cmf$`Mom Freq`,xlim = c(0,0.4),xaxt='n', ylim = c(0,750),breaks=65,yaxt='n',
     main = "B",
     xlab = "Mother Minor Allele Frequencies")
axis(1,at=seq(0,1.0,by=0.01))
axis(2,at=seq(0,750,by=250))

h3$density=h3$counts/sum(h3$counts)
plot(h3,freq=FALSE,xlim = c(0,0.4),xaxt='n',yaxt = 'n',ylim = c(0,0.75),
     main = "B",
     xlab = "Mother Minor Allele Frequencies",ylab = "Normalized Frequency")
axis(1,at=seq(0,1.0,by=0.01))
axis(2,at=seq(0,0.8,by=0.1))

cmf = bloodVCF[!grepl("NA", bloodVCF$`Mom Freq`),]
cmf$`Mom Freq` = as.numeric(cmf$`Mom Freq`) 
View(cmf)
h4 = hist(cmf$`Mom Freq`,xlim = c(0,0.8),xaxt='n',  ylim = c(0,1250),breaks=100,yaxt='n',
     main = "A",
     xlab = "Mother Minor Allele Frequencies")
axis(1,at=seq(0,0.8,by=0.01))
axis(2,at=seq(0,1500,by=250))

h4$density=h4$counts/sum(h4$counts)
plot(h4,freq=FALSE,xlim = c(0,0.4),xaxt='n',yaxt = 'n',ylim = c(0,0.75),
     main = "A",
     xlab = "Mother Minor Allele Frequencies",ylab = "Normalized Frequency")
axis(1,at=seq(0,1.0,by=0.01))
axis(2,at=seq(0,0.8,by=0.1))
