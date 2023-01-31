
# @author : Mª Eugenia Soler (maria.soler@cabimer.es / eusololi@gmail.com)
# @brief : feature_hic_metaplot.R generates HiC data metaplots around a any feature of interest. The input must be a .hic matrix, the output that comes
#from applying the Juicer pipeline (Aiden Lab) to a HiC experiment. straw R package from Aiden Lab (https://github.com/aidenlab/straw/tree/master/R)
#and HiTC package (Servant N, Lajoie BR, Nora EP, Giorgetti L, Chen C, Heard E, Dekker J, Barillot E (2012). “HiTC: Exploration of High-Throughput 'C' 
#experiments.” Bioinformatics. doi: 10.1093/bioinformatics/bts521. were used to read and process the input.hic file.

#Packages needed:
library("Matrix")
library("strawr")
library("Rcpp")
library("HiTC")
options(scipen=999)
chromosomes<-c(paste("chr", 1:22, sep="" ))

feature_hic_metaplot<- function(hic_input_file, unit, binsize, feature_coordinates_file, margins, metahic_matrix_file)
{
import_matrix<-data.frame()
xgi_matrix<-data.frame()
index<-1
for(i in 1:length(chromosomes)
  {
  #Loading each chromosome matrix, using straw package
  hic.data.frame <- strawr::straw(norm="NONE", hic_input_file, as.character(i), as.character(i), unit, binsize, matrix="observed") 
  #Creating a table with the chromosome, coordinates and the names of the bins (xgi_table); 
  #and a matrix with the bin names and contacts between them (import_matrix)
  sorted_hic.data.frame <- hic.data.frame[order(hic.data.frame$x, hic.data.frame$y),]
  start_bins<-seq(1,max(sorted_hic.data.frame$x, sorted_hic.data.frame$y)-9999,by=10000)
  end_bins<-seq(10000,max(sorted_hic.data.frame$x, sorted_hic.data.frame$y),by=10000)
  if(i==1)
    {
    
    bins_names<-seq(index,length(start_bins), by=1)
    }
  else
    {
     index<-index+nrow(bins_coordinates)+1
     bins_names<-seq(index,length(start_bins)+index-1, by=1)
    }
  bin1<-as.vector(sorted_hic.data.frame[,1])
  bin2<-as.vector(sorted_hic.data.frame[,2])
  bin1_name<-c()
  bin2_name<-c()
  m<-1
  for(j in 1:nrow(sorted_hic.data.frame))
    {
    if(bin1[j]==0)
      {
      sorted_hic.data.frame<-sorted_hic.data.frame[-j,]
      }
    else
      {
      bin1_name[m]<-bins_names[which(end_bins==bin1[j])]
      bin2_name[m]<-bins_names[which(end_bins==bin2[j])]
      m<-m+1 
      }
  }
  bins<-cbind(rep(i, length(length(start_bins))),start_bins, end_bins, bins_names)
  sorted_hic.data.frame <-cbind(bin1_name, bin2_name,sorted_hic.data.frame[,3],)
  #Adding the information of new chromosome to import_matrix and xgi_table
  import_matrix <- rbind(import_matrix, sorted_hic.data.frame)
  xgi_table <- rbind(xgi_matrix, bins)
}

matrix<-importC(import_matrix , xgi_matrix, ygi.bed=NULL, allPairwise=FALSE, rm.trans=TRUE)
rm(import_matrix)
feature_coordinates<-read.table(file=feature_coordinates_file)
feature_coordinates <- feature_coordinates[order(feature_coordinates[,1], feature_coordinates[,2], feature_coordinates[,3]),]

##feature midpoint calculation and place of the features in its corresponding bin
midpoint<-c()
for(i in 1:nrow(feature_coordinates))
{
midpoint[i]<-(feature_coordinates[i,3]-feature_coordinates[i,2])/2 + feature_coordinates[i,2]
}
midpoint_dataframe<-cbind(feature_coordinates[,1], midpoint, feature_coordinates[,4:6])

feature_in_bin<-data.frame()
m <- 1
for(i in 1:length(chromosomes))
{
 chr_bins<-subset(xgi_matrix, xgi_matrix[,1]==chromosomes[i])
 chr_feature<-subset(midpoint_dataframe, midpoint_dataframe[,1]==chromosomes[i])
 if(nrow(chr_feature)>1)
  {
  for(j in 1:nrow(chr_feature))
    {
    for(k in 1:nrow(chr_bins))
      {
      if ((as.numeric(chr_feature[j,2]) >= as.numeric(chr_bins[k,2])) && (as.numeric(chr_feature[j,2]) <= as.numeric(chr_bins[k,3])))
        {
        feature_in_bin[m,1]<-as.character(chr_feature[j,1])
        feature_in_bin[m,2]<-as.numeric(chr_feature[j,2])
        feature_in_bin[m,3:5]<-as.numeric(chr_bins[k,2:4])
        feature_in_bin[m,6]<-as.character(chr_feature[j,5])
        m<-m+1
        }
      }
    }
  }
}
##Creating two meta-HiC-matrices, depending on the strand of the feature
rloop_matrix_plus <- matrix(nrow=(margins/binsize)*2+1, ncol=(margins/binsize)*2+1, 0)
rloop_matrix_minus <- matrix(nrow=(margins/binsize)*2+1, ncol=(margins/binsize)*2+1, 0)
##In each feature coordinates +- the margins selected, extraction and normalization (o/e) of its contacts data.
htclist<-do.call(c, matrix)
rm(matrix)
m<-0
n<-0
for(i in 1:22)
{
print(i)
chr_bins<-subset(xgi_matrix, xgi_matrix[,1]==chromosomes[i])
chr_feature<-subset(feature_in_bin, feature_in_bin[,1]==chromosomes[i])
for(j in 1:nrow(chr_feature))
  {
   print(j)
   if((chr_feature[j,3]-margins)>0 && chr_feature[j,4]<max(chr_bins[,2:3]))
    {
    rm(rloop)      
    rloop <- extractRegion(htclist[[i]], c(1,2), chr=chr_feature[j,1], from=(chr_feature[j,3]-margins), to=(chr_feature[j,4]+margins))
    rm(rloop_normperexp)
    rloop_normperexp <- normPerExpected(rloop)
    rm(rloop_matrix_new)
    rloop_matrix_new <- as.matrix(rloop_normperexp@intdata)
    rloop_matrix_new[is.na(rloop_matrix_new)] <- 0
    if(nrow(rloop_matrix_new)==(margins/binsize)*2+1 && ncol(rloop_matrix_new)==(margins/binsize)*2+1)
      {
      if(feature_in_bin[j,6]=="+")
        {
        rloop_matrix_plus <- rloop_matrix_plus + rloop_matrix_new
        m<-m+1
        }
      else
        {
        rloop_matrix_minus <- rloop_matrix_minus + rloop_matrix_new
        n<-n+1
        }
      }
    }
  }
}

##Invertion of the minus strand matrix in order to calculate the average of plus and minus features HiC maps together
rloop_matrix_minus_new <- matrix(nrow=(margins/binsize)*2+1, ncol=(margins/binsize)*2+1, 0)
j <- nrow(rloop_matrix_minus_new)
for(i in 1:ncol(rloop_matrix_minus))
  {
  rloop_matrix_minus_new[,j]<-rev(rloop_matrix_minus[,i])
  j <- j-1
  }
t_matrix_minus_new <- t(rloop_matrix_minus_new)
#Average calculation
rloop_matrix_both <- t_matrix_minus_new+rloop_matrix_plus 
rloop_matrix_both <- rloop_matrix_both/(m+n) ##Dividimos entre el número total de regiones (+ y -)
#Output file
write.table(rloop_matrix_both, file=metahic_matrix_file, quote=FALSE, sep=",")
}
