install.packages("plotly")
install.packages("heatmaply")
install.packages("ggcorrplot")
install.packages("Rcpp")
library("heatmaply")
library("strawr")
library("HiTC")
library("Matrix")
library("Rcpp")
install.packages("HiTC")
remotes::install_github("aidenlab/straw/R")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Matrix")

##  SCIPEN HAY QUE HACERLO SI O SI PARA EVITAR HERRORES POSTERIORES
options(scipen=999)

##Cargamos con straw el subset de la matriz que queramos, con la normalizacion y el tamaño de bin deseado
hic.data.frame <- strawr::straw("NONE", "ENCFF406HHC.hic", "1", "1", "BP", 10000)
head(hic.data.frame)
dim(hic.data.frame)
##Ordenamos la matriz y miramos la ultima fila para hacernos una idea de cual es el ultimo bin
sorted_hic.data.frame <- hic.data.frame[order(hic.data.frame$x, hic.data.frame$y),]
head(sorted_hic.data.frame)
dim(sorted_hic.data.frame)
sorted_hic.data.frame[1027960,]
##Para darle nombre a los bins, como en la matriz por defecto se llaman como el bp final del bin, lo dividimos entre 10000
##(el tamaño del bin) para tener numeros más pequeños, (1,2,3,4,...)
bins_x <- sorted_hic.data.frame[,1]/10000
bins_x[1:10]
bins_y <- sorted_hic.data.frame[,2]/10000
bins_y[1:10]
##Si no, para saber cual es el bin mayor, también podemos hacer esto en lugar de ver la fila 
max(bins_y)

##Creamos una nueva matriz con los nuevos nombres de las bins, y las counts de la matriz. Esa será nuestra matriz de contacto
sorted_bins_hic.data.frame <- data.frame(bins_x, bins_y, sorted_hic.data.frame[,"counts"])
head(sorted_bins_hic.data.frame)
dim(sorted_bins_hic.data.frame)
write.table(sorted_bins_hic.data.frame, file="contact_matrix_22.csv", sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

## Al último bin, le restamos 10000 y lo que nos dé +1 será el último valor de coordenada del start bins
59110000 - 10000
start_bins <- seq(1,51220001, by=10000)
length(start_bins)
end_bins <- seq(10000,51230000, by=10000)
bins_name <- seq(1,5123, by=1)
chr <- rep("chr22", 5123)
##Hacemos la tabla de coordenadas de los bins
bins_coordinates <- data.frame (chr, start_bins, end_bins, bins_name)
is.atomic(sorted_bins_hic.data.frame)
head(bins_coordinates)
write.table(bins_coordinates, file="bins_coordinates_22.bed", sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE, )

##Importamos la matriz de contacto con la tabla que explica las coordenadas de cada bin
chr1 <- importC(con="contact_matrix_22.csv" ,  xgi.bed ="bins_coordinates_22.bed", ygi.bed=NULL, allPairwise=FALSE, rm.trans=TRUE)
CQC(chr16)
chr1$chr1chr1
mapC(chr1$chr1chr1)

##Filtramos las coordenadas de la característica que vamos a mirar (en este caso, conflictos CD o HO) por cromosoma
HO <- read.table(file="Galaxy669-[Cut_on_data_667_HEAD-ON_c4_DRIPc_signal_c5_Rloop_peak_length].tabular")
head(HO)
HO_1<-data.frame()
j <- 1
for (i in 1:nrow(HO))
{
  print(i)
  if (HO[i,1]=="chr22")
  {
  HO_1[j,1] <- HO[i,1]
  HO_1[j,2] <- (HO[i,3] - HO[i,2])/2 + HO[i,2]
  j <- j+1
  }
}
head(HO_1)
dim(HO_1)


head(bins_coordinates)
##Situamos cada coordenada de R loop en su bin correspondiente
rloop_in_bin <- data.frame()
k <- 1
for(i in 1:nrow(HO_1))
{
  print(i)
  for(j in 1:nrow(bins_coordinates))
  {
    if ((as.numeric(HO_1[i,2]) >= as.numeric(bins_coordinates[j,2])) && (as.numeric(HO_1[i,2]) <= as.numeric(bins_coordinates[j,3])))
    {
      rloop_in_bin[k,1] <- HO_1[i,1]
      rloop_in_bin[k,2] <- HO_1[i,2]
      rloop_in_bin[k,3] <- bins_coordinates[j,2]
      rloop_in_bin[k,4] <- bins_coordinates[j,3]
      rloop_in_bin[k,5] <- bins_coordinates[j,4]
      k <- k+1
    }
  }
}
dim(rloop_in_bin)
head(rloop_in_bin)
head(bins_coordinates)

##Creamos una matriz vacía donde se guardarán los datos del meta-HiC
rloop_matrix <- matrix(nrow=401, ncol=401, 0)
dim(rloop_matrix_new)

##Para cada Rloop, extraemos y normalizamos sus datos de contacto (cambiando los NA por ceros pa que se puedan sumar)
##Vamos sumando todos de cada cromosoma, y después lo dividiremos por el número total para calcular la media
for(i in 1:nrow(rloop_in_bin))
{
  print(i)
  rloop <- extractRegion(chr1$chr21chr21, c(1,2), chr="chr21", from=(rloop_in_bin[i,3]-2000000), to=(rloop_in_bin[i,4]+2000000))
  rloop_normperexp <- normPerExpected(rloop)
  rloop_matrix_new <- as.matrix(rloop_normperexp@intdata)
  rloop_matrix_new[is.na(rloop_matrix_new)] <- 0
  rloop_matrix <- rloop_matrix + rloop_matrix_new
}
rloop_matrix[1:5,1:5]

##Dividmos entre el número total y tenemos la matriz con el valor medio de contactos
rloop_matrix<- rloop_matrix/74
write.table(rloop_matrix, file="rloop_cd_metahic_minus.tsv")


##Las volvemos a multiplicar por el número total, las de cada cadena, para poder sumarlas

rloop_matrix_plus <- read.table(file="rloop_cd_metahic_plus.tsv")
rloop_matrix_plus <- rloop_matrix_plus*69

rloop_matrix_minus <- read.table(file="rloop_cd_metahic_minus.tsv")
rloop_matrix_minus <- rloop_matrix_minus*74

##Hacemos una cosa loquísima que no sé ni como sale pero sale, para invertir las matrices y poder sumarlas
rloop_matrix_minus_new <- matrix(nrow=401, ncol=401)
j <- 401
for(i in 1:ncol(rloop_matrix_minus))
{
  print(i)
  rloop_matrix_minus_new[,j]<-rev(rloop_matrix_minus[,i])
  j <- j-1
}
rloop_matrix_minus[1:5,1:5]
rloop_matrix_minus_new[397:401,397:401]
test2 <- t(rloop_matrix_minus_new)
test2[397:401,397:401]
write.table(rloop_matrix_minus_new, file="rloop_ho_metahic_minus_good.tsv") ##Esta es la matriz invertida de la cadena -, para poder sumarla a la +


rloop_matrix_both <- test2+rloop_matrix_plus ##Sumamos la matriz invertida (-) a la positiva (+)
rloop_matrix_both <- rloop_matrix_both/143 ##Dividimos entre el número total de regiones (+ y -)

write.table(rloop_matrix_both, file="rloop_cd_metahic_both.tsv")

##FC
##Cargamos las dos matrices
ho <- read.table("rloop_ho_metahic_both.csv", sep=";", header=TRUE)
dim(ho)
ho<-ho[,2:402]
ho[1:5,1:5]
cd <- read.table(file="rloop_cd_metahic_both.csv",sep=";", header=TRUE)
dim(cd)  
cd<-cd[,2:402]
cd[1:5,1:5]
ho <- as.matrix(ho)
cd <- as.matrix(cd)

fold <- ho/cd
fold[1:5,1:5]
log2fold <- log2(fold)
log2fold[1:5,1:5]

write.table(log2fold, file="log2fold_metahic_hovscd.csv", sep=",", row.names = TRUE, col.names = TRUE)

##Sacar regiones especificas
HO <- read.table(file="Galaxy771-[Cut_on_data_766_CD_CONFLICTS_c4_DRIPc_intensity_c5_length_c6_strand].bed")
chr1 <- importC(con="contact_matrix_11.csv" ,  xgi.bed ="bins_coordinates_11.bed", ygi.bed=NULL, allPairwise=FALSE, rm.trans=TRUE)
HO[276,]
rloop <- extractRegion(chr1$chr11chr11, c(1,2), chr="chr11", from=(HO[276,2]-2000000), to=(HO[276,3]+2000000))
rloop_normperexp <- normPerExpected(rloop)
rloop_matrix_new <- as.matrix(rloop_normperexp@intdata)
rloop_matrix_new[is.na(rloop_matrix_new)] <- 0
rloop_matrix_new[1:5,1:5]
write.table(rloop_matrix_new, file="CD_276.csv", sep=",", row.names=TRUE, col.names=TRUE)
mapC(rloop_normperexp,minrange=0, maxrange=1.5, col.pos=c("blue", NA , "red"), col.na="black")
#####################################################################################################
##Con las matrices iced del paper Barutcu et al smarca4
##leemos la matriz, vemos cual es el ultimo bin, la rehacemos para quitar los nombres de filas y columnas, y viendo las bins creamos
##la tabla de coordenadas
chr <- read.table(file="HiCStein-MCF10a-shBRG1__hg19__chr1__C-40000-iced.matrix", sep="\t")
  dim(chr)
chr[1:5,1:5]
chr[3883,1]
chr <- chr[2:6081,2:6081]
start_bins <- seq(1,155240001, by=40000)
length(start_bins)
155240001+40000
end_bins <- seq(40000,155280001, by=40000)
length(end_bins)
bins_name <- seq(1,3882, by=1)
chr <- rep("chrX", 3882)

bins_coordinates <- data.frame (chr, start_bins, end_bins, bins_name)
head(bins_coordinates)
write.table(bins_coordinates, file="bins_coordinates_chrX_40kb.bed", sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE, )

##Filtramos las coordenadas de la característica que vamos a mirar (en este caso, conflictos CD o HO) por cromosoma
bins_coordinates <- read.table(file="bins_coordinates_chr1_40kb.bed")
head(bins_coordinates)
HO <- read.table(file="brg1_rloopgain_hg19.bed")
head(HO)
HO_1<-data.frame()
j <- 1
for (i in 1:nrow(HO))
{
  print(i)
  if (HO[i,1]=="chrX")
  {
      if (HO[i,6]=="+")
      {
    HO_1[j,1] <- HO[i,1]
    HO_1[j,2] <- (HO[i,3] - HO[i,2])/2 + HO[i,2]
    j <- j+1
    }
  }
}
head(HO_1)
dim(HO_1)

write.table(HO_1, file="brg1gain_chrX_plus_hg19.bed", sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE, )
##Situamos cada coordenada de R loop en su bin correspondiente
    ##Creamos una matriz vacía donde se guardarán los datos del meta-HiC
    rloop_matrix <- matrix(nrow=101, ncol=101, 0)
    dim(rloop_matrix)

chr <- read.table(file="HiCStein-MCF10a-shBRG1__hg19__chrX__C-40000-iced.matrix")
dim(chr)
chr <- chr[2:nrow(chr),2:nrow(chr)]
chr[1:5,1:5]
bins_coordinates <- read.table(file="bins_coordinates_chrX_40kb.bed")
HO_1 <- read.table(file="brg1gain_chrX_plus_hg19.bed", sep="\t")
head(HO_1)
rloop_in_bin <- data.frame()
k <- 1
for(i in 1:nrow(HO_1))
{
  print(i)
  for(j in 1:nrow(bins_coordinates))
  {
    if ((as.numeric(HO_1[i,2]) >= as.numeric(bins_coordinates[j,2])) && (as.numeric(HO_1[i,2]) <= as.numeric(bins_coordinates[j,3])))
    {
      rloop_in_bin[k,1] <- HO_1[i,1]
      rloop_in_bin[k,2] <- HO_1[i,2]
      rloop_in_bin[k,3] <- bins_coordinates[j,2]
      rloop_in_bin[k,4] <- bins_coordinates[j,3]
      rloop_in_bin[k,5] <- bins_coordinates[j,4]
      k <- k+1
    }
  }
}
dim(rloop_in_bin)
head(rloop_in_bin)
head(bins_coordinates)

##Para cada Rloop, extraemos y normalizamos sus datos de contacto (cambiando los NA por ceros pa que se puedan sumar)
##Vamos sumando todos de cada cromosoma, y después lo dividiremos por el número total para calcular la media
k<- 0
for(i in 1:nrow(rloop_in_bin))
{
  print(i)
    if((rloop_in_bin[i,5]>50)&&((rloop_in_bin[i,5]+50)<max(bins_coordinates[,4])))
    {
  bin <- rloop_in_bin[i,5]
  rloop <- chr[(bin-50):(bin+50), (bin-50):(bin+50)]
  rloop_matrix_new <- apply(rloop, 2, as.numeric)
  rloop_matrix_new[is.na(rloop_matrix_new)] <- 0
  rloop_matrix <- rloop_matrix + rloop_matrix_new
  k <- k+1
  }
}
k
rloop_matrix[1:5,1:5]
rloop_matrix_new[1:5,1:5]
dim(rloop_matrix)
rloop_in_bin[5,]
dim(bins_coordinates)
bins_coordinates[3530,]
##Dividmos entre el número total y tenemos la matriz con el valor medio de contactos

130+78+70+30+42+66+77+36+53+44+98+56+15+37+48+94+86+7+150+28+16+53+35 ##brg1 gain 1349

72+56+70+51+69+92+46+47+18+36+39+37+16+41+32+25+19+20+18+8+11+12+22 ##ho - 857

13+19+11+6+9+13+6+7+12+7+10+2+5+15+2+5+11+4+6+5+3+3+4  ##cd + 178
rloop_matrix<- rloop_matrix/1349
write.table(rloop_matrix, file="rloop_brg1gain_plus_shbrg1_metahic_hg19_log10.tsv", sep="\t")
rloop_matrix<- log10(rloop_matrix)


##Las volvemos a multiplicar por el número total, las de cada cadena, para poder sumarlas

rloop_matrix_wt <- read.table(file="rloop_brg1gain_plus_shbrg1_metahic_hg19.tsv")
rloop_matrix_brg1 <-read.table(file="rloop_brg1gain_plus_shscram_metahic_hg19.tsv") 

rloop_matrix_minus <- read.table(file="rloop_cd_metahic_minus.tsv")
rloop_matrix_minus <- rloop_matrix_minus*74

##Hacemos una cosa loquísima que no sé ni como sale pero sale, para invertir las matrices y poder sumarlas
rloop_matrix_minus_new <- matrix(nrow=401, ncol=401)
j <- 401
for(i in 1:ncol(rloop_matrix_minus))
{
  print(i)
  rloop_matrix_minus_new[,j]<-rev(rloop_matrix_minus[,i])
  j <- j-1
}
rloop_matrix_minus[1:5,1:5]
rloop_matrix_minus_new[397:401,397:401]
test2 <- t(rloop_matrix_minus_new)
test2[397:401,397:401]
write.table(rloop_matrix_minus_new, file="rloop_ho_metahic_minus_good.tsv") ##Esta es la matriz invertida de la cadena -, para poder sumarla a la +


rloop_matrix_both <- test2+rloop_matrix_plus ##Sumamos la matriz invertida (-) a la positiva (+)
rloop_matrix_both <- rloop_matrix_both/143 ##Dividimos entre el número total de regiones (+ y -)

write.table(rloop_matrix_both, file="rloop_cd_metahic_both.tsv")

##FC
##Cargamos las dos matrices
ho <- read.table("rloop_brg1gain_plus_shbrg1_metahic_hg19_log10.tsv", sep="\t", header=TRUE)
dim(ho)
ho<-ho[,2:402]
ho[1:5,1:5]
cd <- read.table(file="rloop_brg1gain_plus_shscram_metahic_hg19_log10.tsv",sep="\t", header=TRUE)
dim(cd)  
cd<-cd[,2:402]
cd[1:5,1:5]
ho <- as.matrix(ho)
cd <- as.matrix(cd)

fold <- ho/cd
fold[1:5,1:5]
log2fold <- log2(fold)
log2fold[1:5,1:5]

write.table(log2fold, file="log2fold_metahic_brg1gain_brg1vswt_log10.tsv", sep="\t", row.names = TRUE, col.names = TRUE)
