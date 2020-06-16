###author: Hossein Estiri -- hestiri@mgh.harvard.edu
###this script, pulls in the data, does the prep work, and performs different sequencing schemes.
##preparing the data
seed <- 10000
set.seed(seed)
pattern <- "sequent" #"sequent" for traditional sequential mining (SPM) or for tSPM method use "transitive"
count <- "first" #use the first observation or count "all" for SPM


#-------------------------------
Sys.setenv(R_MAX_NUM_DLLS = 999)
options(java.parameters = "-Xmx8048m")
####  Install and load the required packages
if (!require("easypackages")) install.packages('easypackages', repos = "http://cran.rstudio.com/")
packages("data.table","devtools","backports","Hmisc","dplyr","DT","ggplot2","gridExtra","scales","plyr",
         prompt = F)
##first, load the packages
options(scipen = 999)
if (!require("scales")) install.packages('scales', repos = "http://cran.rstudio.com/")
if (!require("reshape2")) install.packages('reshape2', repos = "http://cran.rstudio.com/")
if (!require("foreach")) install.packages('foreach', repos = "http://cran.rstudio.com/")
if (!require("doParallel")) install.packages('doParallel', repos = "http://cran.rstudio.com/")
#---------------------------------


# you need a data model with 3 columns
# patient number | phenx | date
# 
# phenx is the column that stores all features

# this code is based on i2b2 data model's naming convention
# needs name modifications for other data models

# data is stored in a database I call dbmart
uniqpats <- c(unique(dbmart$patient_num)) 

##-----------------------------------------
print(sprintf("sequencing the data"))

###setup parallel backend
cores<-detectCores()
cl <- makeCluster(cores[1]-2) 
registerDoParallel(cl)
###
# remove labs from sequencing
setDT(dbmart) 

if (count == "first"){
  ## only using the first occurence for each observations 
  dbmart.first <- foreach(p = 1: length(uniqpats),
                          .combine = "rbind",
                          .packages = c("plyr")) %dopar% {
                            tryCatch({
                              pat.dat <- subset(dbmart,dbmart$patient_num == uniqpats[p])
                              #store the first observation for each record
                              first.obser.date <- plyr::ddply(pat.dat,~phenx,summarise,start_date=min(start_date))
                              first.obser.date$start_date <- as.POSIXct(first.obser.date$start_date, "%Y-%m-%d")
                              pat.dat$start_date <- as.POSIXct(pat.dat$start_date, "%Y-%m-%d")
                              pat.dat$unique.key <- paste0(pat.dat$phenx,as.character(pat.dat$start_date,"%Y-%m-%d"))
                              first.obser.date$unique.key <- paste0(first.obser.date$phenx,as.character(first.obser.date$start_date,"%Y-%m-%d"))
                              #only grab data from the first observations
                              pat.dat <- subset(pat.dat, pat.dat$unique.key %in% first.obser.date$unique.key)
                              #remove duplicates
                              pat.dat <- pat.dat[!duplicated(pat.dat$unique.key), ]
                              pat.dat$unique.key <- NULL
                              
                              pat.dat
                            },
                            error = function(fr) {cat("ERROR :",conditionMessage(fr), "\n")})
                          }
  rm(pat.dat);gc()
  #simplify the data
  dbmart.first.sim <- dplyr::select(dbmart.first,patient_num,start_date,phenx)
  rm(dbmart.first)
} else if (count != "first"){
  dbmart.first.sim <- dplyr::select(dbmart,patient_num,start_date,phenx)
}

uniqpats <- c(unique(dbmart.first.sim$patient_num)) 


##now sequence by time this process constructs 2-deep sequences 
# by all observations that happened after 1st observation.

dbseq <- list()
pat.dat.sequence <- list()
seq2deep.agg <- list()


if(pattern == "sequent"){
  for (y in 1:length(uniqpats)) {
    tryCatch({
      print(sprintf(paste0("sequencing data for patient ",y, " of ", length(uniqpats))))
      pat.dat <- subset(dbmart.first.sim,dbmart.first.sim$patient_num == uniqpats[y])
      pat.dat$start_date <- as.POSIXct(pat.dat$start_date, "%Y-%m-%d")
      pat.dat.date <- as.data.table(unique(pat.dat$start_date))
      pat.dat.date[with(pat.dat.date ,order(x)),sequence := .I]
      colnames(pat.dat.date)[1] <- "start_date"
      pat.dat.date$start_date <- as.POSIXct(pat.dat.date$start_date, "%Y-%m-%d")
      pat.dat<- merge(pat.dat,pat.dat.date,by="start_date")
      rm(pat.dat.date)
      max.pat.obs <- max(pat.dat$sequence)
      dbseq[[y]] <- data.frame(foreach(g = 1: (max.pat.obs-1), 
                                       .combine = "rbind") %dopar% {
                                         seq.0 <- subset(pat.dat,pat.dat$sequence == g)
                                         seq.1 <- subset(pat.dat,pat.dat$sequence == g+1)
                                         seq.1$sequence.real=seq.1$sequence
                                         seq.1$sequence=g
                                         pat.dat.sequence <- merge(seq.0,seq.1,by="sequence",allow.cartesian=TRUE)
                                         rm(seq.0,seq.1)
                                         pat.dat.sequence},stringsAsFactors=FALSE)
      pat.dat.sequence = list()
      rm(pat.dat,max.pat.obs)
      if(endsWith(as.character(y),"500")){gc()}
      
    }, 
    error = function(fr) {cat("ERROR :",conditionMessage(fr), "\n")})
  } 
} else if(pattern != "sequent"){
  for (y in 1:length(uniqpats)) {
    tryCatch({
      print(sprintf(paste0("sequencing data for patient ",y, " of ", length(uniqpats))))
      pat.dat <- subset(dbmart.first.sim,dbmart.first.sim$patient_num == uniqpats[y])
      pat.dat$start_date <- as.POSIXct(pat.dat$start_date, "%Y-%m-%d")
      pat.dat.date <- as.data.table(unique(pat.dat$start_date))
      pat.dat.date[with(pat.dat.date ,order(x)),sequence := .I]
      colnames(pat.dat.date)[1] <- "start_date"
      pat.dat.date$start_date <- as.POSIXct(pat.dat.date$start_date, "%Y-%m-%d")
      pat.dat<- merge(pat.dat,pat.dat.date,by="start_date")
      rm(pat.dat.date)
      max.pat.obs <- max(pat.dat$sequence)
      dbseq[[y]] <- data.frame(foreach(g = 1: (max.pat.obs-1), 
                                       .combine = "rbind") %dopar% {
                                         seq.0 <- subset(pat.dat,pat.dat$sequence == g)
                                         seq.1 <- subset(pat.dat,pat.dat$sequence > g)
                                         seq.1$sequence.real=seq.1$sequence
                                         seq.1$sequence=g
                                         pat.dat.sequence <- merge(seq.0,seq.1,by="sequence",allow.cartesian=TRUE)
                                         rm(seq.0,seq.1)
                                         pat.dat.sequence},stringsAsFactors=FALSE)
      
      rm(pat.dat,max.pat.obs,tmp)
      pat.dat.sequence = list()
      if(endsWith(as.character(y),"500")){gc()}
      
    }, 
    error = function(fr) {cat("ERROR :",conditionMessage(fr), "\n")})
  }
}
rm(pat.dat,max.pat.obs)
seq2deep <- data.table::rbindlist(dbseq)
rm(dbseq,pat.dat.sequence);gc()

setDT(seq2deep)
seq2deep[,row := .I]
seq2deep$value.var <- 1
seq2deep <- data.frame(seq2deep)
# seq2deep[,row := .I]
seq2deep$patient_num.y <- NULL
colnames(seq2deep)[3] <- "patient_num" #change this to your patient number column name
# now reducing the 2-deep sequence data to the ones that have more than 1 frequency
seq2deep$seq2deep <- as.character(paste0(seq2deep$phenx.x," -> ",seq2deep$phenx.y))

rm(cl,ablab,dbmart.first.sim,i,covdate,y,cores);gc()

# aggregating by patient number
library(dplyr)
seq2deep.agg <- seq2deep %>% 
  dplyr::group_by(seq2deep) %>%
  dplyr::summarise(distinct_patients=n_distinct(patient_num,na.rm=TRUE))

uniqpats <- c(unique(seq2deep$patient_num)) 

seq2deep.agg.cut <- subset(seq2deep.agg,seq2deep.agg$distinct_patients > round(length(uniqpats)/500)) ##>0.05%
sequences <- c(as.character(unique(seq2deep.agg.cut$seq2deep)))
rm(seq2deep.agg.cut)
gc()
seq2deep2 <- seq2deep[(seq2deep$seq2deep %in% sequences),]
seq2deep2$start_date.y <- seq2deep2$patient_num.y <- NULL
colnames(seq2deep2)[3] <- "patient_num"
setDT(seq2deep2)
seq2deep2[,row := .I]
seq2deep2$value.var <- 1
seq2deep2 <- data.frame(seq2deep2)
seqtrans <- dplyr::select(seq2deep2,seq2deep,phenx.x,phenx.y)
seqtrans <-seqtrans[!duplicated(paste0(seqtrans$seq2deep,seqtrans$phenx.x,seqtrans$phenx.y)), ]
uniqpats <- c(unique(seq2deep2$patient_num)) 
