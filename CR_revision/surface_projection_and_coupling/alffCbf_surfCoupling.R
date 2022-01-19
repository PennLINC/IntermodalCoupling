#------------------------#
# ALFF-CBF SURF COUPLING #
#------------------------#

# read in packages
library(plyr)

# read in demos
alffCbf_subjects <- read.csv("/data/jux/BBL/projects/coupling/subjectsLists/n1400_bblid_datexscanid.csv")
demos <- read.csv("/data/jux/BBL/projects/coupling/subjectsLists/n1601_demographics_go1_20161212.csv")
cbfQa <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/asl/n1601_PcaslQaData_20170403.csv")
restQa <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/rest/n1601_RestQAData_20170714.csv")

# subset csvs
demos <- subset(demos, select = c("bblid","scanid","ageAtScan1","sex"))
cbfQa <- subset(cbfQa, select = c("bblid","scanid","pcaslRelMeanRMSMotion"))
restQa <- subset(restQa, select = c("bblid","scanid","restRelMeanRMSMotion"))

# join csvs
alffCbf_subjDemos_orig <- join(alffCbf_subjects,demos)
alffCbf_subjDemos_cbfQa <- join(alffCbf_subjDemos_orig,cbfQa)
alffCbf_subjDemos_cbfQa_restQa <- join(alffCbf_subjDemos_cbfQa,restQa)
alffCbf_subjDemos <- as.data.frame(alffCbf_subjDemos_cbfQa_restQa)
alffCbf_subjDemos <- alffCbf_subjDemos[complete.cases(alffCbf_subjDemos),]
alffCbf_subjDemos$ageAtScan1 <- (alffCbf_subjDemos$ageAtScan1)/12
alffCbf_subjDemos$sex <- as.factor(alffCbf_subjDemos$sex)

# get all files
setwd("/data/jux/BBL/projects/coupling/couplingSurfaceMaps/alffCbf/lh/stat")
lh_alffCbf_files = list.files(pattern="*.asc")
lh_alffCbf_data = do.call(rbind, lapply(lh_alffCbf_files, function(x) read.table(x, stringsAsFactors = FALSE)))
lh_alffCbf_coupling <- as.data.frame(lh_alffCbf_data$V5)
lh_alffCbf_coupling_n <- t(as.data.frame(split(lh_alffCbf_coupling,1:1400)))

# get all files
setwd("/data/jux/BBL/projects/coupling/couplingSurfaceMaps/alffCbf/rh/stat")
rh_alffCbf_files = list.files(pattern="*.asc")
rh_alffCbf_data = do.call(rbind, lapply(rh_alffCbf_files, function(x) read.table(x, stringsAsFactors = FALSE)))
rh_alffCbf_coupling <- as.data.frame(rh_alffCbf_data$V5)
rh_alffCbf_coupling_n <- t(as.data.frame(split(rh_alffCbf_coupling,1:1400)))

# run model
for (i in 1:10242) {
  lh_alffCbf_agemodel <- lm(lh_alffCbf_coupling_n[,i] ~ ageAtScan1, data=alffCbf_subjDemos)
  lh_alffCbf_ageSexmodel <- lm(lh_alffCbf_coupling_n[,i] ~ ageAtScan1 + sex, data=alffCbf_subjDemos)
  lh_alffCbf_ageSex_qaModel <- lm(lh_alffCbf_coupling_n[,i] ~ ageAtScan1 + sex + pcaslRelMeanRMSMotion + restRelMeanRMSMotion, data=alffCbf_subjDemos)
}

for (i in 1:10242) {
  rh_alffCbf_agemodel <- lm(rh_alffCbf_coupling_n[,i] ~ ageAtScan1, data=alffCbf_subjDemos)
  rh_alffCbf_ageSexmodel <- lm(rh_alffCbf_coupling_n[,i] ~ ageAtScan1 + sex, data=alffCbf_subjDemos)
  rh_alffCbf_ageSex_qaModel <- lm(rh_alffCbf_coupling_n[,i] ~ ageAtScan1 + sex + pcaslRelMeanRMSMotion + restRelMeanRMSMotion, data=alffCbf_subjDemos)
}

# print out summary results
summary(lh_alffCbf_agemodel)
summary(lh_alffCbf_ageSexmodel)
summary(lh_alffCbf_ageSex_qaModel)
summary(rh_alffCbf_agemodel)
summary(rh_alffCbf_ageSexmodel)
summary(rh_alffCbf_ageSex_qaModel)
