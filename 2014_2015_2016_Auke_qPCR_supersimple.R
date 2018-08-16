
library(broom)
library(cowplot)
library(MASS)
library(DHARMa)
library(tidyverse)


setwd("/Users/taaltree/Dropbox/Taal/eDNAproject/2014_2015_2016_Auke_qPCR")
#Read data
Sockeyelong_all <- read_tsv("Sockeyelong_all.tsv")
Sockeyelong_all$Year <- as.factor(Sockeyelong_all$Year)
Coholong_all <- read_tsv("Coholong_all.tsv")
Coholong_all$Year <- as.factor(Coholong_all$Year)

# calculate a new Q_cfs value from scratch
a <-  0.148 # from values calculated in gage_to_Q_Josh_R/
b <-  19.27 # from values calculated in gage_to_Q_Josh_R/ 
c <-  5.236 # from values calculated in gage_to_Q_Josh_R/
slope=-3.329
int=23.614

# Q_cfs = a (Gage_Height - b) ^ c
Sockeyelong_all <- Sockeyelong_all %>% mutate(Q_cfs2 = a*(Gage_Height - b)^c) 
Sockeyelong_all <- Sockeyelong_all %>% mutate(newQ = 10^((CT_mean-int)/slope)) 

Coholong_all <- Coholong_all %>% mutate(Q_cfs2 = a*(Gage_Height - b)^c)

# calculate a new Qcorr_qPCR value from QUANT_mean * Q_cfs2 # because old Qcorr_qPCR values used different translations of gauge height to Qcfs in different years
# this only matters for plotting since the models use QUANT_mean * Q_cfs2 directly
Sockeyelong_all <- Sockeyelong_all %>% mutate(Qcorr_qPCR = QUANT_mean * Q_cfs2) 
Sockeyelong_all <- Sockeyelong_all %>% mutate(newQcorr_qPCR = newQ * Q_cfs2) 

Coholong_all <- Coholong_all %>% mutate(Qcorr_qPCR = QUANT_mean * Q_cfs2)
#Coholong_all <- Coholong_all %>% spread(Cohotype,Count) %>% mutate(Coho_nojacks= Coho_Total-Coho_Jack) %>% gather("Cohotype","Count",19:26)

#####################
#Begin: Subset data to appropriate timescale
#####################
sockeye.adult.start = "06-18"
sockeye.adult.end = "08-01"
coho.adult.start = "08-18"  # first appearance of 1 Coho male and 1 Coho female
coho.adult.end = "10-21"  # 10-21 is last date of 

sockeye.juv.start = "04-15"
sockeye.juv.end = "06-10"
coho.juv.start = "04-15"
coho.juv.end = "06-10"

sockeye_adult <- Sockeyelong_all %>% dplyr::filter(
    (Date >= as.Date(paste0("2014-",sockeye.adult.start)) & Date <= as.Date(paste0("2014-",sockeye.adult.end)))|
        (Date >= as.Date(paste0("2015-",sockeye.adult.start)) & Date <= as.Date(paste0("2015-",sockeye.adult.end)))|
        (Date >= as.Date(paste0("2016-",sockeye.adult.start)) & Date <= as.Date(paste0("2016-",sockeye.adult.end))))

coho_adult <- Coholong_all %>% dplyr::filter(
    (Date >= as.Date(paste0("2014-",coho.adult.start)) & Date <= as.Date(paste0("2014-",coho.adult.end)))|
        (Date >= as.Date(paste0("2015-",coho.adult.start)) & Date <= as.Date(paste0("2015-",coho.adult.end)))|
        (Date >= as.Date(paste0("2016-",coho.adult.start)) & Date <= as.Date(paste0("2016-",coho.adult.end))))

sockeye_juv <- Sockeyelong_all %>% dplyr::filter(
    (Date >= as.Date(paste0("2014-",sockeye.juv.start)) & Date <= as.Date(paste0("2014-",sockeye.juv.end)))|
        (Date >= as.Date(paste0("2015-",sockeye.juv.start)) & Date <= as.Date(paste0("2015-",sockeye.juv.end)))|
        (Date >= as.Date(paste0("2016-",sockeye.juv.start)) & Date <= as.Date(paste0("2016-",sockeye.juv.end))))

coho_juv <- Coholong_all %>% dplyr::filter(
    (Date >= as.Date(paste0("2014-",coho.juv.start)) & Date <= as.Date(paste0("2014-",coho.juv.end)))|
        (Date >= as.Date(paste0("2015-",coho.juv.start)) & Date <= as.Date(paste0("2015-",coho.juv.end)))|
        (Date >= as.Date(paste0("2016-",coho.juv.start)) & Date <= as.Date(paste0("2016-",coho.juv.end))))

#####################
#End: Subset data to appropriate timescale
#####################

#####################
#Begin: Timeline Function
#####################

timelines_full <- function(salmonspecies, start = "2015-01-01", end = "2015-12-31", type = "_adult") { # x == Sockeye, Coho, y == _adult, _juv
    
            salmon_use <- get(paste0(salmonspecies, type)) # Sockeye_adult or Coho_adult
            
            salmon_use <- salmon_use %>% dplyr::filter(Date >= as.Date(start) & Date <= as.Date(end))  # truncate start and stop times
            
            Salmon_type <- str_to_title(paste0(salmonspecies, "type"))
            
            counts <- ggplot(salmon_use, aes(Date, Count, col = get(Salmon_type))) + geom_point(size = 0.5) + geom_line() + xlab("") + ylab("Counts") + background_grid(major = "x", minor="x", size.minor = 0.5, colour.minor = "grey90") #+ theme(legend.position = "none")
            
            qpcrcorrQ <- ggplot(salmon_use, aes(Date, Qcorr_qPCR)) + geom_point(size = 0.5) + geom_line() + xlab("") + ylab("Qcfs*QUANT") + background_grid(major = "x", minor="x", size.minor = 0.5, colour.minor = "grey90")
            
            qpcr <- ggplot(salmon_use, aes(Date, exp(-CT_mean))) + geom_point(size = 0.5) + geom_line() + xlab("") + ylab("exp(qPCR)") + background_grid(major = "x", minor="x", size.minor = 0.5, colour.minor = "grey90")
            
            QUANT <- ggplot(salmon_use, aes(Date, QUANT_mean)) + geom_point(size = 0.5) + geom_line() + xlab("") + ylab("QUANT") + background_grid(major = "x", minor="x", size.minor = 0.5, colour.minor = "grey90")
            
            Q_cfs2 <- ggplot(salmon_use, aes(Date, Q_cfs2)) + geom_line() + xlab("") + ylab("Q_cfs2") + background_grid(major = "x", minor="x", size.minor = 0.5, colour.minor = "grey90")
            
            water <- ggplot(salmon_use, aes(Date, log(Gage_Height))) + geom_line() + xlab("") + ylab("Log(Gage Height)") + background_grid(major = "x", minor="x", size.minor = 0.5, colour.minor = "grey90")
            
            temp <- ggplot(salmon_use, aes(Date, Temp_C)) + geom_line() + xlab("") + ylab("Temp_C") + background_grid(major = "x", minor="x", size.minor = 0.5, colour.minor = "grey90")
            
            #   plot_grid(counts, qpcrcorrQ, qpcr, Q_cfs2, water, temp, nrow = 6, align = "v")
            plot_grid(counts, qpcrcorrQ, QUANT, Q_cfs2, temp, nrow = 5, align = "v")
            
}
#####################
#End: Timeline Function
#####################


#####################
#Begin: Graph Timelines
#####################

timelines_full("sockeye") # use sockeye, not Sockeye 

timelines_full("sockeye", "2015-06-18", "2015-08-01") # after looking at the full timeline, truncating to 01 Aug to prevent measuring dead-salmon DNA in 2016

timelines_full("sockeye", "2016-06-18", "2016-08-01") # after looking at the full timeline, truncating to 01 Aug to prevent measuring dead-salmon DNA in 2016


timelines_full("sockeye", "2015-04-15", "2015-06-10", "_juv") # downstream weir dates only for smolts
timelines_full("sockeye", "2016-04-15", "2016-06-10", "_juv") # downstream weir dates only for smolts

timelines_full("coho") # use coho, not Coho
timelines_full("coho", "2015-06-18", "2015-10-30") # after looking at the full timeline
timelines_full("coho", "2015-04-15", "2015-06-10", "_juv") # downstream weir dates only 
timelines_full("coho", "2016-06-18", "2016-10-30") # after looking at the full timeline
timelines_full("coho", "2016-04-15", "2016-06-10", "_juv") # downstream weir dates only 

#####################
#End: Graph Timelines
#####################


#####################
#Begin:  Models of Counts
#####################
datasubset = subset(sockeye_adult, Sockeyetype == "Sockeye_Adults" &!is.na(QUANT_mean) & Year %in% c("2015", "2016"))
datasubset <- datasubset %>% mutate(QUANTQ = QUANT_mean * Q_cfs2) # to allow plotting of residuals against QUANT_mean * Q_cfs
sockeye.adult.2016.reg <- glm(Count ~ Qcorr_qPCR , family = "quasipoisson", data = datasubset[datasubset$Year == "2016"  , ])
summary(sockeye.adult.2016.reg)
sockeye.adult.2015.reg <- glm(Count ~ Qcorr_qPCR , family = "quasipoisson", data = datasubset[datasubset$Year == "2015"  , ])
summary(sockeye.adult.2015.reg)

sockeye.adult.yearinteraction.reg <- glm(Count ~ Qcorr_qPCR*Year , family = "quasipoisson", data = datasubset)
summary(sockeye.adult.yearinteraction.reg)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           4.2905     0.3690  11.628  5.9e-16 ***
#   Qcorr_qPCR            3.8265     1.1521   3.321  0.00166 ** 
#   Year2016             -1.2439     0.5413  -2.298  0.02570 *  
#   Qcorr_qPCR:Year2016   0.5098     1.3669   0.373  0.71071  

datasubset.coho.adult = subset(coho_adult, Cohotype == "Coho_Adults" &!is.na(QUANT_mean) & Year %in% c("2015", "2016"))
coho.adult.2016.reg <- glm(Count ~ Qcorr_qPCR , family = "quasipoisson", data = datasubset.coho.adult[datasubset.coho.adult$Year == "2016"  , ])
summary(coho.adult.2016.reg)
coho.adult.2015.reg <- glm(Count ~ Qcorr_qPCR , family = "quasipoisson", data = datasubset.coho.adult[datasubset.coho.adult$Year == "2015"  , ])
summary(coho.adult.2015.reg)

coho.adult.yearinteraction.reg <- glm(Count ~ Qcorr_qPCR*Year , family = "quasipoisson", data = datasubset.coho.adult)
summary(coho.adult.yearinteraction.reg) # this should be different since one year had simple mix of life stages and the other was complex

datasubset.coho.total = subset(coho_adult, Cohotype == "Coho_Total" &!is.na(QUANT_mean) & Year %in% c("2015", "2016"))
coho.total.2016.reg <- glm(Count ~ Qcorr_qPCR , family = "quasipoisson", data = datasubset.coho.total[datasubset.coho.total$Year == "2016"  , ])
summary(coho.total.2016.reg)
coho.total.2015.reg <- glm(Count ~ Qcorr_qPCR , family = "quasipoisson", data = datasubset.coho.total[datasubset.coho.total$Year == "2015"  , ])
summary(coho.total.2015.reg)

coho.total.yearinteraction.reg <- glm(Count ~ Qcorr_qPCR*Year , family = "quasipoisson", data = datasubset.coho.total)
summary(coho.total.yearinteraction.reg) # this should be different since one year had simple mix of life stages and the other was complex
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          2.36228    0.16951  13.936  < 2e-16 ***
#   Qcorr_qPCR           3.61484    0.54193   6.670 1.58e-09 ***
#   Year2016             0.03916    0.25310   0.155    0.877    
# Qcorr_qPCR:Year2016  5.25906    2.20471   2.385    0.019 *  

#datasubset.coho.nojacks = subset(coho_adult, Cohotype == "Coho_nojacks" &!is.na(QUANT_mean) & Year %in% c("2015", "2016"))
#coho.nojacks.2016.reg <- glm(Count ~ Qcorr_qPCR , family = "quasipoisson", data = datasubset.coho.nojacks[datasubset.coho.nojacks$Year == "2016"  , ])
#summary(coho.adult.2016.reg)
#coho.nojacks.2015.reg <- glm(Count ~ Qcorr_qPCR , family = "quasipoisson", data = datasubset.coho.nojacks[datasubset.coho.nojacks$Year == "2015"  , ])
#summary(coho.adult.2015.reg)

datasubset.sockeye.juv = subset(sockeye_juv, Sockeyetype == "Sockeye_Smolt" &!is.na(Q_cfs2)& !is.na(QUANT_mean) & Year %in% c("2015", "2016"))
datasubset.coho.juv = subset(coho_juv, Cohotype == "Coho_Smolt" &!is.na(QUANT_mean) & Year %in% c("2015", "2016"))

sockeye.juv.2015.reg <- glm(Count ~ Qcorr_qPCR , family = "quasipoisson", data = datasubset.sockeye.juv[datasubset.sockeye.juv$Year == "2015"  , ])
sockeye.juv.2016.reg <- glm(Count ~ Qcorr_qPCR , family = "quasipoisson", data = datasubset.sockeye.juv[datasubset.sockeye.juv$Year == "2016"  , ])
coho.juv.2015.reg <- glm(Count ~ Qcorr_qPCR , family = "quasipoisson", data = datasubset.coho.juv[datasubset.coho.juv$Year == "2015"  , ])
coho.juv.2016.reg <- glm(Count ~ Qcorr_qPCR , family = "quasipoisson", data = datasubset.coho.juv[datasubset.coho.juv$Year == "2016"  , ])
summary(sockeye.juv.2015.reg)
summary(sockeye.juv.2016.reg)

sockeye.juv.yearinteraction.reg <- glm(Count ~ Qcorr_qPCR*Year , family = "quasipoisson", data = datasubset.sockeye.juv)
summary(sockeye.juv.yearinteraction.reg) 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           3.7486     0.7014   5.344 0.000103 ***
#   Qcorr_qPCR          570.2408   158.4968   3.598 0.002911 ** 
#   Year2016              0.6448     0.8343   0.773 0.452469    
# Qcorr_qPCR:Year2016 111.7506   206.6359   0.541 0.597139    
###########PLOT FITTED RELATIONSHIPS
par(mfrow=c(3,2),mar=c(3,4.5,1,1),oma=c(2,3,1,1))

plot(datasubset[datasubset$Year == "2015"  , ]$Qcorr_qPCR,datasubset[datasubset$Year == "2015"  , ]$Count,type="n",ylab="Sockeye Counts",main="2015",cex.lab=1.5,xlab="")
y_up=exp(augment(sockeye.adult.2015.reg)$.fitted[order(datasubset$Qcorr_qPCR[which(datasubset$Year==2015)])] + 2*augment(sockeye.adult.2015.reg)$.se.fit[order(datasubset$Qcorr_qPCR[which(datasubset$Year==2015)])])
y_down=exp(augment(sockeye.adult.2015.reg)$.fitted[order(datasubset$Qcorr_qPCR[which(datasubset$Year==2015)])] - 2*augment(sockeye.adult.2015.reg)$.se.fit[order(datasubset$Qcorr_qPCR[which(datasubset$Year==2015)])])
polygon(c(rev(datasubset$Qcorr_qPCR[which(datasubset$Year==2015)][order(datasubset$Qcorr_qPCR[which(datasubset$Year==2015)])]),datasubset$Qcorr_qPCR[which(datasubset$Year==2015)][order(datasubset$Qcorr_qPCR[which(datasubset$Year==2015)])]),(c(rev(y_up),y_down)),col="grey80",border=NA)
points(datasubset[datasubset$Year == "2015"  , ]$Qcorr_qPCR,datasubset[datasubset$Year == "2015"  , ]$Count, pch = 19)
lines(datasubset$Qcorr_qPCR[which(datasubset$Year==2015)][order(datasubset$Qcorr_qPCR[which(datasubset$Year==2015)])], exp(augment(sockeye.adult.2015.reg)$.fitted)[order(datasubset$Qcorr_qPCR[which(datasubset$Year==2015)])], lwd=2,lty=2, col="blue")


plot(datasubset[datasubset$Year == "2016"  , ]$Qcorr_qPCR,datasubset[datasubset$Year == "2016"  , ]$Count,main="2016",ylab = "",cex.lab=1.5,xlab="")
y_up=exp(augment(sockeye.adult.2016.reg)$.fitted[order(datasubset$Qcorr_qPCR[which(datasubset$Year==2016)])] + 2*augment(sockeye.adult.2016.reg)$.se.fit[order(datasubset$Qcorr_qPCR[which(datasubset$Year==2016)])])
y_down=exp(augment(sockeye.adult.2016.reg)$.fitted[order(datasubset$Qcorr_qPCR[which(datasubset$Year==2016)])] - 2*augment(sockeye.adult.2016.reg)$.se.fit[order(datasubset$Qcorr_qPCR[which(datasubset$Year==2016)])])
polygon(c(rev(datasubset$Qcorr_qPCR[which(datasubset$Year==2016)][order(datasubset$Qcorr_qPCR[which(datasubset$Year==2016)])]),datasubset$Qcorr_qPCR[which(datasubset$Year==2016)][order(datasubset$Qcorr_qPCR[which(datasubset$Year==2016)])]),(c(rev(y_up),y_down)),col="grey80",border=NA)
points(datasubset[datasubset$Year == "2016"  , ]$Qcorr_qPCR,datasubset[datasubset$Year == "2016"  , ]$Count, pch = 19)
lines(datasubset$Qcorr_qPCR[which(datasubset$Year==2016)][order(datasubset$Qcorr_qPCR[which(datasubset$Year==2016)])], exp(augment(sockeye.adult.2016.reg)$.fitted)[order(datasubset$Qcorr_qPCR[which(datasubset$Year==2016)])], lwd=2,lty=2, col="blue")


plot(datasubset.coho.total[datasubset.coho.total$Year == "2015"  , ]$Qcorr_qPCR,datasubset.coho.total[datasubset.coho.total$Year == "2015"  , ]$Count, ylab = "Coho Counts",cex.lab=1.5,xlab="")
y_up=exp(augment(coho.total.2015.reg)$.fitted[order(datasubset.coho.total$Qcorr_qPCR[which(datasubset.coho.total$Year==2015)])] + 2*augment(coho.total.2015.reg)$.se.fit[order(datasubset.coho.total$Qcorr_qPCR[which(datasubset.coho.total$Year==2015)])])
y_down=exp(augment(coho.total.2015.reg)$.fitted[order(datasubset.coho.total$Qcorr_qPCR[which(datasubset.coho.total$Year==2015)])] - 2*augment(coho.total.2015.reg)$.se.fit[order(datasubset.coho.total$Qcorr_qPCR[which(datasubset.coho.total$Year==2015)])])
polygon(c(rev(datasubset.coho.total$Qcorr_qPCR[which(datasubset.coho.total$Year==2015)][order(datasubset.coho.total$Qcorr_qPCR[which(datasubset.coho.total$Year==2015)])]),datasubset.coho.total$Qcorr_qPCR[which(datasubset.coho.total$Year==2015)][order(datasubset.coho.total$Qcorr_qPCR[which(datasubset.coho.total$Year==2015)])]),(c(rev(y_up),y_down)),col="grey80",border=NA)
points(datasubset.coho.total[datasubset.coho.total$Year == "2015"  , ]$Qcorr_qPCR,datasubset.coho.total[datasubset.coho.total$Year == "2015"  , ]$Count, pch = 19)
lines(datasubset.coho.total$Qcorr_qPCR[which(datasubset.coho.total$Year==2015)][order(datasubset.coho.total$Qcorr_qPCR[which(datasubset.coho.total$Year==2015)])], exp(augment(coho.total.2015.reg)$.fitted)[order(datasubset.coho.total$Qcorr_qPCR[which(datasubset.coho.total$Year==2015)])], lwd=2,lty=2, col="blue")

plot(datasubset.coho.total[datasubset.coho.total$Year == "2016"  , ]$Qcorr_qPCR,datasubset.coho.total[datasubset.coho.total$Year == "2016"  , ]$Count,ylab = "",cex.lab=1.5,xlab="")
y_up=exp(augment(coho.total.2016.reg)$.fitted[order(datasubset.coho.total$Qcorr_qPCR[which(datasubset.coho.total$Year==2016)])] + 2*augment(coho.total.2016.reg)$.se.fit[order(datasubset.coho.total$Qcorr_qPCR[which(datasubset.coho.total$Year==2016)])])
y_down=exp(augment(coho.total.2016.reg)$.fitted[order(datasubset.coho.total$Qcorr_qPCR[which(datasubset.coho.total$Year==2016)])] - 2*augment(coho.total.2016.reg)$.se.fit[order(datasubset.coho.total$Qcorr_qPCR[which(datasubset.coho.total$Year==2016)])])
polygon(c(rev(datasubset.coho.total$Qcorr_qPCR[which(datasubset.coho.total$Year==2016)][order(datasubset.coho.total$Qcorr_qPCR[which(datasubset.coho.total$Year==2016)])]),datasubset.coho.total$Qcorr_qPCR[which(datasubset.coho.total$Year==2016)][order(datasubset.coho.total$Qcorr_qPCR[which(datasubset.coho.total$Year==2016)])]),(c(rev(y_up),y_down)),col="grey80",border=NA)
points(datasubset.coho.total[datasubset.coho.total$Year == "2016"  , ]$Qcorr_qPCR,datasubset.coho.total[datasubset.coho.total$Year == "2016"  , ]$Count, pch = 19)
lines(datasubset.coho.total$Qcorr_qPCR[which(datasubset.coho.total$Year==2016)][order(datasubset.coho.total$Qcorr_qPCR[which(datasubset.coho.total$Year==2016)])], exp(augment(coho.total.2016.reg)$.fitted)[order(datasubset.coho.total$Qcorr_qPCR[which(datasubset.coho.total$Year==2016)])], lwd=2,lty=2, col="blue")

plot(datasubset.sockeye.juv[datasubset.sockeye.juv$Year == "2015"  , ]$Qcorr_qPCR,datasubset.sockeye.juv[datasubset.sockeye.juv$Year == "2015"  , ]$Count,type="n",ylab="Sockeye Smolts",cex.lab=1.5,xlab="")
y_up=exp(augment(sockeye.juv.2015.reg)$.fitted[order(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2015)])] + 2*augment(sockeye.juv.2015.reg)$.se.fit[order(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2015)])])
y_down=exp(augment(sockeye.juv.2015.reg)$.fitted[order(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2015)])] - 2*augment(sockeye.juv.2015.reg)$.se.fit[order(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2015)])])
polygon(c(rev(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2015)][order(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2015)])]),datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2015)][order(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2015)])]),(c(rev(y_up),y_down)),col="grey80",border=NA)
points(datasubset.sockeye.juv[datasubset.sockeye.juv$Year == "2015"  , ]$Qcorr_qPCR,datasubset.sockeye.juv[datasubset.sockeye.juv$Year == "2015"  , ]$Count, pch = 19)
lines(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2015)][order(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2015)])], exp(augment(sockeye.juv.2015.reg)$.fitted)[order(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2015)])], lwd=2,lty=2, col="blue")


plot(datasubset.sockeye.juv[datasubset.sockeye.juv$Year == "2016"  , ]$Qcorr_qPCR,datasubset.sockeye.juv[datasubset.sockeye.juv$Year == "2016"  , ]$Count,ylab = "",cex.lab=1.5,xlab="")
y_up=exp(augment(sockeye.juv.2016.reg)$.fitted[order(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2016)])] + 2*augment(sockeye.juv.2016.reg)$.se.fit[order(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2016)])])
y_down=exp(augment(sockeye.juv.2016.reg)$.fitted[order(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2016)])] - 2*augment(sockeye.juv.2016.reg)$.se.fit[order(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2016)])])
polygon(c(rev(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2016)][order(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2016)])]),datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2016)][order(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2016)])]),(c(rev(y_up),y_down)),col="grey80",border=NA)
points(datasubset.sockeye.juv[datasubset.sockeye.juv$Year == "2016"  , ]$Qcorr_qPCR,datasubset.sockeye.juv[datasubset.sockeye.juv$Year == "2016"  , ]$Count, pch = 19)
lines(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2016)][order(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2016)])], exp(augment(sockeye.juv.2016.reg)$.fitted)[order(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2016)])], lwd=2,lty=2, col="blue")


mtext("Flow x DNA",side=1,cex=1.5,outer=T)

# 
# 
# plot(datasubset.coho.adult[datasubset.coho.adult$Year == "2015"  , ]$Qcorr_qPCR,datasubset.coho.adult[datasubset.coho.adult$Year == "2015"  , ]$Count)
# y_up=exp(augment(coho.adult.2015.reg)$.fitted[order(datasubset.coho.adult$Qcorr_qPCR[which(datasubset.coho.adult$Year==2015)])] + 2*augment(coho.adult.2015.reg)$.se.fit[order(datasubset.coho.adult$Qcorr_qPCR[which(datasubset.coho.adult$Year==2015)])])
# y_down=exp(augment(coho.adult.2015.reg)$.fitted[order(datasubset.coho.adult$Qcorr_qPCR[which(datasubset.coho.adult$Year==2015)])] - 2*augment(coho.adult.2015.reg)$.se.fit[order(datasubset.coho.adult$Qcorr_qPCR[which(datasubset.coho.adult$Year==2015)])])
# polygon(c(rev(datasubset.coho.adult$Qcorr_qPCR[which(datasubset.coho.adult$Year==2015)][order(datasubset.coho.adult$Qcorr_qPCR[which(datasubset.coho.adult$Year==2015)])]),datasubset.coho.adult$Qcorr_qPCR[which(datasubset.coho.adult$Year==2015)][order(datasubset.coho.adult$Qcorr_qPCR[which(datasubset.coho.adult$Year==2015)])]),(c(rev(y_up),y_down)),col="grey80",border=NA)
# points(datasubset.coho.adult[datasubset.coho.adult$Year == "2015"  , ]$Qcorr_qPCR,datasubset.coho.adult[datasubset.coho.adult$Year == "2015"  , ]$Count, pch = 19)
# lines(datasubset.coho.adult$Qcorr_qPCR[which(datasubset.coho.adult$Year==2015)][order(datasubset.coho.adult$Qcorr_qPCR[which(datasubset.coho.adult$Year==2015)])], exp(augment(coho.adult.2015.reg)$.fitted)[order(datasubset.coho.adult$Qcorr_qPCR[which(datasubset.coho.adult$Year==2015)])], lwd=2,lty=2, col="blue")
# 
# 
# plot(datasubset.coho.adult[datasubset.coho.adult$Year == "2016"  , ]$Qcorr_qPCR,datasubset.coho.adult[datasubset.coho.adult$Year == "2016"  , ]$Count)
# y_up=exp(augment(coho.adult.2016.reg)$.fitted[order(datasubset.coho.adult$Qcorr_qPCR[which(datasubset.coho.adult$Year==2016)])] + 2*augment(coho.adult.2016.reg)$.se.fit[order(datasubset.coho.adult$Qcorr_qPCR[which(datasubset.coho.adult$Year==2016)])])
# y_down=exp(augment(coho.adult.2016.reg)$.fitted[order(datasubset.coho.adult$Qcorr_qPCR[which(datasubset.coho.adult$Year==2016)])] - 2*augment(coho.adult.2016.reg)$.se.fit[order(datasubset.coho.adult$Qcorr_qPCR[which(datasubset.coho.adult$Year==2016)])])
# polygon(c(rev(datasubset.coho.adult$Qcorr_qPCR[which(datasubset.coho.adult$Year==2016)][order(datasubset.coho.adult$Qcorr_qPCR[which(datasubset.coho.adult$Year==2016)])]),datasubset.coho.adult$Qcorr_qPCR[which(datasubset.coho.adult$Year==2016)][order(datasubset.coho.adult$Qcorr_qPCR[which(datasubset.coho.adult$Year==2016)])]),(c(rev(y_up),y_down)),col="grey80",border=NA)
# points(datasubset.coho.adult[datasubset.coho.adult$Year == "2016"  , ]$Qcorr_qPCR,datasubset.coho.adult[datasubset.coho.adult$Year == "2016"  , ]$Count, pch = 19)
# lines(datasubset.coho.adult$Qcorr_qPCR[which(datasubset.coho.adult$Year==2016)][order(datasubset.coho.adult$Qcorr_qPCR[which(datasubset.coho.adult$Year==2016)])], exp(augment(coho.adult.2016.reg)$.fitted)[order(datasubset.coho.adult$Qcorr_qPCR[which(datasubset.coho.adult$Year==2016)])], lwd=2,lty=2, col="blue")

###########
#
#Sockeye Adult
#
###########
par(mfrow=c(3,2),mar=c(3,4.5,1,1),oma=c(2,3,1,1))

plot(datasubset$Date[which(datasubset$Year==2015)], datasubset$Count[which(datasubset$Year==2015)], lwd=2, type='n',pch=19,cex=0.25,ylim=c(0,900),ylab="Sockeye Count",xlab="",cex.lab=1.5, main = "2015")
y_up=exp(augment(sockeye.adult.2015.reg)$.fitted + 2*augment(sockeye.adult.2015.reg)$.se.fit)
y_down=exp(augment(sockeye.adult.2015.reg)$.fitted - 2*augment(sockeye.adult.2015.reg)$.se.fit)
polygon(c(rev(datasubset$Date[which(datasubset$Year==2015)]),datasubset$Date[which(datasubset$Year==2015)]),(c(rev(y_up),y_down)),col="grey80",border=NA)
points(datasubset$Date[which(datasubset$Year==2015)], datasubset$Count[which(datasubset$Year==2015)], lwd=2, type='b',pch=19,cex=0.25)
lines(datasubset$Date[which(datasubset$Year==2015)], exp(augment(sockeye.adult.2015.reg)$.fitted), lwd=2,lty=2, col="blue")
# mtext("Sockeye Total Counts Quasipoisson model",side=2,cex=2,outer=T)
#text(datasubset$Date[which(datasubset$Year==2015)][35],800,"2015",cex=1.5)


plot(datasubset$Date[which(datasubset$Year==2016)], datasubset$Count[which(datasubset$Year==2016)], lwd=2, type='n',pch=19,cex=0.25,ylim=c(0,900),ylab="",xlab="",cex.lab=1.5, main = "2016")
y_up=exp(augment(sockeye.adult.2016.reg)$.fitted + 2*augment(sockeye.adult.2016.reg)$.se.fit)
y_down=exp(augment(sockeye.adult.2016.reg)$.fitted - 2*augment(sockeye.adult.2016.reg)$.se.fit)
polygon(c(rev(datasubset$Date[which(datasubset$Year==2016)]),datasubset$Date[which(datasubset$Year==2016)]),(c(rev(y_up),y_down)),col="grey80",border=NA)
points(datasubset$Date[which(datasubset$Year==2016)], datasubset$Count[which(datasubset$Year==2016)], lwd=2, type='b',pch=19,cex=0.25)
lines(datasubset$Date[which(datasubset$Year==2016)], exp(augment(sockeye.adult.2016.reg)$.fitted), lwd=2,lty=2, col="blue")
# mtext("Sockeye Total Counts Quasipoisson model",side=2,cex=2,outer=T)
#text(datasubset$Date[which(datasubset$Year==2016)][35],800,"2016",cex=1.5)

###########
#
#Coho Total
#
###########

plot(datasubset.coho.total$Date[which(datasubset.coho.total$Year==2015)], datasubset.coho.total$Count[which(datasubset.coho.total$Year==2015)], lwd=2, type='n',pch=19,cex=0.25,ylim=c(0,200),ylab="Coho Count",xlab="",cex.lab=1.5)
y_up=exp(augment(coho.total.2015.reg)$.fitted + 2*augment(coho.total.2015.reg)$.se.fit)
y_down=exp(augment(coho.total.2015.reg)$.fitted - 2*augment(coho.total.2015.reg)$.se.fit)
polygon(c(rev(datasubset.coho.total$Date[which(datasubset.coho.total$Year==2015)]),datasubset.coho.total$Date[which(datasubset.coho.total$Year==2015)]),(c(rev(y_up),y_down)),col="grey80",border=NA)
points(datasubset.coho.total$Date[which(datasubset.coho.total$Year==2015)], datasubset.coho.total$Count[which(datasubset.coho.total$Year==2015)], lwd=2, type='b',pch=19,cex=0.25)
lines(datasubset.coho.total$Date[which(datasubset.coho.total$Year==2015)], exp(augment(coho.total.2015.reg)$.fitted), lwd=2,lty=2, col="blue")
#mtext("coho Total Counts Quasipoisson model",side=2,cex=2,outer=T)
#text(datasubset.coho.total$Date[which(datasubset.coho.total$Year==2015)][35],800,"2015",cex=1.5)


plot(datasubset.coho.total$Date[which(datasubset.coho.total$Year==2016)], datasubset.coho.total$Count[which(datasubset.coho.total$Year==2016)], lwd=2, type='n',pch=19,cex=0.25,ylim=c(0,200),ylab="",xlab="",cex.lab=1.5)
y_up=exp(augment(coho.total.2016.reg)$.fitted + 2*augment(coho.total.2016.reg)$.se.fit)
y_down=exp(augment(coho.total.2016.reg)$.fitted - 2*augment(coho.total.2016.reg)$.se.fit)
polygon(c(rev(datasubset.coho.total$Date[which(datasubset.coho.total$Year==2016)]),datasubset.coho.total$Date[which(datasubset.coho.total$Year==2016)]),(c(rev(y_up),y_down)),col="grey80",border=NA)
points(datasubset.coho.total$Date[which(datasubset.coho.total$Year==2016)], datasubset.coho.total$Count[which(datasubset.coho.total$Year==2016)], lwd=2, type='b',pch=19,cex=0.25)
lines(datasubset.coho.total$Date[which(datasubset.coho.total$Year==2016)], exp(augment(coho.total.2016.reg)$.fitted), lwd=2,lty=2, col="blue")
#mtext("coho Total Counts Quasipoisson model",side=2,cex=2,outer=T)
#text(datasubset.coho.total$Date[which(datasubset.coho.total$Year==2016)][35],800,"2016",cex=1.5)

###################
###########
#
#Sockeye Smolt
#
###########
plot(datasubset.sockeye.juv$Date[which(datasubset.sockeye.juv$Year==2015)], datasubset.sockeye.juv$Count[which(datasubset.sockeye.juv$Year==2015)], lwd=2, type='n',pch=19,cex=0.25,ylim=c(0,2000),ylab="Sockeye Smolt",xlab="",cex.lab=1.5)
y_up=exp(augment(sockeye.juv.2015.reg)$.fitted + 2*augment(sockeye.juv.2015.reg)$.se.fit)
y_down=exp(augment(sockeye.juv.2015.reg)$.fitted - 2*augment(sockeye.juv.2015.reg)$.se.fit)
polygon(c(rev(datasubset.sockeye.juv$Date[which(datasubset.sockeye.juv$Year==2015)]),datasubset.sockeye.juv$Date[which(datasubset.sockeye.juv$Year==2015)]),(c(rev(y_up),y_down)),col="grey80",border=NA)
points(datasubset.sockeye.juv$Date[which(datasubset.sockeye.juv$Year==2015)], datasubset.sockeye.juv$Count[which(datasubset.sockeye.juv$Year==2015)], lwd=2, type='b',pch=19,cex=0.25)
lines(datasubset.sockeye.juv$Date[which(datasubset.sockeye.juv$Year==2015)], exp(augment(sockeye.juv.2015.reg)$.fitted), lwd=2,lty=2, col="blue")
# mtext("Sockeye Total Counts Quasipoisson model",side=2,cex=2,outer=T)
#text(datasubset.sockeye.juv$Date[which(datasubset.sockeye.juv$Year==2015)][35],800,"2015",cex=1.5)


plot(datasubset.sockeye.juv$Date[which(datasubset.sockeye.juv$Year==2016)], datasubset.sockeye.juv$Count[which(datasubset.sockeye.juv$Year==2016)], lwd=2, type='n',pch=19,cex=0.25,ylim=c(0,2000),ylab="",xlab="",cex.lab=1.5)
y_up=exp(augment(sockeye.juv.2016.reg)$.fitted + 2*augment(sockeye.juv.2016.reg)$.se.fit)
y_down=exp(augment(sockeye.juv.2016.reg)$.fitted - 2*augment(sockeye.juv.2016.reg)$.se.fit)
polygon(c(rev(datasubset.sockeye.juv$Date[which(datasubset.sockeye.juv$Year==2016)]),datasubset.sockeye.juv$Date[which(datasubset.sockeye.juv$Year==2016)]),(c(rev(y_up),y_down)),col="grey80",border=NA)
points(datasubset.sockeye.juv$Date[which(datasubset.sockeye.juv$Year==2016)], datasubset.sockeye.juv$Count[which(datasubset.sockeye.juv$Year==2016)], lwd=2, type='b',pch=19,cex=0.25)
lines(datasubset.sockeye.juv$Date[which(datasubset.sockeye.juv$Year==2016)], exp(augment(sockeye.juv.2016.reg)$.fitted), lwd=2,lty=2, col="blue")
# mtext("Sockeye Total Counts Quasipoisson model",side=2,cex=2,outer=T)
#text(datasubset.sockeye.juv$Date[which(datasubset.sockeye.juv$Year==2016)][35],800,"2016",cex=1.5)
# 
# 
# ###########
# #
# #Coho Adult
# #
# ###########
# plot(datasubset.coho.adult$Date[which(datasubset.coho.adult$Year==2015)], datasubset.coho.adult$Count[which(datasubset.coho.adult$Year==2015)], lwd=2, type='n',pch=19,cex=0.25,ylim=c(0,150),ylab="Coho Count",xlab="",cex.lab=1.5)
# y_up=exp(augment(coho.adult.2015.reg)$.fitted + 2*augment(coho.adult.2015.reg)$.se.fit)
# y_down=exp(augment(coho.adult.2015.reg)$.fitted - 2*augment(coho.adult.2015.reg)$.se.fit)
# polygon(c(rev(datasubset.coho.adult$Date[which(datasubset.coho.adult$Year==2015)]),datasubset.coho.adult$Date[which(datasubset.coho.adult$Year==2015)]),(c(rev(y_up),y_down)),col="grey80",border=NA)
# points(datasubset.coho.adult$Date[which(datasubset.coho.adult$Year==2015)], datasubset.coho.adult$Count[which(datasubset.coho.adult$Year==2015)], lwd=2, type='b',pch=19,cex=0.25)
# lines(datasubset.coho.adult$Date[which(datasubset.coho.adult$Year==2015)], exp(augment(coho.adult.2015.reg)$.fitted), lwd=2,lty=2, col="blue")
# # mtext("coho Total Counts Quasipoisson model",side=2,cex=2,outer=T)
# text(datasubset.coho.adult$Date[which(datasubset.coho.adult$Year==2015)][35],800,"2015",cex=1.5)
# 
# 
# plot(datasubset.coho.adult$Date[which(datasubset.coho.adult$Year==2016)], datasubset.coho.adult$Count[which(datasubset.coho.adult$Year==2016)], lwd=2, type='n',pch=19,cex=0.25,ylim=c(0,100),ylab="",xlab="",cex.lab=1.5)
# y_up=exp(augment(coho.adult.2016.reg)$.fitted + 2*augment(coho.adult.2016.reg)$.se.fit)
# y_down=exp(augment(coho.adult.2016.reg)$.fitted - 2*augment(coho.adult.2016.reg)$.se.fit)
# polygon(c(rev(datasubset.coho.adult$Date[which(datasubset.coho.adult$Year==2016)]),datasubset.coho.adult$Date[which(datasubset.coho.adult$Year==2016)]),(c(rev(y_up),y_down)),col="grey80",border=NA)
# points(datasubset.coho.adult$Date[which(datasubset.coho.adult$Year==2016)], datasubset.coho.adult$Count[which(datasubset.coho.adult$Year==2016)], lwd=2, type='b',pch=19,cex=0.25)
# lines(datasubset.coho.adult$Date[which(datasubset.coho.adult$Year==2016)], exp(augment(coho.adult.2016.reg)$.fitted), lwd=2,lty=2, col="blue")
# # mtext("coho Total Counts Quasipoisson model",side=2,cex=2,outer=T)
# text(datasubset.coho.adult$Date[which(datasubset.coho.adult$Year==2016)][35],800,"2016",cex=1.5)
# 

###########
#
#Coho nojack
#
# ###########
# 
# plot(datasubset.coho.nojacks$Date[which(datasubset.coho.nojacks$Year==2015)], datasubset.coho.nojacks$Count[which(datasubset.coho.nojacks$Year==2015)], lwd=2, type='n',pch=19,cex=0.25,ylim=c(0,200),ylab="",xlab="",cex.lab=1.5)
# y_up=exp(augment(coho.nojacks.2015.reg)$.fitted + 2*augment(coho.nojacks.2015.reg)$.se.fit)
# y_down=exp(augment(coho.nojacks.2015.reg)$.fitted - 2*augment(coho.nojacks.2015.reg)$.se.fit)
# polygon(c(rev(datasubset.coho.nojacks$Date[which(datasubset.coho.nojacks$Year==2015)]),datasubset.coho.nojacks$Date[which(datasubset.coho.nojacks$Year==2015)]),(c(rev(y_up),y_down)),col="grey80",border=NA)
# points(datasubset.coho.nojacks$Date[which(datasubset.coho.nojacks$Year==2015)], datasubset.coho.nojacks$Count[which(datasubset.coho.nojacks$Year==2015)], lwd=2, type='b',pch=19,cex=0.25)
# lines(datasubset.coho.nojacks$Date[which(datasubset.coho.nojacks$Year==2015)], exp(augment(coho.nojacks.2015.reg)$.fitted), lwd=2,lty=2, col="blue")
# mtext("coho nojacks Counts Quasipoisson model",side=2,cex=2,outer=T)
# text(datasubset.coho.nojacks$Date[which(datasubset.coho.nojacks$Year==2015)][35],800,"2015",cex=1.5)
# 
# 
# plot(datasubset.coho.nojacks$Date[which(datasubset.coho.nojacks$Year==2016)], datasubset.coho.nojacks$Count[which(datasubset.coho.nojacks$Year==2016)], lwd=2, type='n',pch=19,cex=0.25,ylim=c(0,200),ylab="",xlab="",cex.lab=1.5)
# y_up=exp(augment(coho.nojacks.2016.reg)$.fitted + 2*augment(coho.nojacks.2016.reg)$.se.fit)
# y_down=exp(augment(coho.nojacks.2016.reg)$.fitted - 2*augment(coho.nojacks.2016.reg)$.se.fit)
# polygon(c(rev(datasubset.coho.nojacks$Date[which(datasubset.coho.nojacks$Year==2016)]),datasubset.coho.nojacks$Date[which(datasubset.coho.nojacks$Year==2016)]),(c(rev(y_up),y_down)),col="grey80",border=NA)
# points(datasubset.coho.nojacks$Date[which(datasubset.coho.nojacks$Year==2016)], datasubset.coho.nojacks$Count[which(datasubset.coho.nojacks$Year==2016)], lwd=2, type='b',pch=19,cex=0.25)
# lines(datasubset.coho.nojacks$Date[which(datasubset.coho.nojacks$Year==2016)], exp(augment(coho.nojacks.2016.reg)$.fitted), lwd=2,lty=2, col="blue")
# mtext("coho nojacks Counts Quasipoisson model",side=2,cex=2,outer=T)
# text(datasubset.coho.nojacks$Date[which(datasubset.coho.nojacks$Year==2016)][35],800,"2016",cex=1.5)
# 
# # sockeye.adult.2016.reg <- glm(Count ~ QUANT_mean:Q_cfs2 + Qcorr_qpcr.lead+ Temp_C, family = "quasipoisson", data = datasubset[datasubset$Year == "2016" & !is.na(datasubset$Qcorr_qpcr.lead) , ])
# # summary(sockeye.adult.2016.reg)
# # edna.sockeye.adult.2016.reg <- glm(Qcorr_qPCR ~ -1 + Count , data = datasubset[datasubset$Year == "2016"  , ])
# # summary(edna.sockeye.adult.2016.reg)
# 


#####################
#Begin:  Models of eDNA
#####################






#####################
#End:  Models of eDNA
#####################

#############################
#######################################################################################
#############################End of super simple version. All code below is from previous versions
####################################################################################################################
#############################
#######################################################################################
#############################

# plot(datasubset$Date[which(datasubset$Year==2016)],datasubset$Qcorr_qPCR[which(datasubset$Year==2016)], type='b',pch=1)
# points(datasubset$Date[which(datasubset$Year==2016)],(predict(edna.sockeye.adult.2016.reg)), lty=2, col="black",pch=19)
# text(datasubset$Date[which(datasubset$Year==2016)][35],0.6,"2016",cex=1.5) 
# plot((augment(edna.sockeye.adult.2016.reg)$.fitted) ~ augment(edna.sockeye.adult.2016.reg)$Qcorr_qPCR)


# plot(datasubset$Date[which(datasubset$Year==2016)], datasubset$Count[which(datasubset$Year==2016)], type='b', ylim=c(0,900),pch=19)
# points(datasubset$Date[which(datasubset$Year==2016& !is.na(datasubset$Qcorr_qpcr.lead))], exp(predict(sockeye.adult.2016.reg)), lty=2, col="black")
# text(datasubset$Date[which(datasubset$Year==2016& !is.na(datasubset$Qcorr_qpcr.lead))][35],550,"2016",cex=1.5) 




