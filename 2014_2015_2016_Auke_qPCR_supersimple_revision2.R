library(broom)
library(cowplot)
library(MASS)
library(tidyverse)

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

Coholong_all <- Coholong_all %>% mutate(Q_cfs2 = a*(Gage_Height - b)^c)

# calculate a new Qcorr_qPCR value from QUANT_mean * Q_cfs2 # because old Qcorr_qPCR values used different translations of gauge height to Qcfs in different years
# this only matters for plotting since the models use QUANT_mean * Q_cfs2 directly
Sockeyelong_all <- Sockeyelong_all %>% mutate(Qcorr_qPCR = QUANT_mean * Q_cfs2) 

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
#set the only zero in the dataset to the minimum value

#####################
#Begin:  Models of Counts
#####################

#####DATA######
#Sockeye data
datasubset = subset(sockeye_adult, Sockeyetype == "Sockeye_Adults" &!is.na(QUANT_mean) & Year %in% c("2015", "2016"))
#to facilitate use of log(Qcorr_qPCR), set a single zero to the smallest observed Qcorr_qPCR - 1.002627e-05
datasubset$Qcorr_qPCR[2]=min(datasubset[datasubset$Qcorr_qPCR >0  , ]$Qcorr_qPCR)
datasubset.sockeye.juv = subset(sockeye_juv, Sockeyetype == "Sockeye_Smolt" &!is.na(Q_cfs2)& !is.na(QUANT_mean) & Year %in% c("2015", "2016"))

#Coho data
datasubset.coho.adult = subset(coho_adult, Cohotype == "Coho_Adults" &!is.na(QUANT_mean) & Year %in% c("2015", "2016"))
datasubset.coho.total = subset(coho_adult, Cohotype == "Coho_Total" &!is.na(QUANT_mean) & Year %in% c("2015", "2016"))
datasubset.coho.jack = subset(coho_adult, Cohotype == "Coho_Jack" &!is.na(QUANT_mean) & Year %in% c("2015", "2016"))
datasubset.coho.juv = subset(coho_adult, Cohotype == "Coho_Juv" &!is.na(QUANT_mean) & Year %in% c("2015", "2016"))

#datasubset.coho.trans is the same as datasubset.coho.total, but we remove one high leverage outlier from coho in 2016
datasubset.coho.trans=datasubset.coho.adult
datasubset.coho.trans <- datasubset.coho.trans %>% mutate(Count = (datasubset.coho.adult$Count+datasubset.coho.juv$Count+datasubset.coho.jack$Count))
#which(datasubset.coho.jack$Count==62) #exclude outlier 69th observation
datasubset.coho.trans <- datasubset.coho.trans %>% slice(-69)

#####MODELS######

#Sockeye
sockeye.adult.2015.reg <- glm(Count ~ log(Qcorr_qPCR) , family = "quasipoisson", data = datasubset[datasubset$Year == "2015"  , ])
sockeye.adult.2016.reg <- glm(Count ~ log(Qcorr_qPCR) , family = "quasipoisson", data = datasubset[datasubset$Year == "2016"  , ])
sockeye.adult.reg <- glm(Count ~ log(Qcorr_qPCR) , family = "quasipoisson", data = datasubset)

sockeye.juv.2015.reg <- glm(Count ~  log(Qcorr_qPCR) , family = "quasipoisson", data = datasubset.sockeye.juv[datasubset.sockeye.juv$Year == "2015"  , ])
sockeye.juv.2016.reg <- glm(Count ~  log(Qcorr_qPCR) , family = "quasipoisson", data = datasubset.sockeye.juv[datasubset.sockeye.juv$Year == "2016"  , ])
sockeye.juv.reg <- glm(Count ~ log(Qcorr_qPCR) , family = "quasipoisson", data = datasubset.sockeye.juv)

sockeye.adult.reg.int <- glm(Count ~ Year*log(Qcorr_qPCR), family = "quasipoisson", data = datasubset)
sockeye.juv.reg.int <- glm(Count ~ Year*log(Qcorr_qPCR) , family = "quasipoisson", data = datasubset.sockeye.juv)

#Coho
coho.trans.2015.reg<- glm(Count ~ log(Qcorr_qPCR) , family = "quasipoisson", data = datasubset.coho.trans[datasubset.coho.trans$Year == "2015"  , ])
coho.trans.2016.reg<- glm(Count ~ log(Qcorr_qPCR) , family = "quasipoisson", data = datasubset.coho.trans[datasubset.coho.trans$Year == "2016"  , ])
coho.trans.reg <- glm(Count ~ log(Qcorr_qPCR) , family = "quasipoisson", data = datasubset.coho.trans)
coho.trans.reg.int <- glm(Count ~ Year*log(Qcorr_qPCR) , family = "quasipoisson", data = datasubset.coho.trans)


#############################################
#   Models of eDNA vs. coho life stage      #
#############################################
datasubset.coho.trans2=datasubset.coho.adult
datasubset.coho.trans2 <- datasubset.coho.trans2 %>% mutate(CountA = datasubset.coho.adult$Count)%>% mutate(CountJack = datasubset.coho.jack$Count)%>% mutate(CountJuv = datasubset.coho.juv$Count)
#remove same outlier as before (with 62 Jacks)
coho.trans.2016.rev <- lm(Qcorr_qPCR ~ CountA + CountJuv + CountJack, data = datasubset.coho.trans2[datasubset.coho.trans2$Year == "2016" & datasubset.coho.trans2$CountJack <62  , ])
coho.trans.2015.rev <- lm(Qcorr_qPCR ~ CountA  + CountJuv + CountJack , data = datasubset.coho.trans2[datasubset.coho.trans2$Year == "2015"  , ])
coho.trans.rev.int <- lm(Qcorr_qPCR ~ Year*CountA  + Year*CountJuv + Year*CountJack , data = datasubset.coho.trans2[datasubset.coho.trans2$CountJack <62  , ])
coho.trans.rev <- lm(Qcorr_qPCR ~ CountA  + CountJuv + CountJack , data = datasubset.coho.trans2[datasubset.coho.trans2$CountJack <62  , ])



###########PLOT FITTED RELATIONSHIPS

par(mfrow=c(3,3),mar=c(3,4.5,1,0),oma=c(2,3,1,1))
#sockeye adult
plot(datasubset[datasubset$Year == "2015"  , ]$Qcorr_qPCR,datasubset[datasubset$Year == "2015"  , ]$Count,type="n",ylab="Sockeye Counts",main="2015",cex.lab=1.5,xlab="",ylim=c(0,600))
y_up=exp(augment(sockeye.adult.2015.reg)$.fitted[order(datasubset$Qcorr_qPCR[which(datasubset$Year==2015)])] + 2*augment(sockeye.adult.2015.reg)$.se.fit[order(datasubset$Qcorr_qPCR[which(datasubset$Year==2015)])])
y_down=exp(augment(sockeye.adult.2015.reg)$.fitted[order(datasubset$Qcorr_qPCR[which(datasubset$Year==2015)])] - 2*augment(sockeye.adult.2015.reg)$.se.fit[order(datasubset$Qcorr_qPCR[which(datasubset$Year==2015)])])
polygon(c(rev(datasubset$Qcorr_qPCR[which(datasubset$Year==2015)][order(datasubset$Qcorr_qPCR[which(datasubset$Year==2015)])]),datasubset$Qcorr_qPCR[which(datasubset$Year==2015)][order(datasubset$Qcorr_qPCR[which(datasubset$Year==2015)])]),(c(rev(y_up),y_down)),col="grey80",border=NA)
points(datasubset[datasubset$Year == "2015"  , ]$Qcorr_qPCR,datasubset[datasubset$Year == "2015"  , ]$Count, pch = 19)
lines(datasubset$Qcorr_qPCR[which(datasubset$Year==2015)][order(datasubset$Qcorr_qPCR[which(datasubset$Year==2015)])], exp(augment(sockeye.adult.2015.reg)$.fitted)[order(datasubset$Qcorr_qPCR[which(datasubset$Year==2015)])], lwd=2,lty=2, col="blue")

plot(datasubset[datasubset$Year == "2016"  , ]$Qcorr_qPCR,datasubset[datasubset$Year == "2016"  , ]$Count,main="2016",ylab = "",cex.lab=1.5,xlab="",ylim=c(0,600))
y_up=exp(augment(sockeye.adult.2016.reg)$.fitted[order(datasubset$Qcorr_qPCR[which(datasubset$Year==2016)])] + 2*augment(sockeye.adult.2016.reg)$.se.fit[order(datasubset$Qcorr_qPCR[which(datasubset$Year==2016)])])
y_down=exp(augment(sockeye.adult.2016.reg)$.fitted[order(datasubset$Qcorr_qPCR[which(datasubset$Year==2016)])] - 2*augment(sockeye.adult.2016.reg)$.se.fit[order(datasubset$Qcorr_qPCR[which(datasubset$Year==2016)])])
polygon(c(rev(datasubset$Qcorr_qPCR[which(datasubset$Year==2016)][order(datasubset$Qcorr_qPCR[which(datasubset$Year==2016)])]),datasubset$Qcorr_qPCR[which(datasubset$Year==2016)][order(datasubset$Qcorr_qPCR[which(datasubset$Year==2016)])]),(c(rev(y_up),y_down)),col="grey80",border=NA)
points(datasubset[datasubset$Year == "2016"  , ]$Qcorr_qPCR,datasubset[datasubset$Year == "2016"  , ]$Count, pch = 19)
lines(datasubset$Qcorr_qPCR[which(datasubset$Year==2016)][order(datasubset$Qcorr_qPCR[which(datasubset$Year==2016)])], exp(augment(sockeye.adult.2016.reg)$.fitted)[order(datasubset$Qcorr_qPCR[which(datasubset$Year==2016)])], lwd=2,lty=2, col="blue")

#both years
plot(datasubset$Qcorr_qPCR,datasubset$Count,main="Both",ylab = "",cex.lab=1.5,xlab="",ylim=c(0,600))
y_up=exp(augment(sockeye.adult.reg)$.fitted[order(datasubset$Qcorr_qPCR)] + 2*augment(sockeye.adult.reg)$.se.fit[order(datasubset$Qcorr_qPCR)])
y_down=exp(augment(sockeye.adult.reg)$.fitted[order(datasubset$Qcorr_qPCR)] - 2*augment(sockeye.adult.reg)$.se.fit[order(datasubset$Qcorr_qPCR)])
polygon(c(rev(datasubset$Qcorr_qPCR[order(datasubset$Qcorr_qPCR)]),datasubset$Qcorr_qPCR[order(datasubset$Qcorr_qPCR)]),(c(rev(y_up),y_down)),col="grey80",border=NA)
points(datasubset$Qcorr_qPCR,datasubset$Count, pch = 19)
lines(datasubset$Qcorr_qPCR[order(datasubset$Qcorr_qPCR)], exp(augment(sockeye.adult.reg)$.fitted)[order(datasubset$Qcorr_qPCR)], lwd=2,lty=2, col="blue")

#Coho one outlier with lots of jacks removed
plot(datasubset.coho.trans[datasubset.coho.trans$Year == "2015"  , ]$Qcorr_qPCR,datasubset.coho.trans[datasubset.coho.trans$Year == "2015"  , ]$Count, ylab = "Coho Counts",main="",cex.lab=1.5,xlab="",ylim=c(0,120))
y_up=exp(augment(coho.trans.2015.reg)$.fitted[order(datasubset.coho.trans$Qcorr_qPCR[which(datasubset.coho.trans$Year==2015)])] + 2*augment(coho.trans.2015.reg)$.se.fit[order(datasubset.coho.trans$Qcorr_qPCR[which(datasubset.coho.trans$Year==2015)])])
y_down=exp(augment(coho.trans.2015.reg)$.fitted[order(datasubset.coho.trans$Qcorr_qPCR[which(datasubset.coho.trans$Year==2015)])] - 2*augment(coho.trans.2015.reg)$.se.fit[order(datasubset.coho.trans$Qcorr_qPCR[which(datasubset.coho.trans$Year==2015)])])
polygon(c(rev(datasubset.coho.trans$Qcorr_qPCR[which(datasubset.coho.trans$Year==2015)][order(datasubset.coho.trans$Qcorr_qPCR[which(datasubset.coho.trans$Year==2015)])]),datasubset.coho.trans$Qcorr_qPCR[which(datasubset.coho.trans$Year==2015)][order(datasubset.coho.trans$Qcorr_qPCR[which(datasubset.coho.trans$Year==2015)])]),(c(rev(y_up),y_down)),col="grey80",border=NA)
points(datasubset.coho.trans[datasubset.coho.trans$Year == "2015"  , ]$Qcorr_qPCR,datasubset.coho.trans[datasubset.coho.trans$Year == "2015"  , ]$Count, pch = 19)
lines(datasubset.coho.trans$Qcorr_qPCR[which(datasubset.coho.trans$Year==2015)][order(datasubset.coho.trans$Qcorr_qPCR[which(datasubset.coho.trans$Year==2015)])], exp(augment(coho.trans.2015.reg)$.fitted)[order(datasubset.coho.trans$Qcorr_qPCR[which(datasubset.coho.trans$Year==2015)])], lwd=2,lty=2, col="blue")

plot(datasubset.coho.trans[datasubset.coho.trans$Year == "2016"  , ]$Qcorr_qPCR,datasubset.coho.trans[datasubset.coho.trans$Year == "2016"  , ]$Count,ylab = "",main="",cex.lab=1.5,xlab="",ylim=c(0,120))
y_up=exp(augment(coho.trans.2016.reg)$.fitted[order(datasubset.coho.trans$Qcorr_qPCR[which(datasubset.coho.trans$Year==2016)])] + 2*augment(coho.trans.2016.reg)$.se.fit[order(datasubset.coho.trans$Qcorr_qPCR[which(datasubset.coho.trans$Year==2016)])])
y_down=exp(augment(coho.trans.2016.reg)$.fitted[order(datasubset.coho.trans$Qcorr_qPCR[which(datasubset.coho.trans$Year==2016)])] - 2*augment(coho.trans.2016.reg)$.se.fit[order(datasubset.coho.trans$Qcorr_qPCR[which(datasubset.coho.trans$Year==2016)])])
polygon(c(rev(datasubset.coho.trans$Qcorr_qPCR[which(datasubset.coho.trans$Year==2016)][order(datasubset.coho.trans$Qcorr_qPCR[which(datasubset.coho.trans$Year==2016)])]),datasubset.coho.trans$Qcorr_qPCR[which(datasubset.coho.trans$Year==2016)][order(datasubset.coho.trans$Qcorr_qPCR[which(datasubset.coho.trans$Year==2016)])]),(c(rev(y_up),y_down)),col="grey80",border=NA)
points(datasubset.coho.trans[datasubset.coho.trans$Year == "2016"  , ]$Qcorr_qPCR,datasubset.coho.trans[datasubset.coho.trans$Year == "2016"  , ]$Count, pch = 19)
lines(datasubset.coho.trans$Qcorr_qPCR[which(datasubset.coho.trans$Year==2016)][order(datasubset.coho.trans$Qcorr_qPCR[which(datasubset.coho.trans$Year==2016)])], exp(augment(coho.trans.2016.reg)$.fitted)[order(datasubset.coho.trans$Qcorr_qPCR[which(datasubset.coho.trans$Year==2016)])], lwd=2,lty=2, col="blue")

#both years
plot(datasubset.coho.trans$Qcorr_qPCR,datasubset.coho.trans$Count,ylab = "",main="",cex.lab=1.5,xlab="",ylim=c(0,120))
y_up=exp(augment(coho.trans.reg)$.fitted[order(datasubset.coho.trans$Qcorr_qPCR)] + 2*augment(coho.trans.reg)$.se.fit[order(datasubset.coho.trans$Qcorr_qPCR)])
y_down=exp(augment(coho.trans.reg)$.fitted[order(datasubset.coho.trans$Qcorr_qPCR)] - 2*augment(coho.trans.reg)$.se.fit[order(datasubset.coho.trans$Qcorr_qPCR)])
polygon(c(rev(datasubset.coho.trans$Qcorr_qPCR[order(datasubset.coho.trans$Qcorr_qPCR)]),datasubset.coho.trans$Qcorr_qPCR[order(datasubset.coho.trans$Qcorr_qPCR)]),(c(rev(y_up),y_down)),col="grey80",border=NA)
points(datasubset.coho.trans$Qcorr_qPCR,datasubset.coho.trans$Count, pch = 19)
lines(datasubset.coho.trans$Qcorr_qPCR[order(datasubset.coho.trans$Qcorr_qPCR)], exp(augment(coho.trans.reg)$.fitted)[order(datasubset.coho.trans$Qcorr_qPCR)], lwd=2,lty=2, col="blue")



#sockeye smolts
plot(datasubset.sockeye.juv[datasubset.sockeye.juv$Year == "2015"  , ]$Qcorr_qPCR,datasubset.sockeye.juv[datasubset.sockeye.juv$Year == "2015"  , ]$Count,type="n",ylab="Sockeye Smolts",cex.lab=1.5,xlab="",ylim=c(0,1800))
y_up=exp(augment(sockeye.juv.2015.reg)$.fitted[order(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2015)])] + 2*augment(sockeye.juv.2015.reg)$.se.fit[order(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2015)])])
y_down=exp(augment(sockeye.juv.2015.reg)$.fitted[order(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2015)])] - 2*augment(sockeye.juv.2015.reg)$.se.fit[order(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2015)])])
polygon(c(rev(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2015)][order(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2015)])]),datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2015)][order(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2015)])]),(c(rev(y_up),y_down)),col="grey80",border=NA)
points(datasubset.sockeye.juv[datasubset.sockeye.juv$Year == "2015"  , ]$Qcorr_qPCR,datasubset.sockeye.juv[datasubset.sockeye.juv$Year == "2015"  , ]$Count, pch = 19)
lines(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2015)][order(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2015)])], exp(augment(sockeye.juv.2015.reg)$.fitted)[order(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2015)])], lwd=2,lty=2, col="blue")


plot(datasubset.sockeye.juv[datasubset.sockeye.juv$Year == "2016"  , ]$Qcorr_qPCR,datasubset.sockeye.juv[datasubset.sockeye.juv$Year == "2016"  , ]$Count,ylab = "",cex.lab=1.5,xlab="",ylim=c(0,1800))
y_up=exp(augment(sockeye.juv.2016.reg)$.fitted[order(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2016)])] + 2*augment(sockeye.juv.2016.reg)$.se.fit[order(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2016)])])
y_down=exp(augment(sockeye.juv.2016.reg)$.fitted[order(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2016)])] - 2*augment(sockeye.juv.2016.reg)$.se.fit[order(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2016)])])
polygon(c(rev(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2016)][order(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2016)])]),datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2016)][order(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2016)])]),(c(rev(y_up),y_down)),col="grey80",border=NA)
points(datasubset.sockeye.juv[datasubset.sockeye.juv$Year == "2016"  , ]$Qcorr_qPCR,datasubset.sockeye.juv[datasubset.sockeye.juv$Year == "2016"  , ]$Count, pch = 19)
lines(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2016)][order(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2016)])], exp(augment(sockeye.juv.2016.reg)$.fitted)[order(datasubset.sockeye.juv$Qcorr_qPCR[which(datasubset.sockeye.juv$Year==2016)])], lwd=2,lty=2, col="blue")


plot(datasubset.sockeye.juv$Qcorr_qPCR,datasubset.sockeye.juv$Count,ylab = "",cex.lab=1.5,xlab="",ylim=c(0,1800))
y_up=exp(augment(sockeye.juv.reg)$.fitted[order(datasubset.sockeye.juv$Qcorr_qPCR)] + 2*augment(sockeye.juv.reg)$.se.fit[order(datasubset.sockeye.juv$Qcorr_qPCR)])
y_down=exp(augment(sockeye.juv.reg)$.fitted[order(datasubset.sockeye.juv$Qcorr_qPCR)] - 2*augment(sockeye.juv.reg)$.se.fit[order(datasubset.sockeye.juv$Qcorr_qPCR)])
polygon(c(rev(datasubset.sockeye.juv$Qcorr_qPCR[order(datasubset.sockeye.juv$Qcorr_qPCR)]),datasubset.sockeye.juv$Qcorr_qPCR[order(datasubset.sockeye.juv$Qcorr_qPCR)]),(c(rev(y_up),y_down)),col="grey80",border=NA)
points(datasubset.sockeye.juv$Qcorr_qPCR,datasubset.sockeye.juv$Count, pch = 19)
lines(datasubset.sockeye.juv$Qcorr_qPCR[order(datasubset.sockeye.juv$Qcorr_qPCR)], exp(augment(sockeye.juv.reg)$.fitted)[order(datasubset.sockeye.juv$Qcorr_qPCR)], lwd=2,lty=2, col="blue")

mtext("Flow x DNA",side=1,cex=1.5,outer=T)

#######################################################
#######################################################
#
#Model validation - Predicted vs. actual counts
#
#######################################################
#######################################################

###########
#
#Sockeye Adult
#
###########

par(mfrow=c(3,2),mar=c(3,4.5,1,1),oma=c(2,3,1,1))
plot(datasubset$Date[which(datasubset$Year==2015)], datasubset$Count[which(datasubset$Year==2015)], lwd=2, type='n',pch=19,cex=0.25,ylim=c(0,800),ylab="Sockeye Count",xlab="",cex.lab=1.5, main = "2015")
y_up=exp(augment(sockeye.adult.2015.reg)$.fitted + 2*augment(sockeye.adult.2015.reg)$.se.fit)
y_down=exp(augment(sockeye.adult.2015.reg)$.fitted - 2*augment(sockeye.adult.2015.reg)$.se.fit)
polygon(c(rev(datasubset$Date[which(datasubset$Year==2015)]),datasubset$Date[which(datasubset$Year==2015)]),(c(rev(y_up),y_down)),col="grey80",border=NA)
points(datasubset$Date[which(datasubset$Year==2015)], datasubset$Count[which(datasubset$Year==2015)], lwd=2, type='b',pch=19,cex=0.25)
lines(datasubset$Date[which(datasubset$Year==2015)], exp(augment(sockeye.adult.2015.reg)$.fitted), lwd=2,lty=2, col="blue")




plot(datasubset$Date[which(datasubset$Year==2016)], datasubset$Count[which(datasubset$Year==2016)], lwd=2, type='n',pch=19,cex=0.25,ylim=c(0,800),ylab="",xlab="",cex.lab=1.5, main = "2016")
y_up=exp(augment(sockeye.adult.2016.reg)$.fitted + 2*augment(sockeye.adult.2016.reg)$.se.fit)
y_down=exp(augment(sockeye.adult.2016.reg)$.fitted - 2*augment(sockeye.adult.2016.reg)$.se.fit)
polygon(c(rev(datasubset$Date[which(datasubset$Year==2016)]),datasubset$Date[which(datasubset$Year==2016)]),(c(rev(y_up),y_down)),col="grey80",border=NA)
points(datasubset$Date[which(datasubset$Year==2016)], datasubset$Count[which(datasubset$Year==2016)], lwd=2, type='b',pch=19,cex=0.25)
lines(datasubset$Date[which(datasubset$Year==2016)], exp(augment(sockeye.adult.2016.reg)$.fitted), lwd=2,lty=2, col="blue")

###########
#
#Coho Trans
#
###########

plot(datasubset.coho.trans$Date[which(datasubset.coho.trans$Year==2015)], datasubset.coho.trans$Count[which(datasubset.coho.trans$Year==2015)], lwd=2, type='n',pch=19,cex=0.25,ylim=c(0,150),ylab="Coho Count",xlab="",cex.lab=1.5)
y_up=exp(augment(coho.trans.2015.reg)$.fitted + 2*augment(coho.trans.2015.reg)$.se.fit)
y_down=exp(augment(coho.trans.2015.reg)$.fitted - 2*augment(coho.trans.2015.reg)$.se.fit)
polygon(c(rev(datasubset.coho.trans$Date[which(datasubset.coho.trans$Year==2015)]),datasubset.coho.trans$Date[which(datasubset.coho.trans$Year==2015)]),(c(rev(y_up),y_down)),col="grey80",border=NA)
points(datasubset.coho.trans$Date[which(datasubset.coho.trans$Year==2015)], datasubset.coho.trans$Count[which(datasubset.coho.trans$Year==2015)], lwd=2, type='b',pch=19,cex=0.25)
lines(datasubset.coho.trans$Date[which(datasubset.coho.trans$Year==2015)], exp(augment(coho.trans.2015.reg)$.fitted), lwd=2,lty=2, col="blue")

plot(datasubset.coho.trans$Date[which(datasubset.coho.trans$Year==2016)], datasubset.coho.trans$Count[which(datasubset.coho.trans$Year==2016)], lwd=2, type='n',pch=19,cex=0.25,ylim=c(0,150),ylab="",xlab="",cex.lab=1.5)
y_up=exp(augment(coho.trans.2016.reg)$.fitted + 2*augment(coho.trans.2016.reg)$.se.fit)
y_down=exp(augment(coho.trans.2016.reg)$.fitted - 2*augment(coho.trans.2016.reg)$.se.fit)
polygon(c(rev(datasubset.coho.trans$Date[which(datasubset.coho.trans$Year==2016)]),datasubset.coho.trans$Date[which(datasubset.coho.trans$Year==2016)]),(c(rev(y_up),y_down)),col="grey80",border=NA)
points(datasubset.coho.trans$Date[which(datasubset.coho.trans$Year==2016)], datasubset.coho.trans$Count[which(datasubset.coho.trans$Year==2016)], lwd=2, type='b',pch=19,cex=0.25)
lines(datasubset.coho.trans$Date[which(datasubset.coho.trans$Year==2016)], exp(augment(coho.trans.2016.reg)$.fitted), lwd=2,lty=2, col="blue")
#mtext("coho trans Counts Quasipoisson model",side=2,cex=2,outer=T)
#text(datasubset.coho.trans$Date[which(datasubset.coho.trans$Year==2016)][35],800,"2016",cex=1.5)

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


plot(datasubset.sockeye.juv$Date[which(datasubset.sockeye.juv$Year==2016)], datasubset.sockeye.juv$Count[which(datasubset.sockeye.juv$Year==2016)], lwd=2, type='n',pch=19,cex=0.25,ylim=c(0,2000),ylab="",xlab="",cex.lab=1.5)
y_up=exp(augment(sockeye.juv.2016.reg)$.fitted + 2*augment(sockeye.juv.2016.reg)$.se.fit)
y_down=exp(augment(sockeye.juv.2016.reg)$.fitted - 2*augment(sockeye.juv.2016.reg)$.se.fit)
polygon(c(rev(datasubset.sockeye.juv$Date[which(datasubset.sockeye.juv$Year==2016)]),datasubset.sockeye.juv$Date[which(datasubset.sockeye.juv$Year==2016)]),(c(rev(y_up),y_down)),col="grey80",border=NA)
points(datasubset.sockeye.juv$Date[which(datasubset.sockeye.juv$Year==2016)], datasubset.sockeye.juv$Count[which(datasubset.sockeye.juv$Year==2016)], lwd=2, type='b',pch=19,cex=0.25)
lines(datasubset.sockeye.juv$Date[which(datasubset.sockeye.juv$Year==2016)], exp(augment(sockeye.juv.2016.reg)$.fitted), lwd=2,lty=2, col="blue")


#####################
#Begin:  Models of eDNA
#####################
#Quantify lag with sockeye due to simple life history
sockeye.eDNA.reg <- lm(Qcorr_qPCR[6:39] ~ Count[6:39], data = datasubset[datasubset$Year==2016,])
summary(sockeye.eDNA.reg)
sockeye.eDNA.reg2 <- lm(resid(sockeye.eDNA.reg) ~ Count[5:38], data = datasubset[datasubset$Year==2016,])
summary(sockeye.eDNA.reg2)
sockeye.eDNA.reg3 <- lm(resid(sockeye.eDNA.reg2) ~ Count[4:37], data = datasubset[datasubset$Year==2016,])
summary(sockeye.eDNA.reg3)

#quantify life history influence on eDNA with coho
datasubset.eDNA.coho = subset(coho_adult, !is.na(QUANT_mean) & Year %in% c("2015", "2016"))
datasubset.eDNA.coho <- datasubset.eDNA.coho %>% spread(Cohotype,Count)

coho.eDNA2015.reg <- lm(Qcorr_qPCR ~ Coho_Adults + Coho_Jack + Coho_Juv, data = datasubset.eDNA.coho[datasubset.eDNA.coho$Year==2015,])
coho.eDNA2016.reg <- lm(Qcorr_qPCR ~ Coho_Adults + Coho_Jack + Coho_Juv, data = datasubset.eDNA.coho[datasubset.eDNA.coho$Year==2016,])
coho.eDNA.reg <- lm(Qcorr_qPCR ~ Coho_Adults + Coho_Jack + Coho_Juv, data = datasubset.eDNA.coho)

#####################
#End:  Models of eDNA
#####################

