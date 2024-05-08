################################################################################################
# R code to calculate an index of cancer survival as described in the tutorial:
# An index of cancer survival to measure progress in cancer control: A tutorial
# https://doi.org/10.1016/j.canep.2024.102576
# Authors: Manuela Quaresma and F. Javier Rubio
# More details at: https://github.com/ManuelaQuaresma/CSI
################################################################################################


# Delete memory
rm(list=ls())

# Required packages
#options("install.lock"=FALSE)
#install.packages("GJRM", dependencies = TRUE)
#install.packages("knitr", dependencies = TRUE)
#install.packages("tmvnsim", dependencies = TRUE)
#install.packages("numDeriv", dependencies = TRUE)
#install.packages("ggplot2", dependencies = TRUE)

library(knitr)
library(numDeriv)
library(tmvnsim)
library(GJRM)
library(ggplot2)

# Set your Working Directory 
# From RStudio, use the menu to change your working directory under 
# Session > Set Working Directory > Choose Directory.


################################################################################################
# Read the data and weights
# Files available at: https://github.com/ManuelaQuaresma/CSI
################################################################################################

# Read the Replica data set and weights
df <- read.table("Replica_03052024.txt", header = TRUE)
weights <- read.table("weights_cancer_age_sex_specific_03052024.txt", header = TRUE)


# Summaries
dim(df)
head(df)

dim(weights)
head(weights)

table(df$cancer)
table(df$period)

################################################################################################
# Data preparation
################################################################################################

# scaled age at diagnosis
df$agec <- scale(df$age)

# Design matrix for "period" of diagnosis
period <- matrix(0, ncol = 5, nrow = nrow(df))
colnames(period) <- c("period.1", "period.2", "period.3", "period.4", "period.5")

for(i in 1:nrow(df)) period[i,df$period[i]] = 1
df <- as.data.frame(cbind(df, period))

# sample size
n <- nrow(df) 
n

# names and number of cancer sites
can_names <- unique(df$cancer)
ncs <- length(unique(df$cancer))

# Data frame for women
dfw <- df[df$sex == 2,]
can_namesw <- unique(dfw$cancer)
ncsw <- length(can_namesw)

# Data frame for men
dfm <- df[df$sex == 1,]
can_namesm <- unique(dfm$cancer)
ncsm <- length(can_namesm)

# number of age groups
nage <- 5

# age group limits
a11 <- 15; a12 <- 44
a21 <- 45; a22 <- 54
a31 <- 55; a32 <- 64
a41 <- 65; a42 <- 74
a51 <- 75; a52 <- 99

# Names and Number of periods of diagnosis
periods <- sort(unique(df$period))
nperiod <- length(unique(df$period))

################################################################################################
################################################################################################
# Analysis for Women data
################################################################################################
################################################################################################


################################################################################################
# Model fit
################################################################################################
Start <- Sys.time()
Start

model_women <- list() 
ind.minw <- vector()

# Model formula
eq_GJRM = list(time ~ s(log(time), bs = "mpi") + s(agec, bs='cr') + ti(log(time), agec) + 
                 period.2 + period.3 + period.4 + period.5)

for(j in 1:ncsw){
  print(j)
  # Data sets for each cancer site 
  indc <- which(dfw$cancer == can_namesw[j])
  dfc <- dfw[indc,]
  
  # Required quantities
  statusc <- as.vector(dfc$status)
  # Population hazard rates for the uncensored individuals
  hrate.select <- dfc$brate[as.logical(statusc)]
  
  # Optimisation step
  
  # Generalised proportional hazards model
  out_GJRM_1 = gamlss(eq_GJRM, data = dfc, surv = TRUE, margin = 'PH',
                      cens = status, type.cens = "R", hrate = hrate.select)
  
  # Generalised proportional odds model
  out_GJRM_2 = gamlss(eq_GJRM, data = dfc, surv = TRUE, margin = 'PO',
                      cens = status, type.cens = "R", hrate = hrate.select)
  
  # Generalised probit model
  out_GJRM_3 = gamlss(eq_GJRM, data = dfc, surv = TRUE, margin = 'probit',
                      cens = status, type.cens = "R", hrate = hrate.select)
  
  
  # Comparing models 
  ind.minw[j] <- which.min(c(AIC(out_GJRM_1), AIC(out_GJRM_2), AIC(out_GJRM_3)))
  
  if(ind.minw[j] == 1) model_women[[j]] = out_GJRM_1
  if(ind.minw[j] == 2) model_women[[j]] = out_GJRM_2
  if(ind.minw[j] == 3) model_women[[j]] = out_GJRM_3
  
}

################################################################################################
# Weighted net survival
################################################################################################

# Net survival at 1,5,10 years
# For each age group and period of diagnosis

ts <- c(1,5,10)

# Initialising list containing the net survival by period for all cancer sites 
lnames <- paste("NSW",periods,sep="")
NSW <- vector("list", length(lnames))
names(NSW) <- lnames

templist <- list()
MAT <- matrix(0, ncol = length(ts), nrow = nage)
colnames(MAT) <- paste("time",1:length(ts))
rownames(MAT) <- paste("group",1:nage)

for(i in 1:ncsw){
  templist[[i]] <- list(MAT)
  names(templist[[i]]) <- can_namesw[i]
}

for(j in 1:length(lnames)) NSW[[j]] <- templist

rm(templist)

# Populating the net survival list
# period of diagnosis and cancer site for loop
for(i in 1:nperiod){
  print(periods[i])
  for(j in 1:ncsw){
    # Data sets for each cancer site 
    indc <- which(dfw$cancer == can_namesw[j])
    dfc <- dfw[indc,]
    # Data sets for each cancer site and period
    indcy <- which(dfw$cancer == can_namesw[j] & dfw$period == periods[i])
    dfcy <- dfw[indcy,]
    
    # age groups for each cancer site
    inda1 <- which(a11 <= dfcy$age & dfcy$age <= a12)
    inda2 <- which(a21 <= dfcy$age & dfcy$age <= a22)
    inda3 <- which(a31 <= dfcy$age & dfcy$age <= a32)
    inda4 <- which(a41 <= dfcy$age & dfcy$age <= a42)
    inda5 <- which(a51 <= dfcy$age & dfcy$age <= a52)
    
    ind_age <- list(inda1 = inda1, inda2 = inda2, inda3 = inda3, inda4 = inda4, inda5 = inda5)
    
    # Calculating net survival for each age group 
    for(k in 1:nage){
      for(r in 1:length(ts)){
      # NS for cancer site j, period i
      NSW[[i]][[j]][[1]][k,r] <- hazsurv(model_women[[j]], type = 'surv', newdata = dfcy[ind_age[[k]],], t.vec = ts[r],
                                             ls = 1, intervals = FALSE, n.sim = 1000, plot.out = FALSE, print.progress = FALSE)$s
      }  
    }
    
  }
}


# Calculating the weighted net survival
wnsnames <- paste("wnsw",periods,sep="")
wnsw <- vector("list", length(wnsnames))
names(wnsw) <- wnsnames

MAT <- matrix(0, ncol = length(ts), nrow = ncsw)
colnames(MAT) <- paste("time",1:length(ts))
rownames(MAT) <- can_namesw

for(j in 1:nperiod) wnsw[[j]] <- MAT

# weights for women
ww <- weights[weights$sex == 2, ]

for( i in 1:nperiod){
  for(j in 1:ncsw){
    wcy <- ww$stand_weights[ww$cancer == can_namesw[j]]
    wnsw[[i]][j,] <- as.vector(t(NSW[[i]][[j]][[1]])%*%wcy)
  }
}


################################################################################################
################################################################################################
# Analysis for Men data
################################################################################################
################################################################################################


################################################################################################
# Model fit
################################################################################################

model_men <- list() 
ind.minm <- vector()

# Model formula
eq_GJRM = list(time ~ s(log(time), bs = "mpi") + s(agec, bs='cr') + ti(log(time), agec) + 
                 period.2 + period.3 + period.4 + period.5) 

for(j in 1:ncsm){
  print(j)
  # Data sets for each cancer site 
  indc <- which(dfm$cancer == can_namesm[j])
  dfc <- dfm[indc,]
  
  # Required quantities
  statusc <- as.vector(dfc$status)
  # Population hazard rates for the uncensored individuals
  hrate.select <- dfc$brate[as.logical(statusc)]
  
  # Optimisation step
  
  # Generalised proportional hazards model
  out_GJRM_1 = gamlss(eq_GJRM, data = dfc, surv = TRUE, margin = 'PH',
                      cens = status, type.cens = "R", hrate = hrate.select)
  
  # Generalised proportional odds model
  out_GJRM_2 = gamlss(eq_GJRM, data = dfc, surv = TRUE, margin = 'PO',
                      cens = status, type.cens = "R", hrate = hrate.select)
  
  # Generalised probit model
  out_GJRM_3 = gamlss(eq_GJRM, data = dfc, surv = TRUE, margin = 'probit',
                      cens = status, type.cens = "R", hrate = hrate.select)
  
  
  # Comparing models 
  ind.minm[j] <- which.min(c(AIC(out_GJRM_1), AIC(out_GJRM_2), AIC(out_GJRM_3)))
  
  if(ind.minm[j] == 1) model_men[[j]] = out_GJRM_1
  if(ind.minm[j] == 2) model_men[[j]] = out_GJRM_2
  if(ind.minm[j] == 3) model_men[[j]] = out_GJRM_3
  
}

################################################################################################
# Weighted net survival
################################################################################################

# Net survival at 1,5,10 years
# For each age group and period of diagnosis

ts <- c(1,5,10)

# Initialising list containing the net survival by period for all cancer sites 
lnames <- paste("NSM",periods,sep="")
NSM <- vector("list", length(lnames))
names(NSM) <- lnames

templist <- list()
MAT <- matrix(0, ncol = length(ts), nrow = nage)
colnames(MAT) <- paste("time",1:length(ts))
rownames(MAT) <- paste("group",1:nage)

for(i in 1:ncsm){
  templist[[i]] <- list(MAT)
  names(templist[[i]]) <- can_namesm[i]
}

for(j in 1:length(lnames)) NSM[[j]] <- templist

rm(templist)

# Populating the net survival list
# period of diagnosis and cancer site for loop
for(i in 1:nperiod){
  print(periods[i])
  for(j in 1:ncsm){
    # Data sets for each cancer site 
    indc <- which(dfm$cancer == can_namesm[j])
    dfc <- dfm[indc,]
    # Data sets for each cancer site and period
    indcy <- which(dfm$cancer == can_namesm[j] & dfm$period == periods[i])
    dfcy <- dfm[indcy,]
    
    # age groups for each cancer site
    inda1 <- which(a11 <= dfcy$age & dfcy$age <= a12)
    inda2 <- which(a21 <= dfcy$age & dfcy$age <= a22)
    inda3 <- which(a31 <= dfcy$age & dfcy$age <= a32)
    inda4 <- which(a41 <= dfcy$age & dfcy$age <= a42)
    inda5 <- which(a51 <= dfcy$age & dfcy$age <= a52)
    
    ind_age <- list(inda1 = inda1, inda2 = inda2, inda3 = inda3, inda4 = inda4, inda5 = inda5)
    
    # Calculating net survival for each age group 
    for(k in 1:nage){
      for(r in 1:length(ts)){
        # NS for cancer site j, period i
        NSM[[i]][[j]][[1]][k,r] <- hazsurv(model_men[[j]], type = 'surv', newdata = dfcy[ind_age[[k]],], t.vec = ts[r],
                                                ls = 1, intervals = TRUE, n.sim = 1000, plot.out = FALSE, print.progress = FALSE)$s
      }  
    }
    
  }
}


# Calculating the weighted net survival
wnsnames <- paste("wnsm",periods,sep="")
wnsm <- vector("list", length(wnsnames))
names(wnsm) <- wnsnames

MAT <- matrix(0, ncol = length(ts), nrow = ncsm)
colnames(MAT) <- paste("time",1:length(ts))
rownames(MAT) <- can_namesm

for(j in 1:nperiod) wnsm[[j]] <- MAT

# weights for men
wm <- weights[weights$sex == 1, ]

for( i in 1:nperiod){
  for(j in 1:ncsm){
    wcy <- wm$stand_weights[wm$cancer == can_namesm[j]]
    wnsm[[i]][j,] <- as.vector(t(NSM[[i]][[j]][[1]])%*%wcy)
  }
}


###################################################################################
# Cancer index calculation
###################################################################################

cancer_index <- mapply("+", lapply(wnsw,colSums), lapply(wnsm,colSums), SIMPLIFY = FALSE) 

names(cancer_index) <- paste("cancer index", periods)

cancer_index

End <- Sys.time()

Start
End


###################################################################################
# Cancer index visualisation
###################################################################################

# cancer index time series matrix
CIM = cbind(cancer_index$`cancer index 1`,
            cancer_index$`cancer index 2`,
            cancer_index$`cancer index 3`,
            cancer_index$`cancer index 4`,
            cancer_index$`cancer index 5`)

# Cancer index data frame for plotting
CSI = data.frame( NS1 = as.vector(CIM[1,])*100, 
                  NS5 = as.vector(CIM[2,])*100, 
                  NS10 = as.vector(CIM[3,])*100,
                  periods = c("1980-84","1985-89","1990-94","1995-99","2000-04"))

surv = c("1-year Survival", "5-year Survival", "10-year Survival")
surv = factor(surv, levels = surv)




# Create the basic plot 
pcsi = ggplot() +
  geom_point(data = CSI, aes(x = 1:5, y = NS1, color = surv[1]), shape = 18, size = 5) + 
  geom_line(data = CSI, aes(x = 1:5, y = NS1, color = surv[1]), linewidth = 1, colour = "black") +
  geom_point(data = CSI, aes(x = 1:5, y = NS5, color = surv[2]), shape = 19, size = 4) +
  geom_line(data = CSI, aes(x = 1:5, y = NS5, color = surv[2]), linetype = "dashed", linewidth = 1, colour = "black") +
  geom_point(data = CSI, aes(x = 1:5, y = NS10, color = surv[3]), shape = 15, size = 4) +
  geom_line(data = CSI, aes(x = 1:5, y = NS10, color = surv[3]), linetype = "dotted", linewidth = 1, colour = "black") +
  scale_x_discrete(limits = CSI$periods) + 
  ylim(0,100) +
  scale_color_manual(values = c("black", "black", "black"), breaks = c(surv[1], surv[2], surv[3])) +
  guides(color = guide_legend(override.aes=list(shape = c(18,19,15)))) 

# Create the plot showing the values at each point
pcsi2 = ggplot() +
  geom_point(data = CSI, aes(x = 1:5, y = NS1, color = surv[1]), shape = 18, size = 5) + 
  geom_line(data = CSI, aes(x = 1:5, y = NS1, color = surv[1]), linewidth = 1, colour = "black") +
  scale_x_discrete(limits = CSI$periods) + 
  geom_text(aes(x = 1:5, y = CSI$NS1,label = round(CSI$NS1,digits = 2)), nudge_y = 7) +
  geom_point(data = CSI, aes(x = 1:5, y = NS5, color = surv[2]), shape = 19, size = 4) +
  geom_line(data = CSI, aes(x = 1:5, y = NS5, color = surv[2]), linetype = "dashed", linewidth = 1, colour = "black") +
  geom_text(aes(x = 1:5, y = CSI$NS5,label = round(CSI$NS5,digits = 2)), nudge_y = 7) +
  geom_point(data = CSI, aes(x = 1:5, y = NS10, color = surv[3]), shape = 15, size = 4) +
  geom_line(data = CSI, aes(x = 1:5, y = NS10, color = surv[3]), linetype = "dotted", linewidth = 1, colour = "black") +
  geom_text(aes(x = 1:5, y = CSI$NS10,label = round(CSI$NS10,digits = 2)), nudge_y = -7) +
  ylim(0,100) +
  scale_color_manual(values = c("black", "black", "black"), breaks = c(surv[1], surv[2], surv[3])) +
  guides(color = guide_legend(override.aes=list(shape = c(18,19,15)))) 

# Creating labels

csilabs =   labs(x="Period of diagnosis", 
                 y="Index of net survival (%)",
                 title="Trends in the index of net survival",
                 subtitle="Using 'sex-cancer-age weights' from the 2000-2004 cohort")

# Creating theme

csitheme = theme(plot.title = element_text(hjust = 0.5, size = 16), 
                 plot.subtitle = element_text(hjust = 0.5),
                 axis.title.x = element_text(face="bold", color="black", 
                                             size=12, angle=0, vjust = 0.5),
                 axis.title.y = element_text(face="bold", color="black", 
                                             size=12, angle=90, vjust = 0.5),
                 axis.text.x = element_text(face="bold", color="black", 
                                            size=12, angle=45, vjust = 0.5),
                 axis.text.y = element_text(face="bold", color="black", 
                                            size=12, angle=45),
                 panel.background = element_rect(fill = "white",
                                                 colour = "gray",
                                                 size = 0.5, linetype = "solid"),
                 panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                                 colour = "gray"), 
                 panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                                 colour = "gray"),
                 panel.border = element_rect(colour = "black", fill=NA, size=1),
                 legend.position = c(0.85,0.15),
                 legend.title=element_blank(),
                 legend.box.background = element_rect(colour = "black")) 

# Adding plot styles

# Basic
pcsi + csilabs + csitheme 

# Showing values at each point
pcsi2 + csilabs + csitheme 


   
