################################################################################################
# R code to visualise (only) the index of cancer survival as described in the tutorial:
# An index of cancer survival to measure progress in cancer control: A tutorial
# https://doi.org/10.1016/j.canep.2024.102576
# Authors: Manuela Quaresma and F. Javier Rubio
# More details at: https://github.com/ManuelaQuaresma/CSI
################################################################################################


# Delete memory
rm(list=ls())

# Required packages
#options("install.lock"=FALSE)
#install.packages("ggplot2", dependencies = TRUE)


library(ggplot2)

# Set your Working Directory 
# From RStudio, use the menu to change your working directory under 
# Session > Set Working Directory > Choose Directory.

load("cancer_index.RData")

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



