############################################## ROBVIS PLOTS #######################################################

# First we import the ROBINS dataset from the "Data" folder on desktop
ROBINS<-read.csv("C:/Users/sotbi/Desktop/Data/ROBINS.csv", sep=";", header=TRUE)

# Now we will use the ggplot2 & robvis (development) packages to visualize the quality assessment
install.packages("ggplot2")
library(ggplot2)
install.packages("gridExtra")
library(gridExtra)
install.packages("devtools")
library(devtools)
devtools::install_github("mcguinlu/robvis")
install_github("mcguinlu/robvis")
library(robvis)

# Formulate a new data.frame for robvis
rob_dset<-ROBINS[,c("author", "center", "matched", "post_2018", "D1", "D2", "D3", "D4", "D5", "D6", "D7", "robins")]
# Change the name of the column "author" into "study"
colnames(rob_dset)[1] <- "study"


# We create a new folder named "Robvis" on desktop 
folder_name <- "Robvis"
# Specify the desktop path
desktop_path <- file.path("C:/Users/sotbi/Desktop")
# Create the folder path
folder_path <- file.path(desktop_path, folder_name)

# Check if the folder already exists
if (dir.exists(folder_path)) {
  print("Folder already exists!")
} else {
  # Create the folder on the desktop
  dir.create(folder_path)
  print("Folder created successfully!")
}

# Set working space to the desktop folder that was created
setwd("C:/Users/sotbi/Desktop/Robvis/")


# Now we define the subgroups from "rob_dset"
# Then we create the plots and save them in "Robvis" folder

# To capture the size of the plots that will be produced we adjust the width & height parameters from here
# Pooled analysis parameters
# For traffic light plots
W_traf_pool <- 1400
H_traf_pool <- 990
# For summary plots
W_sum_pool <- 1200
H_sum_pool <- 740

# Pooled studies
rob_pooled<-rob_dset[,-c(which(colnames(rob_dset)=="center"), which(colnames(rob_dset)=="matched"), which(colnames(rob_dset)=="post_2018"))]

# Traffic light plot
png("1. ROBINS-I traffic light plot for all studies.png", width = W_traf_pool, height = H_traf_pool)
rob_traffic_light(rob_pooled, tool = "ROBINS-I", psize = 5) +
  ggtitle("ROBINS-I traffic light plot for all studies") +
  theme(plot.title = element_text(size = 18, hjust = 0.5)) +
  theme(legend.text = element_text(size = 15)) +
  theme(plot.margin = margin(10, 10, 10, 15))
dev.off() 
# Summary plot
png("2. ROBINS-I assessment summary for all studies.png", width = W_sum_pool, height = H_sum_pool)
rob_summary(rob_pooled, tool = "ROBINS-I", overall = TRUE, weighted = FALSE) +
  ggtitle("ROBINS-I assessment summary for all studies") +
  theme(plot.title = element_text(size = 20, hjust = 0.5)) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 17)) +
  theme(legend.text = element_text(size = 14)) +
  theme(plot.margin = margin(30, 30, 30, 30)) 
dev.off()
 
# similar plots will be produced so we adjust the width & height parameters from here
# Subgroup analysis parameters
# For traffic light plots
W_traf_sga <- 1300
H_traf_sga <- 920
# For summary plots
W_sum_sga <- 1000
H_sum_sga <- 800

# For studies published after or before 2018
# Studies published after 2018
rob_post_2018<-rob_dset[rob_dset$post_2018==1,]
rob_post_2018<-rob_post_2018[,-c(which(colnames(rob_post_2018)=="center"), which(colnames(rob_post_2018)=="matched"), which(colnames(rob_post_2018)=="post_2018"))]
# Studies published before 2018
rob_pre_2018<-rob_dset[rob_dset$post_2018==0,]
rob_pre_2018<-rob_pre_2018[,-c(which(colnames(rob_pre_2018)=="center"), which(colnames(rob_pre_2018)=="matched"), which(colnames(rob_pre_2018)=="post_2018"))]

# Traffic light plot for studies published after 2018
png("3. ROBINS-I traffic light plot for studies published after 2018.png", width = W_traf_sga, height = H_traf_sga)
rob_traffic_light(rob_post_2018, tool = "ROBINS-I", psize = 7) +
  ggtitle("ROBINS-I traffic light plot for studies published after 2018") +
  theme(plot.title = element_text(size = 18, hjust = 0.5)) +
  theme(legend.text = element_text(size = 15)) +
  theme(plot.margin = margin(10, 10, 10, 15))
dev.off() 
# Traffic light plot for studies published before 2018
png("4. ROBINS-I traffic light plot for studies published before 2018.png", width = W_traf_sga, height = H_traf_sga)
rob_traffic_light(rob_pre_2018, tool = "ROBINS-I", psize = 7) +
  ggtitle("ROBINS-I traffic light plot for studies published before 2018") +
  theme(plot.title = element_text(size = 18, hjust = 0.5)) +
  theme(legend.text = element_text(size = 15)) +
  theme(plot.margin = margin(10, 10, 10, 15))
dev.off()  

png("5. ROBINS-I assessment summary for studies published after or before 2018.png", width = W_sum_sga, height = H_sum_sga)
# Summary plot for studies published after 2018
plot_1 <- rob_summary(rob_post_2018, tool = "ROBINS-I", overall = TRUE, weighted = FALSE) +
  ggtitle("ROBINS-I assessment summary for studies published after 2018") +
  theme(plot.title = element_text(size = 20, hjust = 0.5)) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 17)) +
  theme(legend.text = element_text(size = 14)) +
  theme(plot.margin = margin(30, 30, 30, 30))
# Summary plot for studies published before 2018
plot_2 <- rob_summary(rob_pre_2018, tool = "ROBINS-I", overall = TRUE, weighted = FALSE) +
  ggtitle("ROBINS-I assessment summary for studies published before 2018") +
  theme(plot.title = element_text(size = 20, hjust = 0.5)) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 17)) +
  theme(legend.text = element_text(size = 14)) +
  theme(plot.margin = margin(30, 30, 30, 30)) 
grid.arrange(plot_1, plot_2, ncol = 1)
dev.off()


# Now for studies with or without matching
# Studies with matching
rob_matched<-rob_dset[rob_dset$matched==1,]
rob_matched<-rob_matched[,-c(which(colnames(rob_matched)=="center"), which(colnames(rob_matched)=="matched"), which(colnames(rob_matched)=="post_2018"))]
# Studies without matching
rob_nonmatched<-rob_dset[rob_dset$matched==0,]
rob_nonmatched<-rob_nonmatched[,-c(which(colnames(rob_nonmatched)=="center"), which(colnames(rob_nonmatched)=="matched"), which(colnames(rob_nonmatched)=="post_2018"))]

# Traffic light plot for studies with matching
png("6. ROBINS-I traffic light plot for studies with patient matching.png", width = W_traf_sga, height = H_traf_sga)
rob_traffic_light(rob_matched, tool = "ROBINS-I", psize = 7) +
  ggtitle("ROBINS-I traffic light plot for studies with patient matching") +
  theme(plot.title = element_text(size = 18, hjust = 0.5)) +
  theme(legend.text = element_text(size = 15)) +
  theme(plot.margin = margin(10, 10, 10, 15))
dev.off() 
# Traffic light plot for studies without matching
png("7. ROBINS-I traffic light plot for studies without patient matching.png", width = W_traf_sga, height = H_traf_sga)
rob_traffic_light(rob_nonmatched, tool = "ROBINS-I", psize = 7) +
  ggtitle("ROBINS-I traffic light plot for studies without patient matching") +
  theme(plot.title = element_text(size = 18, hjust = 0.5)) +
  theme(legend.text = element_text(size = 15)) +
  theme(plot.margin = margin(10, 10, 10, 15))
dev.off()  

png("8. ROBINS-I assessment summary for studies with or without patient matching.png", width = W_sum_sga, height = H_sum_sga)
# Summary plot for studies with matching
plot_3 <- rob_summary(rob_matched, tool = "ROBINS-I", overall = TRUE, weighted = FALSE) +
  ggtitle("ROBINS-I assessment summary for studies with patient matching") +
  theme(plot.title = element_text(size = 20, hjust = 0.5)) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 17)) +
  theme(legend.text = element_text(size = 14)) +
  theme(plot.margin = margin(30, 30, 30, 30)) 
# Summary plot for studies without matching
plot_4 <- rob_summary(rob_nonmatched, tool = "ROBINS-I", overall = TRUE, weighted = FALSE) +
  ggtitle("ROBINS-I assessment summary for studies without patient matching") +
  theme(plot.title = element_text(size = 20, hjust = 0.5)) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 17)) +
  theme(legend.text = element_text(size = 14)) +
  theme(plot.margin = margin(30, 30, 30, 30)) 
grid.arrange(plot_3, plot_4, ncol = 1)
dev.off()


# Finally for multi- or single-center studies
# Multicenter studies
rob_multi<-rob_dset[rob_dset$center=="multi",]
rob_multi<-rob_multi[,-c(which(colnames(rob_multi)=="center"), which(colnames(rob_multi)=="matched"), which(colnames(rob_multi)=="post_2018"))]
# Single-center studies
rob_single<-rob_dset[rob_dset$center=="single",]
rob_single<-rob_single[,-c(which(colnames(rob_single)=="center"), which(colnames(rob_single)=="matched"), which(colnames(rob_single)=="post_2018"))]

# Traffic light plot for multicenter studies
png("9. ROBINS-I traffic light plot for multicenter studies.png", width = W_traf_sga, height = H_traf_sga)
rob_traffic_light(rob_multi, tool = "ROBINS-I", psize = 7) +
  ggtitle("ROBINS-I traffic light plot for multicenter studies") +
  theme(plot.title = element_text(size = 18, hjust = 0.5)) +
  theme(legend.text = element_text(size = 15)) +
  theme(plot.margin = margin(10, 10, 10, 15))
dev.off()  
# Traffic light plot for single-center studies
png("10. ROBINS-I traffic light plot for single-center studies.png", width = W_traf_sga, height = H_traf_sga)
rob_traffic_light(rob_single, tool = "ROBINS-I", psize = 7) +
  ggtitle("ROBINS-I traffic light plot for single-center studies") +
  theme(plot.title = element_text(size = 18, hjust = 0.5)) +
  theme(legend.text = element_text(size = 15)) +
  theme(plot.margin = margin(10, 10, 10, 15))
dev.off()  


png("11. ROBINS-I assessment summary for multi- or single-center studies.png", width = W_sum_sga, height = H_sum_sga)
# Summary plot for multicenter studies
plot_5 <- rob_summary(rob_multi, tool = "ROBINS-I", overall = TRUE, weighted = FALSE) +
  ggtitle("ROBINS-I assessment summary for multicenter studies") +
  theme(plot.title = element_text(size = 20, hjust = 0.5)) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 17)) +
  theme(legend.text = element_text(size = 14)) +
  theme(plot.margin = margin(30, 30, 30, 30)) 
# Summary plot for single-center studies
plot_6 <- rob_summary(rob_single, tool = "ROBINS-I", overall = TRUE, weighted = FALSE) +
  ggtitle("ROBINS-I assessment summary for single-center studies") +
  theme(plot.title = element_text(size = 20, hjust = 0.5)) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 17)) +
  theme(legend.text = element_text(size = 14)) +
  theme(plot.margin = margin(30, 30, 30, 30))
grid.arrange(plot_5, plot_6, ncol = 1)
dev.off()


################################################ NOS PLOTS #########################################################

# Formulate a new data.frame for NOS (Newcastle - Ottawa Scale)
nos_dset<-ROBINS[,c("author", "center", "matched", "post_2018", "stars")]

# We create a new folder named "NOS" on desktop 
folder_name <- "NOS"
# Specify the desktop path
desktop_path <- file.path("C:/Users/sotbi/Desktop")
# Create the folder path
folder_path <- file.path(desktop_path, folder_name)

# Check if the folder already exists
if (dir.exists(folder_path)) {
  print("Folder already exists!")
} else {
  # Create the folder on the desktop
  dir.create(folder_path)
  print("Folder created successfully!")
}

# Set working space to the desktop folder that was created
setwd("C:/Users/sotbi/Desktop/NOS/")

# Plot parameters
# Pooled studies
W_pool <- 1100
H_pool <- 780

# Subgroups
W_sga <- 1000
H_sga <- 1250

# Pooled studies
png("1. NOS rating histogram in pooled studies.png", width = W_pool, height = H_pool)
par(mar = c(5, 6, 4, 0) + 0.1, mgp = c(3.25, 1, 0))
# Create a histogram with each unique value in NOS rating as a bin, omitting the first 3 bins
hist_data <- hist(nos_dset$stars, 
                  breaks = seq(3, 9, by = 1),
                  plot = FALSE)
# Convert frequencies to percentages
percentages <- prop.table(hist_data$counts) * 100
# Plot the histogram with percentages on the y-axis
barplot(percentages, 
        names.arg = hist_data$breaks[-1], 
        col = "brown1", 
        main = "NOS rating in pooled studies", 
        xlab = "Quality stars", 
        ylab = "Percentage of data (%)",
        xlim = c(0, 9),
        ylim = c(0, 60),
        cex.main = 1.7,
        cex.lab = 1.5,
        cex.axis = 1.5,
        cex.names = 1.5)
dev.off()


# Studies published before or after 2018
png("2. NOS rating histogram in studies published after or before 2018.png", width = W_sga, height = H_sga)
par(mfrow = c(2, 1), mar = c(5, 6, 4, 0) + 0.1, mgp = c(3.25, 1, 0))
hist_data_post_2018 <- hist(nos_dset[nos_dset$post_2018==1,]$stars, 
                            breaks = seq(3, 9, by = 1),
                            plot = FALSE)
percentages_post_2018 <- prop.table(hist_data_post_2018$counts) * 100
barplot(percentages_post_2018, 
        names.arg = hist_data_post_2018$breaks[-1], 
        col = "darkgoldenrod4", 
        main = "NOS rating in studies published after 2018", 
        xlab = "Quality stars", 
        ylab = "Percentage of data (%)",
        xlim = c(0, 9),
        ylim = c(0, 60),
        cex.main = 2,
        cex.lab = 1.5,
        cex.axis = 1.5,
        cex.names = 1.5)

hist_data_pre_2018 <- hist(nos_dset[nos_dset$post_2018==0,]$stars, 
                           breaks = seq(3, 9, by = 1),
                           plot = FALSE)
percentages_pre_2018 <- prop.table(hist_data_pre_2018$counts) * 100
barplot(percentages_pre_2018, 
        names.arg = hist_data_pre_2018$breaks[-1], 
        col = "darkgoldenrod2", 
        main = "NOS rating in studies published before 2018", 
        xlab = "Quality stars", 
        ylab = "Percentage of data (%)",
        xlim = c(0, 9),
        ylim = c(0, 60),
        cex.main = 2,
        cex.lab = 1.5,
        cex.axis = 1.5,
        cex.names = 1.5)
dev.off()


# Studies with or without patient matching
png("3. NOS rating histogram in studies with or without patient matching.png", width = W_sga, height = H_sga)
par(mfrow = c(2, 1), mar = c(5, 6, 4, 0) + 0.1, mgp = c(3.25, 1, 0))
hist_data_matched <- hist(nos_dset[nos_dset$matched==1,]$stars, 
                          breaks = seq(3, 9, by = 1),
                          plot = FALSE)
percentages_matched <- prop.table(hist_data_matched$counts) * 100
barplot(percentages_matched, 
        names.arg = hist_data_matched$breaks[-1], 
        col = "chartreuse4", 
        main = "NOS rating in studies with patient matching", 
        xlab = "Quality stars", 
        ylab = "Percentage of data (%)",
        xlim = c(0, 9),
        ylim = c(0, 60),
        cex.main = 2,
        cex.lab = 1.5,
        cex.axis = 1.5,
        cex.names = 1.5)

hist_data_non_matched <- hist(nos_dset[nos_dset$matched==0,]$stars, 
                              breaks = seq(3, 9, by = 1),
                              plot = FALSE)
percentages_non_matched <- prop.table(hist_data_non_matched$counts) * 100
barplot(percentages_non_matched, 
        names.arg = hist_data_non_matched$breaks[-1], 
        col = "chartreuse2", 
        main = "NOS rating in studies without patient matching", 
        xlab = "Quality stars", 
        ylab = "Percentage of data (%)",
        xlim = c(0, 9),
        ylim = c(0, 60),
        cex.main = 2,
        cex.lab = 1.5,
        cex.axis = 1.5,
        cex.names = 1.5)
dev.off()


# Multicenter or single-center studies
png("4. NOS rating histogram in multicenter or single-center studies.png", width = W_sga, height = H_sga)
par(mfrow = c(2, 1), mar = c(5, 6, 4, 0) + 0.1, mgp = c(3.25, 1, 0))
hist_data_multi <- hist(nos_dset[nos_dset$center=="multi",]$stars, 
                        breaks = seq(3, 9, by = 1),
                        plot = FALSE)
percentages_multi <- prop.table(hist_data_multi$counts) * 100
barplot(percentages_multi, 
        names.arg = hist_data_multi$breaks[-1], 
        col = "deepskyblue4", 
        main = "NOS rating in multicenter studies", 
        xlab = "Quality stars", 
        ylab = "Percentage of data (%)",
        xlim = c(0, 9),
        ylim = c(0, 60),
        cex.main = 2,
        cex.lab = 1.5,
        cex.axis = 1.5,
        cex.names = 1.5)

hist_data_single <- hist(nos_dset[nos_dset$center=="single",]$stars, 
                         breaks = seq(3, 9, by = 1),
                         plot = FALSE)
percentages_single <- prop.table(hist_data_single$counts) * 100
barplot(percentages_single, 
        names.arg = hist_data_single$breaks[-1], 
        col = "deepskyblue2", 
        main = "NOS rating in single-center studies", 
        xlab = "Quality stars", 
        ylab = "Percentage of data (%)",
        xlim = c(0, 9),
        ylim = c(0, 60),
        cex.main = 2,
        cex.lab = 1.5,
        cex.axis = 1.5,
        cex.names = 1.5)
dev.off()



