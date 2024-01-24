############################################## PACKAGES & LIBRARIES ################################################

# Package installation and library loading
install.packages("meta")
install.packages("metafor")
install.packages("metasens")

library(meta)
library(metafor)
library(metasens)


################################################### DATASETS #######################################################

# First we import the required datasets
Q_MC<-read.csv("C:/Users/sotbi/Desktop/Data/Q_MC.csv", sep=";", header=TRUE)
EBL<-read.csv("C:/Users/sotbi/Desktop/Data/EBL.csv", sep=";", header=TRUE)
OT<-read.csv("C:/Users/sotbi/Desktop/Data/OT.csv", sep=";", header=TRUE)

# Then we inspect the datasets
head(Q_MC)
head(EBL)
head(OT)

#Now we define the dataset of interest
dset<-Q_MC

# We create a new folder named "Plots_MC" on desktop 
folder_name <- "Plots_MC"
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

# And now we set the working directory to the new folder 
# This is where the plots are to be saved
setwd("C:/Users/sotbi/Desktop/Plots_MC")


################################################ DATASET FACTORS ###################################################

# First we define the order of reporting the different subgroups
order_post2018 <- c("Studies published after 2018", "Studies published before 2018")
order_matching <- c("Studies with patient matching", "Studies without patient matching")
order_center <- c("Multicenter studies", "Single-center studies")
order_robins <- c("Low Risk of Bias", "Moderate Risk of Bias", "Serious Risk of Bias")
 
# Then we apply the predefined order of reporting to the dataset
dset$post2018 <- factor(dset$post2018, levels = order_post2018)
dset$matching <- factor(dset$matching, levels = order_matching)
dset$center <- factor(dset$center, levels = order_center)
dset$robins <- factor(dset$robins, levels = order_robins)


############################################ RANDOM EFFECTS MODEL ##################################################

# First we define the width & height parameters for the different plots
# The same parameters for MRA plots are defined in the same way below
# Pooled analysis forest plot parameters
W_for_pool <- 800
H_for_pool <- 800
# Funnel & Radial plot parameters
W_fr <- 1200 
H_fr <- 780
# Subgroup analysis forest plot parameters
W_for_sga <- 800
H_for_sga <- 950

# Pooled meta analysis
png("1. Forest plot (pooled).png", width = W_for_pool, height = H_for_pool)
mc0re<-metacont(nexp,mexp,sexp,nctrl,mctrl,sctrl,data = dset, studlab = study, comb.fixed = FALSE, hakn = TRUE)
forestmc0re<-forest(mc0re, text.addline1 = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlim = c(min(round(mc0re$lower,3))-0.1, max(round(mc0re$upper,3))+0.1), showweights = TRUE, digits = 3, digits.se = 3, fs.xlab = 1, label.e = "RPN / RAPN", label.c = "OPN", label.right = "Favours OPN", label.left = "Favours RPN / RAPN", just = "center")
weights_mc <- mc0re$w.random
dev.off()

png("2. Funnel plot.png", width = W_fr, height = H_fr)
# We set the margins to add space
# mar = c(bottom, left, top, right)
# mgp1: distance between the axis title and labels, mgp2: distance between the axis ticks and labels
par(mar = c(6, 6, 5, 4), mgp = c(3.75, 1, 0))  
funnelmc0re<-funnel(mc0re, pch = 20,cex = 1, contour = c(0.9,0.95,0.99), col.contour = c("darkgray", "gray", "lightgray"), xlim = c(min(round(mc0re$TE,3))-0.1, max(round(mc0re$TE,3))+2.5), xlab = "", ylab = "", axes = FALSE)
# First we draw a rectangle by obtaining the coordinates of the plot
coords <- par("usr")
rect(coords[1], coords[3], coords[2], coords[4], border = "black", lwd = 1)
axis(2, cex.axis = 1.5, pos = coords[1])
axis(1, cex.axis = 1.5, pos = coords[3])
legend(min(round(mc0re$TE,3)), 0.8, c("0.1 > p > 0.05", "0.05 > p > 0.01", "p < 0.01"), fill = c("darkgray", "gray", "lightgray"), bty = "n", cex = 1.5)
title("Funnel plot with contours of statistical significance", cex.main = 2, xlab = "Mean Difference (MD)", ylab = "Standard Error (SE)", cex.lab = 1.5)
dev.off()

png("3. Radial plot with Egger's test.png", width = W_fr, height = H_fr)
radialmc0re<-radial(mc0re)
eggersmc0re<-metabias(mc0re, method.bias = "linreg", plotit = TRUE)
eggersmc0re
regmc0re<-lm(I(mc0re$TE/mc0re$seTE) ~ I(1/mc0re$seTE))
par(mar = c(6, 6, 5, 4), mgp = c(3.75, 1, 0)) 
radial(mc0re, cex = 1.5, cex.lab = 1.5, axes = FALSE)
coords <- par("usr")
rect(coords[1], coords[3], coords[2], coords[4], border = "black", lwd = 1)
axis(2, cex.axis = 1.5, pos = coords[1])
axis(1, cex.axis = 1.5, pos = coords[3])
abline(regmc0re)
title("Radial plot with a solid regression line for Egger's test", cex.main = 2)
dev.off()

# Small study Effects
png("4. Funnel plot with small study effects.png", width = W_fr, height = H_fr)
par(mar = c(6, 6, 5, 4), mgp = c(3.75, 1, 0)) 
l2<-limitmeta(mc0re)
funnel(l2, cex = 1.5, xlim = c(min(round(mc0re$TE,3))-0.1, max(round(mc0re$TE,3))+2.5), xlab = "", ylab = "", axes = FALSE)
coords <- par("usr")
rect(coords[1], coords[3], coords[2], coords[4], border = "black", lwd = 1)
axis(2, cex.axis = 1.5, pos = coords[1])
axis(1, cex.axis = 1.5, pos = coords[3])
print(l2, digits = 3)
title("Funnel plot with a curved regression line for small study effects", cex.main = 2, xlab = "Mean Difference (MD)", ylab = "Standard Error (SE)", cex.lab = 1.5)
dev.off()


# subgroup analysis (post2018)
png("5. Forest plot with subgroups according to publication year.png", width = W_for_sga, height = H_for_sga)
mc0repost2018<-metacont(nexp,mexp,sexp,nctrl,mctrl,sctrl,data = dset, comb.fixed = FALSE, hakn = TRUE, studlab = study, byvar = post2018, print.byvar = FALSE)
forestmc0repost2018<-forest(mc0repost2018, text.addline1 = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlim = c(min(round(mc0repost2018$lower,3))-0.1, max(round(mc0repost2018$upper,3))+0.1), showweights = TRUE, digits = 3, digits.se = 3, fs.xlab = 1, label.e = "RPN / RAPN", label.c = "OPN", label.right = "Favours OPN", label.left = "Favours RPN / RAPN", just = "center")
dev.off()
# Heterogeneity assessment
round(mc0repost2018$I2.w,3) # Expected Values of I²
round(mc0repost2018$lower.I2.w,3) # CI95% of I² (Lower Bounds) 
round(mc0repost2018$upper.I2.w,3) # CI95% of I² (Upper Bounds) 

# subgroup analysis (matching)
png("6. Forest plot with subgroups according to patient matching.png", width = W_for_sga, height = H_for_sga)
mc0rematched<-metacont(nexp,mexp,sexp,nctrl,mctrl,sctrl,data = dset, comb.fixed = FALSE, hakn = TRUE, studlab = study, byvar = matching, print.byvar = FALSE)
forestmc0rematched<-forest(mc0rematched, text.addline1 = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlim = c(min(round(mc0rematched$lower,3))-0.1, max(round(mc0rematched$upper,3))+0.1), showweights = TRUE, digits = 3, digits.se = 3, fs.xlab = 1, label.e = "RPN / RAPN", label.c = "OPN", label.right = "Favours OPN", label.left = "Favours RPN / RAPN", just = "center")
dev.off()
# Heterogeneity assessment
round(mc0rematched$I2.w,3) # Expected Values of I²
round(mc0rematched$lower.I2.w,3) # CI95% of I² (Lower Bounds) 
round(mc0rematched$upper.I2.w,3) # CI95% of I² (Upper Bounds) 

# subgroup analysis (center)
png("7. Forest plot with subgroups according to the number of centers.png", width = W_for_sga, height = H_for_sga)
mc0recenter<-metacont(nexp,mexp,sexp,nctrl,mctrl,sctrl,data = dset, comb.fixed = FALSE, hakn = TRUE, studlab = study, byvar = center, print.byvar = FALSE)
forestmc0recenter<-forest(mc0recenter, text.addline1 = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlim = c(min(round(mc0recenter$lower,3))-0.1, max(round(mc0recenter$upper,3))+0.1), showweights = TRUE, digits = 3, digits.se = 3, fs.xlab = 1, label.e = "RPN / RAPN", label.c = "OPN", label.right = "Favours OPN", label.left = "Favours RPN / RAPN", just = "center")
dev.off()
# Heterogeneity assessment
round(mc0recenter$I2.w,3) # Expected Values of I²
round(mc0recenter$lower.I2.w,3) # CI95% of I² (Lower Bounds) 
round(mc0recenter$upper.I2.w,3) # CI95% of I² (Upper Bounds)

# subgroup analysis (robins)
png("8. Forest plot with subgroups according to ROBINS-I.png", width = W_for_sga, height = H_for_sga)
mc0rerobins<-metacont(nexp,mexp,sexp,nctrl,mctrl,sctrl,data = dset, comb.fixed = FALSE, hakn = TRUE, studlab = study, byvar = robins, print.byvar = FALSE)
forestmc0rerobins<-forest(mc0rerobins, text.addline1 = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlim = c(min(round(mc0rerobins$lower,3))-0.1, max(round(mc0rerobins$upper,3))+0.1), showweights = TRUE, digits = 3, digits.se = 3, fs.xlab = 1, label.e = "RPN / RAPN", label.c = "OPN", label.right = "Favours OPN", label.left = "Favours RPN / RAPN", just = "center")
dev.off()
# Heterogeneity assessment
round(mc0rerobins$I2.w,3) # Expected Values of I²
round(mc0rerobins$lower.I2.w,3) # CI95% of I² (Lower Bounds) 
round(mc0rerobins$upper.I2.w,3) # CI95% of I² (Upper Bounds)

# We have similar patterns of MRA plots from now on, so we can adjust their width & height parameters from here
# Pooled
W_pool <- 1200
H_pool <- 780
# Matching
W_match <- 1000
H_match <- 1250
# Center
W_center <- 1000
H_center <- 1250
# Robins
W_robins <- 1000 
H_robins <- 1500

# Meta Regression Analysis for publication year in all studies (using metafor)
png("9. MRA plot according to publication year (pooled).png", width = W_pool, height = H_pool)
par(mar = c(5, 5.5, 4, 2) + 0.1)
mc0remtf<-escalc(measure = "MD",n1i = nexp, m1i = mexp, sd1i = sexp, n2i = nctrl,m2i = mctrl, sd2i = sctrl, data = dset, append = TRUE)
mc0remtfma<-rma(yi, vi, data = mc0remtf, method = "REML", mods = ~year)
yearvec<-seq(min(dset$year), max(dset$year), 1)
preds2<-predict(mc0remtfma, newmods = yearvec)
wi2<-1/sqrt(mc0remtf$vi+mc0remtfma$tau2)
size2<-1+2*(wi2 - min(wi2)) / (max(wi2) - min(wi2))
plot(mc0remtf$year, mc0remtf$yi, pch = 1, cex = size2, xlim = c(min(mc0remtf$year)-1, max(mc0remtf$year)+1), ylim = c(min(mc0remtf$yi)-1, max(mc0remtf$yi)+1), ylab = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlab = "Year of publication", main = "Meta-regression of the comparative effect in all studies", cex.main = 2, cex.lab = 1.5, xaxt = "n", yaxt = "n")
lines(yearvec, preds2$pred, lwd = 2, col = "red")
lines(yearvec, preds2$ci.lb, lty = "dashed")
lines(yearvec, preds2$ci.ub, lty = "dashed")
abline(h=0, lwd = 2, lty = "dotted")
axis(side = 2, at = seq(round(min(mc0remtf$yi)-1,0), round(max(mc0remtf$yi)+1,0), 1), cex.axis = 1.5)
axis(side = 1, cex.axis = 1.5)
dev.off()


png("10. MRA plot according to publication year in studies with or without patient matching.png", width = W_match, height = H_match)
par(mfrow = c(2, 1), mar = c(5, 5.5, 4, 2) + 0.1)
# Meta Regression Analysis for publication year in studies with patient matching (using metafor)
mc0remtfmatched<-escalc(measure = "MD",n1i = nexp, m1i = mexp, sd1i = sexp, n2i = nctrl,m2i = mctrl, sd2i = sctrl, data = dset[dset$matching=="Studies with patient matching",], append = TRUE)
mc0remtfmamatched<-rma(yi, vi, data = mc0remtfmatched, method = "REML", mods = ~year)
yearvecmatched<-seq(min(dset[dset$matching=="Studies with patient matching",]$year), max(dset[dset$matching=="Studies with patient matching",]$year), 1)
predsmatched2<-predict(mc0remtfmamatched, newmods = yearvecmatched)
wimatched2<-1/sqrt(mc0remtfmatched$vi+mc0remtfmamatched$tau2)
sizematched2<- 1+2*(wimatched2 - min(wimatched2)) / (max(wimatched2) - min(wimatched2))
plot(mc0remtfmatched$year, mc0remtfmatched$yi, pch = 1, cex = sizematched2, xlim = c(min(mc0remtf$year)-1, max(mc0remtf$year)+1), ylim = c(min(mc0remtf$yi)-1, max(mc0remtf$yi)+1), ylab = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlab = "Year of publication", main = "Meta-regression of the comparative effect in studies with patient matching", cex.main = 2, cex.lab = 1.5, xaxt = "n", yaxt = "n")
lines(yearvecmatched, predsmatched2$pred, lwd = 2, col = "darkgreen")
lines(yearvecmatched, predsmatched2$ci.lb, lty = "dashed")
lines(yearvecmatched, predsmatched2$ci.ub, lty = "dashed")
abline(h=0, lwd = 2, lty = "dotted")
axis(side = 2, at = seq(round(min(mc0remtf$yi)-1,0), round(max(mc0remtf$yi)+1,0), 1), cex.axis = 1.5)
axis(side = 1, cex.axis = 1.5)
# Meta Regression Analysis for publication year in studies without patient matching (using metafor)
mc0remtfnonmatched<-escalc(measure = "MD",n1i = nexp, m1i = mexp, sd1i = sexp, n2i = nctrl,m2i = mctrl, sd2i = sctrl, data = dset[dset$matching=="Studies without patient matching",], append = TRUE)
mc0remtfmanonmatched<-rma(yi, vi, data = mc0remtfnonmatched, method = "REML", mods = ~year)
yearvecnonmatched<-seq(min(dset[dset$matching=="Studies without patient matching",]$year), max(dset[dset$matching=="Studies without patient matching",]$year), 1)
predsnonmatched2<-predict(mc0remtfmanonmatched, newmods = yearvecnonmatched)
winonmatched2<-1/sqrt(mc0remtfnonmatched$vi+mc0remtfmanonmatched$tau2)
sizenonmatched2<- 1+2*(winonmatched2 - min(winonmatched2)) / (max(winonmatched2) - min(winonmatched2))
plot(mc0remtfnonmatched$year, mc0remtfnonmatched$yi, pch = 1, cex = sizenonmatched2, xlim = c(min(mc0remtf$year)-1, max(mc0remtf$year)+1), ylim = c(min(mc0remtf$yi)-1, max(mc0remtf$yi)+1), ylab = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlab = "Year of publication", main = "Meta-regression of the comparative effect in studies without patient matching", cex.main = 2, cex.lab = 1.5, xaxt = "n", yaxt = "n")
lines(yearvecnonmatched, predsnonmatched2$pred, lwd = 2, col = "lightgreen")
lines(yearvecnonmatched, predsnonmatched2$ci.lb, lty = "dashed")
lines(yearvecnonmatched, predsnonmatched2$ci.ub, lty = "dashed")
abline(h=0, lwd = 2, lty = "dotted")
axis(side = 2, at = seq(round(min(mc0remtf$yi)-1,0), round(max(mc0remtf$yi)+1,0), 1), cex.axis = 1.5)
axis(side = 1, cex.axis = 1.5)
dev.off()


png("11. MRA plot according to publication year in multi- or single-center studies.png", width = W_center, height = H_center)
par(mfrow = c(2, 1), mar = c(5, 5.5, 4, 2) + 0.1)
# Meta Regression Analysis for publication year in multicenter studies (using metafor)
mc0remtfmulti<-escalc(measure = "MD",n1i = nexp, m1i = mexp, sd1i = sexp, n2i = nctrl,m2i = mctrl, sd2i = sctrl, data = dset[dset$center=="Multicenter studies",], append = TRUE)
mc0remtfmamulti<-rma(yi, vi, data = mc0remtfmulti, method = "REML", mods = ~year)
yearvecmulti<-seq(min(data = dset[dset$center=="Multicenter studies",]$year), max(data = dset[dset$center=="Multicenter studies",]$year), 1)
predsmulti2<-predict(mc0remtfmamulti, newmods = yearvecmulti)
wimulti2<-1/sqrt(mc0remtfmulti$vi+mc0remtfmamulti$tau2)
sizemulti2<- 1+2*(wimulti2 - min(wimulti2)) / (max(wimulti2) - min(wimulti2))
plot(mc0remtfmulti$year, mc0remtfmulti$yi, pch = 1, cex = sizemulti2, xlim = c(min(mc0remtf$year)-1, max(mc0remtf$year)+1), ylim = c(min(mc0remtf$yi)-1, max(mc0remtf$yi)+1), ylab = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlab = "Year of publication", main = "Meta-regression of the comparative effect in multicenter studies", cex.main = 2, cex.lab = 1.5, xaxt = "n", yaxt = "n")
lines(yearvecmulti, predsmulti2$pred, lwd = 2, col = "darkblue")
lines(yearvecmulti, predsmulti2$ci.lb, lty = "dashed")
lines(yearvecmulti, predsmulti2$ci.ub, lty = "dashed")
abline(h=0, lwd = 2, lty = "dotted")
axis(side = 2, at = seq(round(min(mc0remtf$yi)-1,0), round(max(mc0remtf$yi)+1,0), 1), cex.axis = 1.5)
axis(side = 1, cex.axis = 1.5)
# Meta Regression Analysis for publication year in single-center studies (using metafor)
mc0remtfsingle<-escalc(measure = "MD",n1i = nexp, m1i = mexp, sd1i = sexp, n2i = nctrl,m2i = mctrl, sd2i = sctrl, data = dset[dset$center=="Single-center studies",], append = TRUE)
mc0remtfmasingle<-rma(yi, vi, data = mc0remtfsingle, method = "REML", mods = ~year)
yearvecsingle<-seq(min(data = dset[dset$center=="Single-center studies",]$year), max(data = dset[dset$center=="Single-center studies",]$year), 1)
predssingle2<-predict(mc0remtfmasingle, newmods = yearvecsingle)
wisingle2<-1/sqrt(mc0remtfsingle$vi+mc0remtfmasingle$tau2)
sizesingle2<- 1+2*(wisingle2 - min(wisingle2)) / (max(wisingle2) - min(wisingle2))
plot(mc0remtfsingle$year, mc0remtfsingle$yi, pch = 1, cex = sizesingle2, xlim = c(min(mc0remtf$year)-1, max(mc0remtf$year)+1), ylim = c(min(mc0remtf$yi)-1, max(mc0remtf$yi)+1), ylab = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlab = "Year of publication", main = "Meta-regression of the comparative effect in single-center studies", cex.main = 2, cex.lab = 1.5, xaxt = "n", yaxt = "n")
lines(yearvecsingle, predssingle2$pred, lwd = 2, col = "lightblue")
lines(yearvecsingle, predssingle2$ci.lb, lty = "dashed")
lines(yearvecsingle, predssingle2$ci.ub, lty = "dashed")
abline(h=0, lwd = 2, lty = "dotted")
axis(side = 2, at = seq(round(min(mc0remtf$yi)-1,0), round(max(mc0remtf$yi)+1,0), 1), cex.axis = 1.5)
axis(side = 1, cex.axis = 1.5)
dev.off()


png("12. MRA plot according to publication year in studies with low - moderate - serious risk of bias.png", width = W_robins, height = H_robins)
par(mfrow = c(3, 1), mar = c(5, 5.5, 4, 2) + 0.1)
# Meta Regression Analysis for publication year in ROBINS-I: "Low Risk of Bias" studies (using metafor)
mc0remtflow<-escalc(measure = "MD",n1i = nexp, m1i = mexp, sd1i = sexp, n2i = nctrl,m2i = mctrl, sd2i = sctrl, data = dset[dset$robins=="Low Risk of Bias",], append = TRUE)
mc0remtfmalow<-rma(yi, vi, data = mc0remtflow, method = "REML", mods = ~year)
yearveclow<-seq(min(data = dset[dset$robins=="Low Risk of Bias",]$year), max(data = dset[dset$robins=="Low Risk of Bias",]$year), 1)
predslow2<-predict(mc0remtfmalow, newmods = yearveclow)
wilow2<-1/sqrt(mc0remtflow$vi+mc0remtfmalow$tau2)
sizelow2<- 1+2*(wilow2 - min(wilow2)) / (max(wilow2) - min(wilow2))
plot(mc0remtflow$year, mc0remtflow$yi, pch = 1, cex = sizelow2, xlim = c(min(mc0remtf$year)-1, max(mc0remtf$year)+1), ylim = c(min(mc0remtf$yi)-1, max(mc0remtf$yi)+1), ylab = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlab = "Year of publication", main = "Meta-regression of the comparative effect in studies with ROBINS-I: Low", cex.main = 2.2, cex.lab = 1.7, xaxt = "n", yaxt = "n")
lines(yearveclow, predslow2$pred, lwd = 2, col = "purple")
lines(yearveclow, predslow2$ci.lb, lty = "dashed")
lines(yearveclow, predslow2$ci.ub, lty = "dashed")
abline(h=0, lwd = 2, lty = "dotted")
axis(side = 2, at = seq(round(min(mc0remtf$yi)-1,0), round(max(mc0remtf$yi)+1,0), 1), cex.axis = 1.7)
axis(side = 1, cex.axis = 1.7)
# Meta Regression Analysis for publication year in ROBINS-I: "Moderate Risk of Bias" studies (using metafor)
mc0remtfmoderate<-escalc(measure = "MD",n1i = nexp, m1i = mexp, sd1i = sexp, n2i = nctrl,m2i = mctrl, sd2i = sctrl, data = dset[dset$robins=="Moderate Risk of Bias",], append = TRUE)
mc0remtfmamoderate<-rma(yi, vi, data = mc0remtfmoderate, method = "REML", mods = ~year)
yearvecmoderate<-seq(min(data = dset[dset$robins=="Moderate Risk of Bias",]$year), max(data = dset[dset$robins=="Moderate Risk of Bias",]$year), 1)
predsmoderate2<-predict(mc0remtfmamoderate, newmods = yearvecmoderate)
wimoderate2<-1/sqrt(mc0remtfmoderate$vi+mc0remtfmamoderate$tau2)
sizemoderate2<- 1+2*(wimoderate2 - min(wimoderate2)) / (max(wimoderate2) - min(wimoderate2))
plot(mc0remtfmoderate$year, mc0remtfmoderate$yi, pch = 1, cex = sizemoderate2, xlim = c(min(mc0remtf$year)-1, max(mc0remtf$year)+1), ylim = c(min(mc0remtf$yi)-1, max(mc0remtf$yi)+1), ylab = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlab = "Year of publication", main = "Meta-regression of the comparative effect in studies with ROBINS-I: Moderate", cex.main = 2.2, cex.lab = 1.7, xaxt = "n", yaxt = "n")
lines(yearvecmoderate, predsmoderate2$pred, lwd = 2, col = "magenta")
lines(yearvecmoderate, predsmoderate2$ci.lb, lty = "dashed")
lines(yearvecmoderate, predsmoderate2$ci.ub, lty = "dashed")
abline(h=0, lwd = 2, lty = "dotted")
axis(side = 2, at = seq(round(min(mc0remtf$yi)-1,0), round(max(mc0remtf$yi)+1,0), 1), cex.axis = 1.7)
axis(side = 1, cex.axis = 1.7)
# Meta Regression Analysis for publication year in ROBINS-I: "Serious Risk of Bias" studies (using metafor)
mc0remtfserious<-escalc(measure = "MD",n1i = nexp, m1i = mexp, sd1i = sexp, n2i = nctrl,m2i = mctrl, sd2i = sctrl, data = dset[dset$robins=="Serious Risk of Bias",], append = TRUE)
mc0remtfmaserious<-rma(yi, vi, data = mc0remtfserious, method = "REML", mods = ~year)
yearvecserious<-seq(min(data = dset[dset$robins=="Serious Risk of Bias",]$year), max(data = dset[dset$robins=="Serious Risk of Bias",]$year), 1)
predsserious2<-predict(mc0remtfmaserious, newmods = yearvecserious)
wiserious2<-1/sqrt(mc0remtfserious$vi+mc0remtfmaserious$tau2)
sizeserious2<- 1+2*(wiserious2 - min(wiserious2)) / (max(wiserious2) - min(wiserious2))
plot(mc0remtfserious$year, mc0remtfserious$yi, pch = 1, cex = sizeserious2, xlim = c(min(mc0remtf$year)-1, max(mc0remtf$year)+1), ylim = c(min(mc0remtf$yi)-1, max(mc0remtf$yi)+1), ylab = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlab = "Year of publication", main = "Meta-regression of the comparative effect in studies with ROBINS-I: Serious", cex.main = 2.2, cex.lab = 1.7, xaxt = "n", yaxt = "n")
lines(yearvecserious, predsserious2$pred, lwd = 2, col = "pink")
lines(yearvecserious, predsserious2$ci.lb, lty = "dashed")
lines(yearvecserious, predsserious2$ci.ub, lty = "dashed")
abline(h=0, lwd = 2, lty = "dotted")
axis(side = 2, at = seq(round(min(mc0remtf$yi)-1,0), round(max(mc0remtf$yi)+1,0), 1), cex.axis = 1.7)
axis(side = 1, cex.axis = 1.7)
dev.off()



png("13. MRA plot according to NOS quality stars (pooled).png", width = W_pool, height = H_pool)
par(mar = c(5, 5.5, 4, 2) + 0.1)
# Meta Regression Analysis for Newcastle-Ottawa Scale (NOS) quality stars in all studies (using metafor)
mc0remtfstars<-escalc(measure = "MD",n1i = nexp, m1i = mexp, sd1i = sexp, n2i = nctrl,m2i = mctrl, sd2i = sctrl, data = dset, append = TRUE)
mc0remtfmastars<-rma(yi, vi, data = mc0remtfstars, method = "REML", mods = ~stars)
starsvec<-seq(min(dset$stars), max(dset$stars), 1)
predsstars2<-predict(mc0remtfmastars, newmods = starsvec)
wistars2<-1/sqrt(mc0remtfstars$vi+mc0remtfmastars$tau2)
sizestars2<-1+2*(wistars2 - min(wistars2)) / (max(wistars2) - min(wistars2))
plot(mc0remtfstars$stars, mc0remtfstars$yi, pch = 1, cex = sizestars2, xlim = c(min(mc0remtfstars$stars)-1, max(mc0remtfstars$stars)+1), ylim = c(min(mc0remtfstars$yi)-1, max(mc0remtfstars$yi)+1), ylab = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlab = "NOS quality stars", main = "Meta-regression of the comparative effect in all studies", cex.main = 2, cex.lab = 1.5, xaxt = "n", yaxt = "n")
lines(starsvec, predsstars2$pred, lwd = 2, col = "red")
lines(starsvec, predsstars2$ci.lb, lty = "dashed")
lines(starsvec, predsstars2$ci.ub, lty = "dashed")
abline(h=0, lwd = 2, lty = "dotted")
axis(side = 2, at = seq(round(min(mc0remtfstars$yi)-1,0), round(max(mc0remtfstars$yi)+1,0), 1), cex.axis = 1.5)
axis(side = 1, cex.axis = 1.5)
dev.off()


png("14. MRA plot according to NOS quality stars in studies with or without patient matching.png", W_match, height = H_match)
par(mfrow = c(2, 1), mar = c(5, 5.5, 4, 2) + 0.1)
# Meta Regression Analysis for Newcastle-Ottawa Scale (NOS) quality stars in studies with patient matching (using metafor)
mc0remtfstarsmatched<-escalc(measure = "MD",n1i = nexp, m1i = mexp, sd1i = sexp, n2i = nctrl,m2i = mctrl, sd2i = sctrl, data = dset[dset$matching=="Studies with patient matching",], append = TRUE)
mc0remtfmastarsmatched<-rma(yi, vi, data = mc0remtfstarsmatched, method = "REML", mods = ~stars)
starsvecmatched<-seq(min(dset[dset$matching=="Studies with patient matching",]$stars), max(dset[dset$matching=="Studies with patient matching",]$stars), 1)
predsstarsmatched2<-predict(mc0remtfmastarsmatched, newmods = starsvecmatched)
wistarsmatched2<-1/sqrt(mc0remtfstarsmatched$vi+mc0remtfmastarsmatched$tau2)
sizestarsmatched2<-1+2*(wistarsmatched2 - min(wistarsmatched2)) / (max(wistarsmatched2) - min(wistarsmatched2))
plot(mc0remtfstarsmatched$stars, mc0remtfstarsmatched$yi, pch = 1, cex = sizestarsmatched2, xlim = c(min(mc0remtfstars$stars)-1, max(mc0remtfstars$stars)+1), ylim = c(min(mc0remtfstars$yi)-1, max(mc0remtfstars$yi)+1), ylab = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlab = "NOS quality stars", main = "Meta-regression of the comparative effect in studies with patient matching", cex.main = 2, cex.lab = 1.5, xaxt = "n", yaxt = "n")
lines(starsvecmatched, predsstarsmatched2$pred, lwd = 2, col = "darkgreen")
lines(starsvecmatched, predsstarsmatched2$ci.lb, lty = "dashed")
lines(starsvecmatched, predsstarsmatched2$ci.ub, lty = "dashed")
abline(h=0, lwd = 2, lty = "dotted")
axis(side = 2, at = seq(round(min(mc0remtfstars$yi)-1,0), round(max(mc0remtfstars$yi)+1,0), 1), cex.axis = 1.5)
axis(side = 1, cex.axis = 1.5)
# Meta Regression Analysis for Newcastle-Ottawa Scale (NOS) quality stars in studies without patient matching (using metafor)
mc0remtfstarsnonmatched<-escalc(measure = "MD",n1i = nexp, m1i = mexp, sd1i = sexp, n2i = nctrl,m2i = mctrl, sd2i = sctrl, data = dset[dset$matching=="Studies without patient matching",], append = TRUE)
mc0remtfmastarsnonmatched<-rma(yi, vi, data = mc0remtfstarsnonmatched, method = "REML", mods = ~stars)
starsvecnonmatched<-seq(min(dset[dset$matching=="Studies without patient matching",]$stars), max(dset[dset$matching=="Studies without patient matching",]$stars), 1)
predsstarsnonmatched2<-predict(mc0remtfmastarsnonmatched, newmods = starsvecnonmatched)
wistarsnonmatched2<-1/sqrt(mc0remtfstarsnonmatched$vi+mc0remtfmastarsnonmatched$tau2)
sizestarsnonmatched2<-1+2*(wistarsnonmatched2 - min(wistarsnonmatched2)) / (max(wistarsnonmatched2) - min(wistarsnonmatched2))
plot(mc0remtfstarsnonmatched$stars, mc0remtfstarsnonmatched$yi, pch = 1, cex = sizestarsnonmatched2, xlim = c(min(mc0remtfstars$stars)-1, max(mc0remtfstars$stars)+1), ylim = c(min(mc0remtfstars$yi)-1, max(mc0remtfstars$yi)+1), ylab = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlab = "NOS quality stars", main = "Meta-regression of the comparative effect in studies without patient matching", cex.main = 2, cex.lab = 1.5, xaxt = "n", yaxt = "n")
lines(starsvecnonmatched, predsstarsnonmatched2$pred, lwd = 2, col = "lightgreen")
lines(starsvecnonmatched, predsstarsnonmatched2$ci.lb, lty = "dashed")
lines(starsvecnonmatched, predsstarsnonmatched2$ci.ub, lty = "dashed")
abline(h=0, lwd = 2, lty = "dotted")
axis(side = 2, at = seq(round(min(mc0remtfstars$yi)-1,0), round(max(mc0remtfstars$yi)+1,0), 1), cex.axis = 1.5)
axis(side = 1, cex.axis = 1.5)
dev.off()


png("15. MRA plot according to NOS quality stars in single- or multicenter studies.png", width = W_center, height = H_center)
par(mfrow = c(2, 1), mar = c(5, 5.5, 4, 2) + 0.1)
# Meta Regression Analysis for Newcastle-Ottawa Scale (NOS) quality stars in multicenter studies (using metafor)
mc0remtfstarsmulti<-escalc(measure = "MD",n1i = nexp, m1i = mexp, sd1i = sexp, n2i = nctrl,m2i = mctrl, sd2i = sctrl, data = dset[dset$center=="Multicenter studies",], append = TRUE)
mc0remtfmastarsmulti<-rma(yi, vi, data = mc0remtfstarsmulti, method = "REML", mods = ~stars)
starsvecmulti<-seq(min(dset[dset$center=="Multicenter studies",]$stars), max(dset[dset$center=="Multicenter studies",]$stars), 1)
predsstarsmulti2<-predict(mc0remtfmastarsmulti, newmods = starsvecmulti)
wistarsmulti2<-1/sqrt(mc0remtfstarsmulti$vi+mc0remtfmastarsmulti$tau2)
sizestarsmulti2<-1+2*(wistarsmulti2 - min(wistarsmulti2)) / (max(wistarsmulti2) - min(wistarsmulti2))
plot(mc0remtfstarsmulti$stars, mc0remtfstarsmulti$yi, pch = 1, cex = sizestarsmulti2, xlim = c(min(mc0remtfstars$stars)-1, max(mc0remtfstars$stars)+1), ylim = c(min(mc0remtfstars$yi)-1, max(mc0remtfstars$yi)+1), ylab = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlab = "NOS quality stars", main = "Meta-regression of the comparative effect in multicenter studies", cex.main = 2, cex.lab = 1.5, xaxt = "n", yaxt = "n")
lines(starsvecmulti, predsstarsmulti2$pred, lwd = 2, col = "darkblue")
lines(starsvecmulti, predsstarsmulti2$ci.lb, lty = "dashed")
lines(starsvecmulti, predsstarsmulti2$ci.ub, lty = "dashed")
abline(h=0, lwd = 2, lty = "dotted")
axis(side = 2, at = seq(round(min(mc0remtfstars$yi)-1,0), round(max(mc0remtfstars$yi)+1,0), 1), cex.axis = 1.5)
axis(side = 1, cex.axis = 1.5)
# Meta Regression Analysis for Newcastle-Ottawa Scale (NOS) quality stars in single-center studies (using metafor)
mc0remtfstarssingle<-escalc(measure = "MD",n1i = nexp, m1i = mexp, sd1i = sexp, n2i = nctrl,m2i = mctrl, sd2i = sctrl, data = dset[dset$center=="Single-center studies",], append = TRUE)
mc0remtfmastarssingle<-rma(yi, vi, data = mc0remtfstarssingle, method = "REML", mods = ~stars)
starsvecsingle<-seq(min(dset[dset$center=="Single-center studies",]$stars), max(dset[dset$center=="Single-center studies",]$stars), 1)
predsstarssingle2<-predict(mc0remtfmastarssingle, newmods = starsvecsingle)
wistarssingle2<-1/sqrt(mc0remtfstarssingle$vi+mc0remtfmastarssingle$tau2)
sizestarssingle2<-1+2*(wistarssingle2 - min(wistarssingle2)) / (max(wistarssingle2) - min(wistarssingle2))
plot(mc0remtfstarssingle$stars, mc0remtfstarssingle$yi, pch = 1, cex = sizestarssingle2, xlim = c(min(mc0remtfstars$stars)-1, max(mc0remtfstars$stars)+1), ylim = c(min(mc0remtfstars$yi)-1, max(mc0remtfstars$yi)+1), ylab = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlab = "NOS quality stars", main = "Meta-regression of the comparative effect in single-center studies", cex.main = 2, cex.lab = 1.5, xaxt = "n", yaxt = "n")
lines(starsvecsingle, predsstarssingle2$pred, lwd = 2, col = "lightblue")
lines(starsvecsingle, predsstarssingle2$ci.lb, lty = "dashed")
lines(starsvecsingle, predsstarssingle2$ci.ub, lty = "dashed")
abline(h=0, lwd = 2, lty = "dotted")
axis(side = 2, at = seq(round(min(mc0remtfstars$yi)-1,0), round(max(mc0remtfstars$yi)+1,0), 1), cex.axis = 1.5)
axis(side = 1, cex.axis = 1.5)
dev.off()



#################################### MEAN EBL & OT (RANDOM EFFECTS MODEL) ############################################

# Splitting the exp & ctrl arms
expcols<-c("study","year","country","center","matching","stars","robins","post2018","start","end","nexp","mexp","sexp")
ctrlcols<-c("study","year","country","center","matching","stars","robins","post2018","start","end","nctrl","mctrl","sctrl")


# MEAN Q in exp & ctrl
expq<-Q_MC[, expcols]
ctrlq<-Q_MC[, ctrlcols]
expqmm<-metamean(nexp,mexp,sexp,data = expq, studlab = study, comb.fixed = FALSE, hakn = TRUE)
ctrlqmm<-metamean(nctrl,mctrl,sctrl,data = ctrlq, studlab = study, comb.fixed = FALSE, hakn = TRUE)
png("16. Forest plot for Q in RPN - RAPN (pooled).png", width = W_for_pool-100, height = H_for_pool)
forestexpqmm<-forest(expqmm, xlab = "Estimated per minute blood loss in ml/min (RPN / RAPN)", xlim = c(min(round(expqmm$lower,3))-0.1, max(round(expqmm$upper,3))+0.1), showweights = TRUE, digits = 3, digits.se = 3, just = "center")
dev.off()
png("17. Forest plot for Q in OPN (pooled).png", width = W_for_pool-100, height = H_for_pool)
forestctrlqmm<-forest(ctrlqmm, xlab = "Estimated per minute blood loss in ml/min (OPN)", xlim = c(min(round(ctrlqmm$lower,3))-0.1, max(round(ctrlqmm$upper,3))+0.1), showweights = TRUE, digits = 3, digits.se = 3, just = "center")
dev.off()


# MEAN EBL in exp & ctrl
expebl<-EBL[, expcols]
ctrlebl<-EBL[, ctrlcols]
expeblmm<-metamean(nexp,mexp,sexp,data = expebl, studlab = study, comb.fixed = FALSE, hakn = TRUE)
ctrleblmm<-metamean(nctrl,mctrl,sctrl,data = ctrlebl, studlab = study, comb.fixed = FALSE, hakn = TRUE)
png("18. Forest plot for EBL in RPN - RAPN (pooled).png", width = W_for_pool-50, height = H_for_pool)
forestexpeblmm<-forest(expeblmm, xlab = "Estimated blood loss in ml (RPN / RAPN)", xlim = c(min(round(expeblmm$lower,3))-0.1, max(round(expeblmm$upper,3))+0.1), showweights = TRUE, digits = 3, digits.se = 3, just = "center")
dev.off()
png("19. Forest plot for EBL in OPN (pooled).png", width = W_for_pool-50, height = H_for_pool)
forestctrleblmm<-forest(ctrleblmm, xlab = "Estimated blood loss in ml (OPN)", xlim = c(min(round(ctrleblmm$lower,3))-0.1, max(round(ctrleblmm$upper,3))+0.1), showweights = TRUE, digits = 3, digits.se = 3, just = "center")
dev.off()


# MEAN OT in exp & ctrl
expot<-OT[, expcols]
ctrlot<-OT[, ctrlcols]
expotmm<-metamean(nexp,mexp,sexp,data = expot, studlab = study, comb.fixed = FALSE, hakn = TRUE)
ctrlotmm<-metamean(nctrl,mctrl,sctrl,data = ctrlot, studlab = study, comb.fixed = FALSE, hakn = TRUE)
png("20. Forest plot for OT in RPN - RAPN (pooled).png", width = W_for_pool-50, height = H_for_pool)
forestexpotmm<-forest(expotmm, xlab = "Operative time in min (RPN / RAPN)", xlim = c(min(round(expotmm$lower,3))-0.1, max(round(expotmm$upper,3))+0.1), showweights = TRUE, digits = 3, digits.se = 3, just = "center")
dev.off()
png("21. Forest plot for OT in OPN (pooled).png", width = W_for_pool-50, height = H_for_pool)
forestctrlotmm<-forest(ctrlotmm, xlab = "Operative time in min (OPN)", xlim = c(min(round(ctrlotmm$lower,3))-0.1, max(round(ctrlotmm$upper,3))+0.1), showweights = TRUE, digits = 3, digits.se = 3, just = "center")
dev.off()


# Now we explore the mean difference in EBL
eblre<-metacont(nexp,mexp,sexp,nctrl,mctrl,sctrl,data = EBL, studlab = study, comb.fixed = FALSE, hakn = TRUE)
png("22. Forest plot for the MD in EBL (pooled).png", width = W_for_pool+100, height = H_for_pool+20)
foresteblre<-forest(eblre, text.addline1 = "Mean difference in estimated blood loss (RPN / RAPN vs. OPN) in ml", xlim = c(min(round(eblre$lower,3))-0.1, max(round(eblre$upper,3))+0.1), showweights = TRUE, digits = 3, digits.se = 3, fs.xlab = 1, label.e = "RPN / RAPN", label.c = "OPN", label.right = "Favours OPN", label.left = "Favours RPN / RAPN", just = "center")
dev.off()


# And now we explore the mean difference in OT
otre<-metacont(nexp,mexp,sexp,nctrl,mctrl,sctrl,data = OT, studlab = study, comb.fixed = FALSE, hakn = TRUE)
png("23. Forest plot for the MD in OT (pooled).png", width = W_for_pool+100, height = H_for_pool+20)
forestotre<-forest(otre, text.addline1 = "Mean difference in operative time (RPN / RAPN vs. OPN) in min", xlim = c(min(round(otre$lower,3))-0.1, max(round(otre$upper,3))+0.1), showweights = TRUE, digits = 3, digits.se = 3, fs.xlab = 1, label.e = "RPN / RAPN", label.c = "OPN", label.right = "Favours OPN", label.left = "Favours RPN / RAPN", just = "center")
dev.off()


############################################## SA OF OUTLIERS (OSA) #################################################

# We create a new folder named "Plots_OSA" on desktop 
folder_name <- "Plots_OSA"
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

# And now we set the working directory to the new folder 
# This is where the plots are to be saved
setwd("C:/Users/sotbi/Desktop/Plots_OSA")

# We will proceed with the MA on the new dataset: dset_osa
# First we initialize the dset data.frame
dset<-Q_MC 

# We define the order of reporting the different subgroups
order_post2018 <- c("Studies published after 2018", "Studies published before 2018")
order_matching <- c("Studies with patient matching", "Studies without patient matching")
order_center <- c("Multicenter studies", "Single-center studies")
order_robins <- c("Low Risk of Bias", "Moderate Risk of Bias", "Serious Risk of Bias")

# Then we apply the predefined order of reporting to the dataset
dset$post2018 <- factor(dset$post2018, levels = order_post2018)
dset$matching <- factor(dset$matching, levels = order_matching)
dset$center <- factor(dset$center, levels = order_center)
dset$robins <- factor(dset$robins, levels = order_robins)


# We are going to identify outliers in dset based on the range between upper & lower bounds of the individual estimates
# First we compute these ranges
d_ranges <- mc0re$upper - mc0re$lower
# Then we define a reference range by adopting a precision threshold (typical values for outliers threshold: 2, 3, 3.5)
threshold <- 2
d_ref <- threshold * sd(d_ranges)
# Now we identify those outliers
ind <- which(d_ranges>d_ref)
dset_osa <- dset[-ind,]

# Now we define the width & height parameters for the new forest plots
# Pooled analysis forest plot parameters
W_for_pool <- 800
H_for_pool <- 750

# Subgroup analysis forest plot parameters
W_for_sga <- 800
H_for_sga <- 850


# Pooled meta analysis
png("1. Forest plot (pooled) - OSA.png", width = W_for_pool, height = H_for_pool)
mc0re_osa<-metacont(nexp,mexp,sexp,nctrl,mctrl,sctrl,data = dset_osa, studlab = study, comb.fixed = FALSE, hakn = TRUE)
forestmc0re_osa<-forest(mc0re_osa, text.addline1 = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlim = c(min(round(mc0re_osa$lower,3))-0.1, max(round(mc0re_osa$upper,3))+0.1), showweights = TRUE, digits = 3, digits.se = 3, fs.xlab = 1, label.e = "RPN / RAPN", label.c = "OPN", label.right = "Favours OPN", label.left = "Favours RPN / RAPN", just = "center")
dev.off()

png("2. Funnel plot - OSA.png", width = W_fr, height = H_fr)
par(mar = c(6, 6, 5, 4), mgp = c(3.75, 1, 0))
funnelmc0re_osa<-funnel(mc0re_osa, pch = 20,cex = 1, contour = c(0.9,0.95,0.99), col.contour = c("darkgray", "gray", "lightgray"), xlim = c(min(round(mc0re_osa$TE,3))-0.1, max(round(mc0re_osa$TE,3))+1.25), xlab = "", ylab = "", axes = FALSE)
coords <- par("usr")
rect(coords[1], coords[3], coords[2], coords[4], border = "black", lwd = 1)
axis(2, cex.axis = 1.5, pos = coords[1])
axis(1, cex.axis = 1.5, pos = coords[3])
legend(min(round(mc0re_osa$TE,3)), 0.025, c("0.1 > p > 0.05", "0.05 > p > 0.01", "p < 0.01"), fill = c("darkgray", "gray", "lightgray"), bty = "n", cex = 1.5)
title("Funnel plot with contours of statistical significance", cex.main = 2, xlab = "Mean Difference (MD)", ylab = "Standard Error (SE)", cex.lab = 1.5)
dev.off()

png("3. Radial plot with Egger's test - OSA.png", width = W_fr, height = H_fr)
radialmc0re_osa<-radial(mc0re_osa)
eggersmc0re_osa<-metabias(mc0re_osa, method.bias = "linreg", plotit = TRUE)
eggersmc0re_osa
regmc0re_osa<-lm(I(mc0re_osa$TE/mc0re_osa$seTE) ~ I(1/mc0re_osa$seTE))
par(mar = c(6, 6, 5, 4), mgp = c(3.75, 1, 0)) 
radial(mc0re_osa, cex = 1.5, cex.lab = 1.5, axes = FALSE)
coords <- par("usr")
rect(coords[1], coords[3], coords[2], coords[4], border = "black", lwd = 1)
axis(2, cex.axis = 1.5, pos = coords[1])
axis(1, cex.axis = 1.5, pos = coords[3])
abline(regmc0re_osa)
title("Radial plot with a solid regression line for Egger's test", cex.main = 2)
dev.off()


# Small study Effects
png("4. Funnel plot with small study effects - OSA.png", width = W_fr, height = H_fr)
par(mar = c(6, 6, 5, 4), mgp = c(3.75, 1, 0))
l2_osa<-limitmeta(mc0re_osa)
funnel(l2_osa, cex = 1.5, xlim = c(min(round(mc0re_osa$TE,3))-0.1, max(round(mc0re_osa$TE,3))+1.25), xlab = "", ylab = "", axes = FALSE)
coords <- par("usr")
rect(coords[1], coords[3], coords[2], coords[4], border = "black", lwd = 1)
axis(2, cex.axis = 1.5, pos = coords[1])
axis(1, cex.axis = 1.5, pos = coords[3])
print(l2_osa, digits = 3)
title("Funnel plot with a curved regression line for small study effects", cex.main = 2, xlab = "Mean Difference (MD)", ylab = "Standard Error (SE)", cex.lab = 1.5)
dev.off()


# subgroup analysis (post2018)
png("5. Forest plot with subgroups according to publication year - OSA.png", width = W_for_sga, height = H_for_sga)
mc0repost2018_osa<-metacont(nexp,mexp,sexp,nctrl,mctrl,sctrl,data = dset_osa, comb.fixed = FALSE, hakn = TRUE, studlab = study, byvar = post2018, print.byvar = FALSE)
forestmc0repost2018_osa<-forest(mc0repost2018_osa, text.addline1 = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlim = c(min(round(mc0repost2018_osa$lower,3))-0.1, max(round(mc0repost2018_osa$upper,3))+0.1), showweights = TRUE, digits = 3, digits.se = 3, fs.xlab = 1, label.e = "RPN / RAPN", label.c = "OPN", label.right = "Favours OPN", label.left = "Favours RPN / RAPN", just = "center")
dev.off()
# Heterogeneity assessment
round(mc0repost2018_osa$I2.w,3) # Expected Values of I²
round(mc0repost2018_osa$lower.I2.w,3) # CI95% of I² (Lower Bounds) 
round(mc0repost2018_osa$upper.I2.w,3) # CI95% of I² (Upper Bounds) 

# subgroup analysis (matching)
png("6. Forest plot with subgroups according to patient matching - OSA.png", width = W_for_sga, height = H_for_sga)
mc0rematched_osa<-metacont(nexp,mexp,sexp,nctrl,mctrl,sctrl,data = dset_osa, comb.fixed = FALSE, hakn = TRUE, studlab = study, byvar = matching, print.byvar = FALSE)
forestmc0rematched_osa<-forest(mc0rematched_osa, text.addline1 = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlim = c(min(round(mc0rematched_osa$lower,3))-0.1, max(round(mc0rematched_osa$upper,3))+0.1), showweights = TRUE, digits = 3, digits.se = 3, fs.xlab = 1, label.e = "RPN / RAPN", label.c = "OPN", label.right = "Favours OPN", label.left = "Favours RPN / RAPN", just = "center")
dev.off()
# Heterogeneity assessment
round(mc0rematched_osa$I2.w,3) # Expected Values of I²
round(mc0rematched_osa$lower.I2.w,3) # CI95% of I² (Lower Bounds) 
round(mc0rematched_osa$upper.I2.w,3) # CI95% of I² (Upper Bounds) 

# subgroup analysis (center)
png("7. Forest plot with subgroups according to the number of centers - OSA.png", width = W_for_sga, height = H_for_sga)
mc0recenter_osa<-metacont(nexp,mexp,sexp,nctrl,mctrl,sctrl,data = dset_osa, comb.fixed = FALSE, hakn = TRUE, studlab = study, byvar = center, print.byvar = FALSE)
forestmc0recenter_osa<-forest(mc0recenter_osa, text.addline1 = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlim = c(min(round(mc0recenter_osa$lower,3))-0.1, max(round(mc0recenter_osa$upper,3))+0.1), showweights = TRUE, digits = 3, digits.se = 3, fs.xlab = 1, label.e = "RPN / RAPN", label.c = "OPN", label.right = "Favours OPN", label.left = "Favours RPN / RAPN", just = "center")
dev.off()
# Heterogeneity assessment
round(mc0recenter_osa$I2.w,3) # Expected Values of I²
round(mc0recenter_osa$lower.I2.w,3) # CI95% of I² (Lower Bounds) 
round(mc0recenter_osa$upper.I2.w,3) # CI95% of I² (Upper Bounds) 

# subgroup analysis (robins)
png("8. Forest plot with subgroups according to ROBINS-I - OSA.png", width = W_for_sga, height = H_for_sga)
mc0rerobins_osa<-metacont(nexp,mexp,sexp,nctrl,mctrl,sctrl,data = dset_osa, comb.fixed = FALSE, hakn = TRUE, studlab = study, byvar = robins, print.byvar = FALSE)
forestmc0rerobins_osa<-forest(mc0rerobins_osa, text.addline1 = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlim = c(min(round(mc0rerobins_osa$lower,3))-0.1, max(round(mc0rerobins_osa$upper,3))+0.1), showweights = TRUE, digits = 3, digits.se = 3, fs.xlab = 1, label.e = "RPN / RAPN", label.c = "OPN", label.right = "Favours OPN", label.left = "Favours RPN / RAPN", just = "center")
dev.off()
# Heterogeneity assessment
round(mc0rerobins_osa$I2.w,3) # Expected Values of I²
round(mc0rerobins_osa$lower.I2.w,3) # CI95% of I² (Lower Bounds) 
round(mc0rerobins_osa$upper.I2.w,3) # CI95% of I² (Upper Bounds)


# Now we will proceed with the MRA on the new dataset: dset_osa

# Meta Regression Analysis for publication year in all studies (using metafor)
png("9. MRA plot according to publication year (pooled) - OSA.png", width = W_pool, height = H_pool)
par(mar = c(5, 5.5, 4, 2) + 0.1)
mc0remtf_osa<-escalc(measure = "MD",n1i = nexp, m1i = mexp, sd1i = sexp, n2i = nctrl,m2i = mctrl, sd2i = sctrl, data = dset_osa, append = TRUE)
mc0remtfma_osa<-rma(yi, vi, data = mc0remtf_osa, method = "REML", mods = ~year)
yearvec_osa<-seq(min(dset_osa$year), max(dset_osa$year), 1)
preds2_osa<-predict(mc0remtfma_osa, newmods = yearvec_osa)
wi2_osa<-1/sqrt(mc0remtf_osa$vi+mc0remtfma_osa$tau2)
size2_osa<-1+2*(wi2_osa - min(wi2_osa)) / (max(wi2_osa) - min(wi2_osa))
plot(mc0remtf_osa$year, mc0remtf_osa$yi, pch = 1, cex = size2_osa, xlim = c(min(mc0remtf_osa$year)-1, max(mc0remtf_osa$year)+1), ylim = c(min(mc0remtf_osa$yi)-1, max(mc0remtf_osa$yi)+1), ylab = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlab = "Year of publication", main = "Meta-regression of the comparative effect in all studies", cex.main = 2, cex.lab = 1.5, xaxt = "n", yaxt = "n")
lines(yearvec_osa, preds2_osa$pred, lwd = 2, col = "red")
lines(yearvec_osa, preds2_osa$ci.lb, lty = "dashed")
lines(yearvec_osa, preds2_osa$ci.ub, lty = "dashed")
abline(h=0, lwd = 2, lty = "dotted")
axis(side = 2, at = seq(round(min(mc0remtf_osa$yi)-1,0), round(max(mc0remtf_osa$yi)+1,0), 1), cex.axis = 1.5)
axis(side = 1, cex.axis = 1.5)
dev.off()


png("10. MRA plot according to publication year in studies with or without patient matching - OSA.png", width = W_match, height = H_match)
par(mfrow = c(2, 1), mar = c(5, 5.5, 4, 2) + 0.1)
# Meta Regression Analysis for publication year in studies with patient matching (using metafor)
mc0remtfmatched_osa<-escalc(measure = "MD",n1i = nexp, m1i = mexp, sd1i = sexp, n2i = nctrl,m2i = mctrl, sd2i = sctrl, data = dset_osa[dset_osa$matching=="Studies with patient matching",], append = TRUE)
mc0remtfmamatched_osa<-rma(yi, vi, data = mc0remtfmatched_osa, method = "REML", mods = ~year)
yearvecmatched_osa<-seq(min(dset_osa[dset_osa$matching=="Studies with patient matching",]$year), max(dset_osa[dset_osa$matching=="Studies with patient matching",]$year), 1)
predsmatched2_osa<-predict(mc0remtfmamatched_osa, newmods = yearvecmatched_osa)
wimatched2_osa<-1/sqrt(mc0remtfmatched_osa$vi+mc0remtfmamatched_osa$tau2)
sizematched2_osa<- 1+2*(wimatched2_osa - min(wimatched2_osa)) / (max(wimatched2_osa) - min(wimatched2_osa))
plot(mc0remtfmatched_osa$year, mc0remtfmatched_osa$yi, pch = 1, cex = sizematched2_osa, xlim = c(min(mc0remtf_osa$year)-1, max(mc0remtf_osa$year)+1), ylim = c(min(mc0remtf_osa$yi)-1, max(mc0remtf_osa$yi)+1), ylab = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlab = "Year of publication", main = "Meta-regression of the comparative effect in studies with patient matching", cex.main = 2, cex.lab = 1.5, xaxt = "n", yaxt = "n")
lines(yearvecmatched_osa, predsmatched2_osa$pred, lwd = 2, col = "darkgreen")
lines(yearvecmatched_osa, predsmatched2_osa$ci.lb, lty = "dashed")
lines(yearvecmatched_osa, predsmatched2_osa$ci.ub, lty = "dashed")
abline(h=0, lwd = 2, lty = "dotted")
axis(side = 2, at = seq(round(min(mc0remtf_osa$yi)-1,0), round(max(mc0remtf_osa$yi)+1,0), 1), cex.axis = 1.5)
axis(side = 1, cex.axis = 1.5)
# Meta Regression Analysis for publication year in studies without patient matching (using metafor)
mc0remtfnonmatched_osa<-escalc(measure = "MD",n1i = nexp, m1i = mexp, sd1i = sexp, n2i = nctrl,m2i = mctrl, sd2i = sctrl, data = dset_osa[dset_osa$matching=="Studies without patient matching",], append = TRUE)
mc0remtfmanonmatched_osa<-rma(yi, vi, data = mc0remtfnonmatched_osa, method = "REML", mods = ~year)
yearvecnonmatched_osa<-seq(min(dset_osa[dset_osa$matching=="Studies without patient matching",]$year), max(dset_osa[dset_osa$matching=="Studies without patient matching",]$year), 1)
predsnonmatched2_osa<-predict(mc0remtfmanonmatched_osa, newmods = yearvecnonmatched_osa)
winonmatched2_osa<-1/sqrt(mc0remtfnonmatched_osa$vi+mc0remtfmanonmatched_osa$tau2)
sizenonmatched2_osa<- 1+2*(winonmatched2_osa - min(winonmatched2_osa)) / (max(winonmatched2_osa) - min(winonmatched2_osa))
plot(mc0remtfnonmatched_osa$year, mc0remtfnonmatched_osa$yi, pch = 1, cex = sizenonmatched2_osa, xlim = c(min(mc0remtf_osa$year)-1, max(mc0remtf_osa$year)+1), ylim = c(min(mc0remtf_osa$yi)-1, max(mc0remtf_osa$yi)+1), ylab = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlab = "Year of publication", main = "Meta-regression of the comparative effect in studies without patient matching", cex.main = 2, cex.lab = 1.5, xaxt = "n", yaxt = "n")
lines(yearvecnonmatched_osa, predsnonmatched2_osa$pred, lwd = 2, col = "lightgreen")
lines(yearvecnonmatched_osa, predsnonmatched2_osa$ci.lb, lty = "dashed")
lines(yearvecnonmatched_osa, predsnonmatched2_osa$ci.ub, lty = "dashed")
abline(h=0, lwd = 2, lty = "dotted")
axis(side = 2, at = seq(round(min(mc0remtf_osa$yi)-1,0), round(max(mc0remtf_osa$yi)+1,0), 1), cex.axis = 1.5)
axis(side = 1, cex.axis = 1.5)
dev.off()


png("11. MRA plot according to publication year in multi- or single-center studies - OSA.png", width = W_center, height = H_center)
par(mfrow = c(2, 1), mar = c(5, 5.5, 4, 2) + 0.1)
# Meta Regression Analysis for publication year in multicenter studies (using metafor)
mc0remtfmulti_osa<-escalc(measure = "MD",n1i = nexp, m1i = mexp, sd1i = sexp, n2i = nctrl,m2i = mctrl, sd2i = sctrl, data = dset_osa[dset_osa$center=="Multicenter studies",], append = TRUE)
mc0remtfmamulti_osa<-rma(yi, vi, data = mc0remtfmulti_osa, method = "REML", mods = ~year)
yearvecmulti_osa<-seq(min(data = dset_osa[dset_osa$center=="Multicenter studies",]$year), max(data = dset_osa[dset_osa$center=="Multicenter studies",]$year), 1)
predsmulti2_osa<-predict(mc0remtfmamulti_osa, newmods = yearvecmulti_osa)
wimulti2_osa<-1/sqrt(mc0remtfmulti_osa$vi+mc0remtfmamulti_osa$tau2)
sizemulti2_osa<- 1+2*(wimulti2_osa - min(wimulti2_osa)) / (max(wimulti2_osa) - min(wimulti2_osa))
plot(mc0remtfmulti_osa$year, mc0remtfmulti_osa$yi, pch = 1, cex = sizemulti2_osa, xlim = c(min(mc0remtf_osa$year)-1, max(mc0remtf_osa$year)+1), ylim = c(min(mc0remtf_osa$yi)-1, max(mc0remtf_osa$yi)+1), ylab = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlab = "Year of publication", main = "Meta-regression of the comparative effect in multicenter studies", cex.main = 2, cex.lab = 1.5, xaxt = "n", yaxt = "n")
lines(yearvecmulti_osa, predsmulti2_osa$pred, lwd = 2, col = "darkblue")
lines(yearvecmulti_osa, predsmulti2_osa$ci.lb, lty = "dashed")
lines(yearvecmulti_osa, predsmulti2_osa$ci.ub, lty = "dashed")
abline(h=0, lwd = 2, lty = "dotted")
axis(side = 2, at = seq(round(min(mc0remtf_osa$yi)-1,0), round(max(mc0remtf_osa$yi)+1,0), 1), cex.axis = 1.5)
axis(side = 1, cex.axis = 1.5)
# Meta Regression Analysis for publication year in single-center studies (using metafor)
mc0remtfsingle_osa<-escalc(measure = "MD",n1i = nexp, m1i = mexp, sd1i = sexp, n2i = nctrl,m2i = mctrl, sd2i = sctrl, data = dset_osa[dset_osa$center=="Single-center studies",], append = TRUE)
mc0remtfmasingle_osa<-rma(yi, vi, data = mc0remtfsingle_osa, method = "REML", mods = ~year)
yearvecsingle_osa<-seq(min(data = dset_osa[dset_osa$center=="Single-center studies",]$year), max(data = dset_osa[dset_osa$center=="Single-center studies",]$year), 1)
predssingle2_osa<-predict(mc0remtfmasingle_osa, newmods = yearvecsingle_osa)
wisingle2_osa<-1/sqrt(mc0remtfsingle_osa$vi+mc0remtfmasingle_osa$tau2)
sizesingle2_osa<- 1+2*(wisingle2_osa - min(wisingle2_osa)) / (max(wisingle2_osa) - min(wisingle2_osa))
plot(mc0remtfsingle_osa$year, mc0remtfsingle_osa$yi, pch = 1, cex = sizesingle2_osa, xlim = c(min(mc0remtf_osa$year)-1, max(mc0remtf_osa$year)+1), ylim = c(min(mc0remtf_osa$yi)-1, max(mc0remtf_osa$yi)+1), ylab = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlab = "Year of publication", main = "Meta-regression of the comparative effect in single-center studies", cex.main = 2, cex.lab = 1.5, xaxt = "n", yaxt = "n")
lines(yearvecsingle_osa, predssingle2_osa$pred, lwd = 2, col = "lightblue")
lines(yearvecsingle_osa, predssingle2_osa$ci.lb, lty = "dashed")
lines(yearvecsingle_osa, predssingle2_osa$ci.ub, lty = "dashed")
abline(h=0, lwd = 2, lty = "dotted")
axis(side = 2, at = seq(round(min(mc0remtf_osa$yi)-1,0), round(max(mc0remtf_osa$yi)+1,0), 1), cex.axis = 1.5)
axis(side = 1, cex.axis = 1.5)
dev.off()


png("12. MRA plot according to publication year in studies with low - moderate - serious risk of bias - OSA.png", width = W_robins, height = H_robins)
par(mfrow = c(3, 1), mar = c(5, 5.5, 4, 2) + 0.1)
# Meta Regression Analysis for publication year in ROBINS-I: "Low Risk of Bias" studies (using metafor)
mc0remtflow_osa<-escalc(measure = "MD",n1i = nexp, m1i = mexp, sd1i = sexp, n2i = nctrl,m2i = mctrl, sd2i = sctrl, data = dset_osa[dset_osa$robins=="Low Risk of Bias",], append = TRUE)
mc0remtfmalow_osa<-rma(yi, vi, data = mc0remtflow_osa, method = "REML", mods = ~year)
yearveclow_osa<-seq(min(data = dset_osa[dset_osa$robins=="Low Risk of Bias",]$year), max(data = dset_osa[dset_osa$robins=="Low Risk of Bias",]$year), 1)
predslow2_osa<-predict(mc0remtfmalow_osa, newmods = yearveclow_osa)
wilow2_osa<-1/sqrt(mc0remtflow_osa$vi+mc0remtfmalow_osa$tau2)
sizelow2_osa<- 1+2*(wilow2_osa - min(wilow2_osa)) / (max(wilow2_osa) - min(wilow2_osa))
plot(mc0remtflow_osa$year, mc0remtflow_osa$yi, pch = 1, cex = sizelow2_osa, xlim = c(min(mc0remtf_osa$year)-1, max(mc0remtf_osa$year)+1), ylim = c(min(mc0remtf_osa$yi)-1, max(mc0remtf_osa$yi)+1), ylab = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlab = "Year of publication", main = "Meta-regression of the comparative effect in studies with ROBINS-I: Low", cex.main = 2.2, cex.lab = 1.7, xaxt = "n", yaxt = "n")
lines(yearveclow_osa, predslow2_osa$pred, lwd = 2, col = "purple")
lines(yearveclow_osa, predslow2_osa$ci.lb, lty = "dashed")
lines(yearveclow_osa, predslow2_osa$ci.ub, lty = "dashed")
abline(h=0, lwd = 2, lty = "dotted")
axis(side = 2, at = seq(round(min(mc0remtf_osa$yi)-1,0), round(max(mc0remtf_osa$yi)+1,0), 1), cex.axis = 1.7)
axis(side = 1, cex.axis = 1.7)
# Meta Regression Analysis for publication year in ROBINS-I: "Moderate Risk of Bias" studies (using metafor)
mc0remtfmoderate_osa<-escalc(measure = "MD",n1i = nexp, m1i = mexp, sd1i = sexp, n2i = nctrl,m2i = mctrl, sd2i = sctrl, data = dset_osa[dset_osa$robins=="Moderate Risk of Bias",], append = TRUE)
mc0remtfmamoderate_osa<-rma(yi, vi, data = mc0remtfmoderate_osa, method = "REML", mods = ~year)
yearvecmoderate_osa<-seq(min(data = dset_osa[dset_osa$robins=="Moderate Risk of Bias",]$year), max(data = dset_osa[dset_osa$robins=="Moderate Risk of Bias",]$year), 1)
predsmoderate2_osa<-predict(mc0remtfmamoderate_osa, newmods = yearvecmoderate_osa)
wimoderate2_osa<-1/sqrt(mc0remtfmoderate_osa$vi+mc0remtfmamoderate_osa$tau2)
sizemoderate2_osa<- 1+2*(wimoderate2_osa - min(wimoderate2_osa)) / (max(wimoderate2_osa) - min(wimoderate2_osa))
plot(mc0remtfmoderate_osa$year, mc0remtfmoderate_osa$yi, pch = 1, cex = sizemoderate2_osa, xlim = c(min(mc0remtf_osa$year)-1, max(mc0remtf_osa$year)+1), ylim = c(min(mc0remtf_osa$yi)-1, max(mc0remtf_osa$yi)+1), ylab = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlab = "Year of publication", main = "Meta-regression of the comparative effect in studies with ROBINS-I: Moderate", cex.main = 2.2, cex.lab = 1.7, xaxt = "n", yaxt = "n")
lines(yearvecmoderate_osa, predsmoderate2_osa$pred, lwd = 2, col = "magenta")
lines(yearvecmoderate_osa, predsmoderate2_osa$ci.lb, lty = "dashed")
lines(yearvecmoderate_osa, predsmoderate2_osa$ci.ub, lty = "dashed")
abline(h=0, lwd = 2, lty = "dotted")
axis(side = 2, at = seq(round(min(mc0remtf_osa$yi)-1,0), round(max(mc0remtf_osa$yi)+1,0), 1), cex.axis = 1.7)
axis(side = 1, cex.axis = 1.7)
# Meta Regression Analysis for publication year in ROBINS-I: "Serious Risk of Bias" studies (using metafor)
mc0remtfserious_osa<-escalc(measure = "MD",n1i = nexp, m1i = mexp, sd1i = sexp, n2i = nctrl,m2i = mctrl, sd2i = sctrl, data = dset_osa[dset_osa$robins=="Serious Risk of Bias",], append = TRUE)
mc0remtfmaserious_osa<-rma(yi, vi, data = mc0remtfserious_osa, method = "REML", mods = ~year)
yearvecserious_osa<-seq(min(data = dset_osa[dset_osa$robins=="Serious Risk of Bias",]$year), max(data = dset_osa[dset_osa$robins=="Serious Risk of Bias",]$year), 1)
predsserious2_osa<-predict(mc0remtfmaserious_osa, newmods = yearvecserious_osa)
wiserious2_osa<-1/sqrt(mc0remtfserious_osa$vi+mc0remtfmaserious_osa$tau2)
sizeserious2_osa<- 1+2*(wiserious2_osa - min(wiserious2_osa)) / (max(wiserious2_osa) - min(wiserious2_osa))
plot(mc0remtfserious_osa$year, mc0remtfserious_osa$yi, pch = 1, cex = sizeserious2_osa, xlim = c(min(mc0remtf_osa$year)-1, max(mc0remtf_osa$year)+1), ylim = c(min(mc0remtf_osa$yi)-1, max(mc0remtf_osa$yi)+1), ylab = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlab = "Year of publication", main = "Meta-regression of the comparative effect in studies with ROBINS-I: Serious", cex.main = 2.2, cex.lab = 1.7, xaxt = "n", yaxt = "n")
lines(yearvecserious_osa, predsserious2_osa$pred, lwd = 2, col = "pink")
lines(yearvecserious_osa, predsserious2_osa$ci.lb, lty = "dashed")
lines(yearvecserious_osa, predsserious2_osa$ci.ub, lty = "dashed")
abline(h=0, lwd = 2, lty = "dotted")
axis(side = 2, at = seq(round(min(mc0remtf_osa$yi)-1,0), round(max(mc0remtf_osa$yi)+1,0), 1), cex.axis = 1.7)
axis(side = 1, cex.axis = 1.7)
dev.off()



png("13. MRA plot according to NOS quality stars (pooled) - OSA.png", width = W_pool, height = H_pool)
par(mar = c(5, 5.5, 4, 2) + 0.1)
# Meta Regression Analysis for Newcastle-Ottawa Scale (NOS) quality stars in all studies (using metafor)
mc0remtfstars_osa<-escalc(measure = "MD",n1i = nexp, m1i = mexp, sd1i = sexp, n2i = nctrl,m2i = mctrl, sd2i = sctrl, data = dset_osa, append = TRUE)
mc0remtfmastars_osa<-rma(yi, vi, data = mc0remtfstars_osa, method = "REML", mods = ~stars)
starsvec_osa<-seq(min(dset_osa$stars), max(dset_osa$stars), 1)
predsstars2_osa<-predict(mc0remtfmastars_osa, newmods = starsvec_osa)
wistars2_osa<-1/sqrt(mc0remtfstars_osa$vi+mc0remtfmastars_osa$tau2)
sizestars2_osa<-1+2*(wistars2_osa - min(wistars2_osa)) / (max(wistars2_osa) - min(wistars2_osa))
plot(mc0remtfstars_osa$stars, mc0remtfstars_osa$yi, pch = 1, cex = sizestars2_osa, xlim = c(min(mc0remtfstars_osa$stars)-1, max(mc0remtfstars_osa$stars)+1), ylim = c(min(mc0remtfstars_osa$yi)-1, max(mc0remtfstars_osa$yi)+1), ylab = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlab = "NOS quality stars", main = "Meta-regression of the comparative effect in all studies", cex.main = 2, cex.lab = 1.5, xaxt = "n", yaxt = "n")
lines(starsvec_osa, predsstars2_osa$pred, lwd = 2, col = "red")
lines(starsvec_osa, predsstars2_osa$ci.lb, lty = "dashed")
lines(starsvec_osa, predsstars2_osa$ci.ub, lty = "dashed")
abline(h=0, lwd = 2, lty = "dotted")
axis(side = 2, at = seq(round(min(mc0remtfstars_osa$yi)-1,0), round(max(mc0remtfstars_osa$yi)+1,0), 1), cex.axis = 1.5)
axis(side = 1, cex.axis = 1.5)
dev.off()


png("14. MRA plot according to NOS quality stars in studies with or without patient matching - OSA.png", W_match, height = H_match)
par(mfrow = c(2, 1), mar = c(5, 5.5, 4, 2) + 0.1)
# Meta Regression Analysis for Newcastle-Ottawa Scale (NOS) quality stars in studies with patient matching (using metafor)
mc0remtfstarsmatched_osa<-escalc(measure = "MD",n1i = nexp, m1i = mexp, sd1i = sexp, n2i = nctrl,m2i = mctrl, sd2i = sctrl, data = dset_osa[dset_osa$matching=="Studies with patient matching",], append = TRUE)
mc0remtfmastarsmatched_osa<-rma(yi, vi, data = mc0remtfstarsmatched_osa, method = "REML", mods = ~stars)
starsvecmatched_osa<-seq(min(dset_osa[dset_osa$matching=="Studies with patient matching",]$stars), max(dset_osa[dset_osa$matching=="Studies with patient matching",]$stars), 1)
predsstarsmatched2_osa<-predict(mc0remtfmastarsmatched_osa, newmods = starsvecmatched_osa)
wistarsmatched2_osa<-1/sqrt(mc0remtfstarsmatched_osa$vi+mc0remtfmastarsmatched_osa$tau2)
sizestarsmatched2_osa<-1+2*(wistarsmatched2_osa - min(wistarsmatched2_osa)) / (max(wistarsmatched2_osa) - min(wistarsmatched2_osa))
plot(mc0remtfstarsmatched_osa$stars, mc0remtfstarsmatched_osa$yi, pch = 1, cex = sizestarsmatched2_osa, xlim = c(min(mc0remtfstars_osa$stars)-1, max(mc0remtfstars_osa$stars)+1), ylim = c(min(mc0remtfstars_osa$yi)-1, max(mc0remtfstars_osa$yi)+1), ylab = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlab = "NOS quality stars", main = "Meta-regression of the comparative effect in studies with patient matching", cex.main = 2, cex.lab = 1.5, xaxt = "n", yaxt = "n")
lines(starsvecmatched_osa, predsstarsmatched2_osa$pred, lwd = 2, col = "darkgreen")
lines(starsvecmatched_osa, predsstarsmatched2_osa$ci.lb, lty = "dashed")
lines(starsvecmatched_osa, predsstarsmatched2_osa$ci.ub, lty = "dashed")
abline(h=0, lwd = 2, lty = "dotted")
axis(side = 2, at = seq(round(min(mc0remtfstars_osa$yi)-1,0), round(max(mc0remtfstars_osa$yi)+1,0), 1), cex.axis = 1.5)
axis(side = 1, cex.axis = 1.5)
# Meta Regression Analysis for Newcastle-Ottawa Scale (NOS) quality stars in studies without patient matching (using metafor)
mc0remtfstarsnonmatched_osa<-escalc(measure = "MD",n1i = nexp, m1i = mexp, sd1i = sexp, n2i = nctrl,m2i = mctrl, sd2i = sctrl, data = dset_osa[dset_osa$matching=="Studies without patient matching",], append = TRUE)
mc0remtfmastarsnonmatched_osa<-rma(yi, vi, data = mc0remtfstarsnonmatched_osa, method = "REML", mods = ~stars)
starsvecnonmatched_osa<-seq(min(dset_osa[dset_osa$matching=="Studies without patient matching",]$stars), max(dset_osa[dset_osa$matching=="Studies without patient matching",]$stars), 1)
predsstarsnonmatched2_osa<-predict(mc0remtfmastarsnonmatched_osa, newmods = starsvecnonmatched_osa)
wistarsnonmatched2_osa<-1/sqrt(mc0remtfstarsnonmatched_osa$vi+mc0remtfmastarsnonmatched_osa$tau2)
sizestarsnonmatched2_osa<-1+2*(wistarsnonmatched2_osa - min(wistarsnonmatched2_osa)) / (max(wistarsnonmatched2_osa) - min(wistarsnonmatched2_osa))
plot(mc0remtfstarsnonmatched_osa$stars, mc0remtfstarsnonmatched_osa$yi, pch = 1, cex = sizestarsnonmatched2_osa, xlim = c(min(mc0remtfstars_osa$stars)-1, max(mc0remtfstars_osa$stars)+1), ylim = c(min(mc0remtfstars_osa$yi)-1, max(mc0remtfstars_osa$yi)+1), ylab = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlab = "NOS quality stars", main = "Meta-regression of the comparative effect in studies without patient matching", cex.main = 2, cex.lab = 1.5, xaxt = "n", yaxt = "n")
lines(starsvecnonmatched_osa, predsstarsnonmatched2_osa$pred, lwd = 2, col = "lightgreen")
lines(starsvecnonmatched_osa, predsstarsnonmatched2_osa$ci.lb, lty = "dashed")
lines(starsvecnonmatched_osa, predsstarsnonmatched2_osa$ci.ub, lty = "dashed")
abline(h=0, lwd = 2, lty = "dotted")
axis(side = 2, at = seq(round(min(mc0remtfstars_osa$yi)-1,0), round(max(mc0remtfstars_osa$yi)+1,0), 1), cex.axis = 1.5)
axis(side = 1, cex.axis = 1.5)
dev.off()


png("15. MRA plot according to NOS quality stars in single- or multicenter studies - OSA.png", width = W_center, height = H_center)
par(mfrow = c(2, 1), mar = c(5, 5.5, 4, 2) + 0.1)
# Meta Regression Analysis for Newcastle-Ottawa Scale (NOS) quality stars in multicenter studies (using metafor)
mc0remtfstarsmulti_osa<-escalc(measure = "MD",n1i = nexp, m1i = mexp, sd1i = sexp, n2i = nctrl,m2i = mctrl, sd2i = sctrl, data = dset_osa[dset_osa$center=="Multicenter studies",], append = TRUE)
mc0remtfmastarsmulti_osa<-rma(yi, vi, data = mc0remtfstarsmulti_osa, method = "REML", mods = ~stars)
starsvecmulti_osa<-seq(min(dset_osa[dset_osa$center=="Multicenter studies",]$stars), max(dset_osa[dset_osa$center=="Multicenter studies",]$stars), 1)
predsstarsmulti2_osa<-predict(mc0remtfmastarsmulti_osa, newmods = starsvecmulti_osa)
wistarsmulti2_osa<-1/sqrt(mc0remtfstarsmulti_osa$vi+mc0remtfmastarsmulti_osa$tau2)
sizestarsmulti2_osa<-1+2*(wistarsmulti2_osa - min(wistarsmulti2_osa)) / (max(wistarsmulti2_osa) - min(wistarsmulti2_osa))
plot(mc0remtfstarsmulti_osa$stars, mc0remtfstarsmulti_osa$yi, pch = 1, cex = sizestarsmulti2_osa, xlim = c(min(mc0remtfstars_osa$stars)-1, max(mc0remtfstars_osa$stars)+1), ylim = c(min(mc0remtfstars_osa$yi)-1, max(mc0remtfstars_osa$yi)+1), ylab = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlab = "NOS quality stars", main = "Meta-regression of the comparative effect in multicenter studies", cex.main = 2, cex.lab = 1.5, xaxt = "n", yaxt = "n")
lines(starsvecmulti_osa, predsstarsmulti2_osa$pred, lwd = 2, col = "darkblue")
lines(starsvecmulti_osa, predsstarsmulti2_osa$ci.lb, lty = "dashed")
lines(starsvecmulti_osa, predsstarsmulti2_osa$ci.ub, lty = "dashed")
abline(h=0, lwd = 2, lty = "dotted")
axis(side = 2, at = seq(round(min(mc0remtfstars_osa$yi)-1,0), round(max(mc0remtfstars_osa$yi)+1,0), 1), cex.axis = 1.5)
axis(side = 1, cex.axis = 1.5)
# Meta Regression Analysis for Newcastle-Ottawa Scale (NOS) quality stars in single-center studies (using metafor)
mc0remtfstarssingle_osa<-escalc(measure = "MD",n1i = nexp, m1i = mexp, sd1i = sexp, n2i = nctrl,m2i = mctrl, sd2i = sctrl, data = dset_osa[dset_osa$center=="Single-center studies",], append = TRUE)
mc0remtfmastarssingle_osa<-rma(yi, vi, data = mc0remtfstarssingle_osa, method = "REML", mods = ~stars)
starsvecsingle_osa<-seq(min(dset_osa[dset_osa$center=="Single-center studies",]$stars), max(dset_osa[dset_osa$center=="Single-center studies",]$stars), 1)
predsstarssingle2_osa<-predict(mc0remtfmastarssingle_osa, newmods = starsvecsingle_osa)
wistarssingle2_osa<-1/sqrt(mc0remtfstarssingle_osa$vi+mc0remtfmastarssingle_osa$tau2)
sizestarssingle2_osa<-1+2*(wistarssingle2_osa - min(wistarssingle2_osa)) / (max(wistarssingle2_osa) - min(wistarssingle2_osa))
plot(mc0remtfstarssingle_osa$stars, mc0remtfstarssingle_osa$yi, pch = 1, cex = sizestarssingle2_osa, xlim = c(min(mc0remtfstars_osa$stars)-1, max(mc0remtfstars_osa$stars)+1), ylim = c(min(mc0remtfstars_osa$yi)-1, max(mc0remtfstars_osa$yi)+1), ylab = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlab = "NOS quality stars", main = "Meta-regression of the comparative effect in single-center studies", cex.main = 2, cex.lab = 1.5, xaxt = "n", yaxt = "n")
lines(starsvecsingle_osa, predsstarssingle2_osa$pred, lwd = 2, col = "lightblue")
lines(starsvecsingle_osa, predsstarssingle2_osa$ci.lb, lty = "dashed")
lines(starsvecsingle_osa, predsstarssingle2_osa$ci.ub, lty = "dashed")
abline(h=0, lwd = 2, lty = "dotted")
axis(side = 2, at = seq(round(min(mc0remtfstars_osa$yi)-1,0), round(max(mc0remtfstars_osa$yi)+1,0), 1), cex.axis = 1.5)
axis(side = 1, cex.axis = 1.5)
dev.off()


################################ SA of matched & ROBINS-I: Low (ML) studies ########################################

# We create a new folder named "Plots_ML" on desktop 
folder_name <- "Plots_ML"
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

# And now we set the working directory to the new folder 
# This is where the plots are to be saved
setwd("C:/Users/sotbi/Desktop/Plots_ML")

# We will proceed with the MA on the new dataset: dset_ml
# First we initialize the dset data.frame
dset<-Q_MC 

# We define the order of reporting the different subgroups
order_post2018 <- c("Studies published after 2018", "Studies published before 2018")
order_matching <- c("Studies with patient matching", "Studies without patient matching")
order_center <- c("Multicenter studies", "Single-center studies")
order_robins <- c("Low Risk of Bias", "Moderate Risk of Bias", "Serious Risk of Bias")

# Then we apply the predefined order of reporting to the dataset
dset$post2018 <- factor(dset$post2018, levels = order_post2018)
dset$matching <- factor(dset$matching, levels = order_matching)
dset$center <- factor(dset$center, levels = order_center)
dset$robins <- factor(dset$robins, levels = order_robins)


# Then we are going to transform the dset to meet the requirements of concomitant patient matching & ROBINS-I: Low
dset_ml <- dset[dset$matching=="Studies with patient matching" & dset$robins=="Low Risk of Bias",]

# Now we define the width & height parameters for the new forest plots
# Pooled analysis forest plot parameters
W_for_pool <- 800
H_for_pool <- 400


# Pooled meta analysis
png("1. Forest plot - ML.png", width = W_for_pool, height = H_for_pool)
mc0re_ml<-metacont(nexp,mexp,sexp,nctrl,mctrl,sctrl,data = dset_ml, studlab = study, comb.fixed = FALSE, hakn = TRUE)
forestmc0re_ml<-forest(mc0re_ml, text.addline1 = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlim = c(min(round(mc0re_ml$lower,3))-0.1, max(round(mc0re_ml$upper,3))+0.1), showweights = TRUE, digits = 3, digits.se = 3, fs.xlab = 1, label.e = "RPN / RAPN", label.c = "OPN", label.right = "Favours OPN", label.left = "Favours RPN / RAPN", just = "center")
dev.off()

png("2. Funnel plot - ML.png", width = W_fr, height = H_fr)
par(mar = c(6, 6, 5, 4), mgp = c(3.75, 1, 0)) 
funnelmc0re_ml<-funnel(mc0re_ml, pch = 20,cex = 1, contour = c(0.9,0.95,0.99), col.contour = c("darkgray", "gray", "lightgray"), xlim = c(min(round(mc0re_ml$TE,3))-0.1, max(round(mc0re_ml$TE,3))+2.2), xlab = "", ylab = "", axes = FALSE)
coords <- par("usr")
rect(coords[1], coords[3], coords[2], coords[4], border = "black", lwd = 1)
axis(2, cex.axis = 1.5, pos = coords[1])
axis(1, cex.axis = 1.5, pos = coords[3])
legend(min(round(mc0re_ml$TE,3)), 0.025, c("0.1 > p > 0.05", "0.05 > p > 0.01", "p < 0.01"), fill = c("darkgray", "gray", "lightgray"), bty = "n", cex = 1.5)
title("Funnel plot with contours of statistical significance", cex.main = 2, xlab = "Mean Difference (MD)", ylab = "Standard Error (SE)", cex.lab = 1.5)
dev.off()

png("3. Radial plot with Egger's test - ML.png", width = W_fr, height = H_fr)
radialmc0re_ml<-radial(mc0re_ml)
eggersmc0re_ml<-metabias(mc0re_ml, method.bias = "linreg", plotit = TRUE, k.min = 5)
eggersmc0re_ml
regmc0re_ml<-lm(I(mc0re_ml$TE/mc0re_ml$seTE) ~ I(1/mc0re_ml$seTE))
par(mar = c(6, 6, 5, 4), mgp = c(3.75, 1, 0))
radial(mc0re_ml, cex = 1.5, cex.lab = 1.5, axes = FALSE)
coords <- par("usr")
rect(coords[1], coords[3], coords[2], coords[4], border = "black", lwd = 1)
axis(2, cex.axis = 1.5, pos = coords[1])
axis(1, cex.axis = 1.5, pos = coords[3])
abline(regmc0re_ml)
title("Radial plot with a solid regression line for Egger's test", cex.main = 2)
dev.off()

# Small study Effects
png("4. Funnel plot with small study effects - ML.png", width = W_fr, height = H_fr)
par(mar = c(6, 6, 5, 4), mgp = c(3.75, 1, 0)) 
l2_ml<-limitmeta(mc0re_ml)
funnel(l2_ml, cex = 1.5, xlim = c(min(round(mc0re_ml$TE,3))-0.1, max(round(mc0re_ml$TE,3))+1), xlab = "", ylab = "", axes = FALSE)
coords <- par("usr")
rect(coords[1], coords[3], coords[2], coords[4], border = "black", lwd = 1)
axis(2, cex.axis = 1.5, pos = coords[1])
axis(1, cex.axis = 1.5, pos = coords[3])
print(l2_ml, digits = 3)
title("Funnel plot with a curved regression line for small study effects", cex.main = 2, xlab = "Mean Difference (MD)", ylab = "Standard Error (SE)", cex.lab = 1.5)
dev.off()


# Now we will proceed with the MRA on the new dataset: dset_ml

# Meta Regression Analysis for publication year in all studies (using metafor)
png("5. MRA plot according to publication year - ML.png", width = W_pool, height = H_pool)
par(mar = c(5, 5.5, 4, 2) + 0.1)
mc0remtf_ml<-escalc(measure = "MD",n1i = nexp, m1i = mexp, sd1i = sexp, n2i = nctrl,m2i = mctrl, sd2i = sctrl, data = dset_ml, append = TRUE)
mc0remtfma_ml<-rma(yi, vi, data = mc0remtf_ml, method = "REML", mods = ~year)
yearvec_ml<-seq(min(dset_ml$year), max(dset_ml$year), 1)
preds2_ml<-predict(mc0remtfma_ml, newmods = yearvec_ml)
wi2_ml<-1/sqrt(mc0remtf_ml$vi+mc0remtfma_ml$tau2)
size2_ml<-1+2*(wi2_ml - min(wi2_ml)) / (max(wi2_ml) - min(wi2_ml))
plot(mc0remtf_ml$year, mc0remtf_ml$yi, pch = 1, cex = size2_ml, xlim = c(min(mc0remtf_ml$year)-1, max(mc0remtf_ml$year)+1), ylim = c(min(mc0remtf_ml$yi)-1, max(mc0remtf_ml$yi)+1), ylab = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlab = "Year of publication", main = "Meta-regression of the comparative effect in all studies", cex.main = 2, cex.lab = 1.5, xaxt = "n", yaxt = "n")
lines(yearvec_ml, preds2_ml$pred, lwd = 2, col = "red")
lines(yearvec_ml, preds2_ml$ci.lb, lty = "dashed")
lines(yearvec_ml, preds2_ml$ci.ub, lty = "dashed")
abline(h=0, lwd = 2, lty = "dotted")
axis(side = 2, at = seq(round(min(mc0remtf_ml$yi)-1,0), round(max(mc0remtf_ml$yi)+1,0), 1), cex.axis = 1.5)
axis(side = 1, cex.axis = 1.5)
dev.off()


png("6. MRA plot according to NOS quality stars - ML.png", width = W_pool, height = H_pool)
par(mar = c(5, 5.5, 4, 2) + 0.1)
# Meta Regression Analysis for Newcastle-Ottawa Scale (NOS) quality stars in all studies (using metafor)
mc0remtfstars_ml<-escalc(measure = "MD",n1i = nexp, m1i = mexp, sd1i = sexp, n2i = nctrl,m2i = mctrl, sd2i = sctrl, data = dset_ml, append = TRUE)
mc0remtfmastars_ml<-rma(yi, vi, data = mc0remtfstars_ml, method = "REML", mods = ~stars)
starsvec_ml<-seq(min(dset_ml$stars), max(dset_ml$stars), 1)
predsstars2_ml<-predict(mc0remtfmastars_ml, newmods = starsvec_ml)
wistars2_ml<-1/sqrt(mc0remtfstars_ml$vi+mc0remtfmastars_ml$tau2)
sizestars2_ml<-1+2*(wistars2_ml - min(wistars2_ml)) / (max(wistars2_ml) - min(wistars2_ml))
plot(mc0remtfstars_ml$stars, mc0remtfstars_ml$yi, pch = 1, cex = sizestars2_ml, xlim = c(min(mc0remtfstars_ml$stars)-1, max(mc0remtfstars_ml$stars)+1), ylim = c(min(mc0remtfstars_ml$yi)-1, max(mc0remtfstars_ml$yi)+1), ylab = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlab = "NOS quality stars", main = "Meta-regression of the comparative effect in all studies", cex.main = 2, cex.lab = 1.5, xaxt = "n", yaxt = "n")
lines(starsvec_ml, predsstars2_ml$pred, lwd = 2, col = "red")
lines(starsvec_ml, predsstars2_ml$ci.lb, lty = "dashed")
lines(starsvec_ml, predsstars2_ml$ci.ub, lty = "dashed")
abline(h=0, lwd = 2, lty = "dotted")
axis(side = 2, at = seq(round(min(mc0remtfstars_ml$yi)-1,0), round(max(mc0remtfstars_ml$yi)+1,0), 1), cex.axis = 1.5)
axis(side = 1, cex.axis = 1.5)
dev.off()


######################################### SA of Large Studies (LS) ############################################

# We create a new folder named "Plots_LS" on desktop 
folder_name <- "Plots_LS"
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

# And now we set the working directory to the new folder 
# This is where the plots are to be saved
setwd("C:/Users/sotbi/Desktop/Plots_LS")

# We will proceed with the MA on the new dataset: dset_ls
# First we initialize the dset data.frame
dset<-Q_MC 

# We define the order of reporting the different subgroups
order_post2018 <- c("Studies published after 2018", "Studies published before 2018")
order_matching <- c("Studies with patient matching", "Studies without patient matching")
order_center <- c("Multicenter studies", "Single-center studies")
order_robins <- c("Low Risk of Bias", "Moderate Risk of Bias", "Serious Risk of Bias")

# Then we apply the predefined order of reporting to the dataset
dset$post2018 <- factor(dset$post2018, levels = order_post2018)
dset$matching <- factor(dset$matching, levels = order_matching)
dset$center <- factor(dset$center, levels = order_center)
dset$robins <- factor(dset$robins, levels = order_robins)


# We are going to identify large studies based on their populations
# First we compute the populations
d_pop <- dset$nexp + dset$nctrl
# We seek those studies with population above the avarage
ind <- which(d_pop > mean(d_pop))
dset_ls <- dset[ind,]

# Now we define the width & height parameters for the new forest plots
# Pooled analysis forest plot parameters
W_for_pool <- 800
H_for_pool <- 350


# Pooled meta analysis
png("1. Forest plot (pooled) - LS.png", width = W_for_pool, height = H_for_pool)
mc0re_ls<-metacont(nexp,mexp,sexp,nctrl,mctrl,sctrl,data = dset_ls, studlab = study, comb.fixed = FALSE, hakn = TRUE)
forestmc0re_ls<-forest(mc0re_ls, text.addline1 = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlim = c(min(round(mc0re_ls$lower,3))-0.1, max(round(mc0re_ls$upper,3))+0.1), showweights = TRUE, digits = 3, digits.se = 3, fs.xlab = 1, label.e = "RPN / RAPN", label.c = "OPN", label.right = "Favours OPN", label.left = "Favours RPN / RAPN", just = "center")
dev.off()

png("2. Funnel plot - LS.png", width = W_fr, height = H_fr)
par(mar = c(6, 6, 5, 4), mgp = c(3.75, 1, 0))
funnelmc0re_ls<-funnel(mc0re_ls, pch = 20,cex = 1, contour = c(0.9,0.95,0.99), col.contour = c("darkgray", "gray", "lightgray"), xlim = c(min(round(mc0re_ls$TE,3))-0.1, max(round(mc0re_ls$TE,3))+1.2), xlab = "", ylab = "", axes = FALSE)
coords <- par("usr")
rect(coords[1], coords[3], coords[2], coords[4], border = "black", lwd = 1)
axis(2, cex.axis = 1.5, pos = coords[1])
axis(1, cex.axis = 1.5, pos = coords[3])
legend(min(round(mc0re_ls$TE,3))-0.1, 0.025, c("0.1 > p > 0.05", "0.05 > p > 0.01", "p < 0.01"), fill = c("darkgray", "gray", "lightgray"), bty = "n", cex = 1.5)
title("Funnel plot with contours of statistical significance", cex.main = 2, xlab = "Mean Difference (MD)", ylab = "Standard Error (SE)", cex.lab = 1.5)
dev.off()

png("3. Radial plot with Egger's test - LS.png", width = W_fr, height = H_fr)
radialmc0re_ls<-radial(mc0re_ls)
eggersmc0re_ls<-metabias(mc0re_ls, method.bias = "linreg", plotit = TRUE, k.min = 5)
eggersmc0re_ls
regmc0re_ls<-lm(I(mc0re_ls$TE/mc0re_ls$seTE) ~ I(1/mc0re_ls$seTE))
par(mar = c(6, 6, 5, 4), mgp = c(3.75, 1, 0)) 
radial(mc0re_ls, cex = 1.5, cex.lab = 1.5, axes = FALSE)
coords <- par("usr")
rect(coords[1], coords[3], coords[2], coords[4], border = "black", lwd = 1)
axis(2, cex.axis = 1.5, pos = coords[1])
axis(1, cex.axis = 1.5, pos = coords[3])
abline(regmc0re_ls)
title("Radial plot with a solid regression line for Egger's test", cex.main = 2)
dev.off()


# Small study Effects
png("4. Funnel plot with small study effects - LS.png", width = W_fr, height = H_fr)
par(mar = c(6, 6, 5, 4), mgp = c(3.75, 1, 0))
l2_ls<-limitmeta(mc0re_ls)
funnel(l2_ls, cex = 1.5, xlim = c(min(round(mc0re_ls$TE,3))-0.1, max(round(mc0re_ls$TE,3))+0.5), xlab = "", ylab = "", axes = FALSE)
coords <- par("usr")
rect(coords[1], coords[3], coords[2], coords[4], border = "black", lwd = 1)
axis(2, cex.axis = 1.5, pos = coords[1])
axis(1, cex.axis = 1.5, pos = coords[3])
print(l2_ls, digits = 3)
title("Funnel plot with a curved regression line for small study effects", cex.main = 2, xlab = "Mean Difference (MD)", ylab = "Standard Error (SE)", cex.lab = 1.5)
dev.off()


# Now we will proceed with the MRA on the new dataset: dset_ls

# Meta Regression Analysis for publication year in all studies (using metafor)
png("5. MRA plot according to publication year (pooled) - LS.png", width = W_pool, height = H_pool)
par(mar = c(5, 5.5, 4, 2) + 0.1)
mc0remtf_ls<-escalc(measure = "MD",n1i = nexp, m1i = mexp, sd1i = sexp, n2i = nctrl,m2i = mctrl, sd2i = sctrl, data = dset_ls, append = TRUE)
mc0remtfma_ls<-rma(yi, vi, data = mc0remtf_ls, method = "REML", mods = ~year)
yearvec_ls<-seq(min(dset_ls$year), max(dset_ls$year), 1)
preds2_ls<-predict(mc0remtfma_ls, newmods = yearvec_ls)
wi2_ls<-1/sqrt(mc0remtf_ls$vi+mc0remtfma_ls$tau2)
size2_ls<-1+2*(wi2_ls - min(wi2_ls)) / (max(wi2_ls) - min(wi2_ls))
plot(mc0remtf_ls$year, mc0remtf_ls$yi, pch = 1, cex = size2_ls, xlim = c(min(mc0remtf_ls$year)-1, max(mc0remtf_ls$year)+1), ylim = c(min(mc0remtf_ls$yi)-1, max(mc0remtf_ls$yi)+1), ylab = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlab = "Year of publication", main = "Meta-regression of the comparative effect in all studies", cex.main = 2, cex.lab = 1.5, xaxt = "n", yaxt = "n")
lines(yearvec_ls, preds2_ls$pred, lwd = 2, col = "red")
lines(yearvec_ls, preds2_ls$ci.lb, lty = "dashed")
lines(yearvec_ls, preds2_ls$ci.ub, lty = "dashed")
abline(h=0, lwd = 2, lty = "dotted")
axis(side = 2, at = seq(round(min(mc0remtf_ls$yi)-1,0), round(max(mc0remtf_ls$yi)+1,0), 1), cex.axis = 1.5)
axis(side = 1, cex.axis = 1.5)
dev.off()


png("6. MRA plot according to NOS quality stars (pooled) - LS.png", width = W_pool, height = H_pool)
par(mar = c(5, 5.5, 4, 2) + 0.1)
# Meta Regression Analysis for Newcastle-Ottawa Scale (NOS) quality stars in all studies (using metafor)
mc0remtfstars_ls<-escalc(measure = "MD",n1i = nexp, m1i = mexp, sd1i = sexp, n2i = nctrl,m2i = mctrl, sd2i = sctrl, data = dset_ls, append = TRUE)
mc0remtfmastars_ls<-rma(yi, vi, data = mc0remtfstars_ls, method = "REML", mods = ~stars)
starsvec_ls<-seq(min(dset_ls$stars), max(dset_ls$stars), 1)
predsstars2_ls<-predict(mc0remtfmastars_ls, newmods = starsvec_ls)
wistars2_ls<-1/sqrt(mc0remtfstars_ls$vi+mc0remtfmastars_ls$tau2)
sizestars2_ls<-1+2*(wistars2_ls - min(wistars2_ls)) / (max(wistars2_ls) - min(wistars2_ls))
plot(mc0remtfstars_ls$stars, mc0remtfstars_ls$yi, pch = 1, cex = sizestars2_ls, xlim = c(min(mc0remtfstars_ls$stars)-1, max(mc0remtfstars_ls$stars)+1), ylim = c(min(mc0remtfstars_ls$yi)-1, max(mc0remtfstars_ls$yi)+1), ylab = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlab = "NOS quality stars", main = "Meta-regression of the comparative effect in all studies", cex.main = 2, cex.lab = 1.5, xaxt = "n", yaxt = "n")
lines(starsvec_ls, predsstars2_ls$pred, lwd = 2, col = "red")
lines(starsvec_ls, predsstars2_ls$ci.lb, lty = "dashed")
lines(starsvec_ls, predsstars2_ls$ci.ub, lty = "dashed")
abline(h=0, lwd = 2, lty = "dotted")
axis(side = 2, at = seq(round(min(mc0remtfstars_ls$yi)-1,0), round(max(mc0remtfstars_ls$yi)+1,0), 1), cex.axis = 1.5)
axis(side = 1, cex.axis = 1.5)
dev.off()


################################ SA of matched & multicenter (MM) studies ########################################

# We create a new folder named "Plots_MM" on desktop 
folder_name <- "Plots_MM"
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

# And now we set the working directory to the new folder 
# This is where the plots are to be saved
setwd("C:/Users/sotbi/Desktop/Plots_MM")

# We will proceed with the MA on the new dataset: dset_mm
# First we initialize the dset data.frame
dset<-Q_MC 

# We define the order of reporting the different subgroups
order_post2018 <- c("Studies published after 2018", "Studies published before 2018")
order_matching <- c("Studies with patient matching", "Studies without patient matching")
order_center <- c("Multicenter studies", "Single-center studies")
order_robins <- c("Low Risk of Bias", "Moderate Risk of Bias", "Serious Risk of Bias")

# Then we apply the predefined order of reporting to the dataset
dset$post2018 <- factor(dset$post2018, levels = order_post2018)
dset$matching <- factor(dset$matching, levels = order_matching)
dset$center <- factor(dset$center, levels = order_center)
dset$robins <- factor(dset$robins, levels = order_robins)


# Then we are going to transform the dset to meet the requirements of multicenter studies with patient matching
dset_mm <- dset[dset$matching=="Studies with patient matching" & dset$center=="Multicenter studies",]

# Now we define the width & height parameters for the new forest plots
# Pooled analysis forest plot parameters
W_for_pool <- 800
H_for_pool <- 300


# Pooled meta analysis
png("1. Forest plot - MM.png", width = W_for_pool, height = H_for_pool)
mc0re_mm<-metacont(nexp,mexp,sexp,nctrl,mctrl,sctrl,data = dset_mm, studlab = study, comb.fixed = FALSE, hakn = TRUE)
forestmc0re_mm<-forest(mc0re_mm, text.addline1 = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlim = c(min(round(mc0re_mm$lower,3))-0.1, max(round(mc0re_mm$upper,3))+0.1), showweights = TRUE, digits = 3, digits.se = 3, fs.xlab = 1, label.e = "RPN / RAPN", label.c = "OPN", label.right = "Favours OPN", label.left = "Favours RPN / RAPN", just = "center")
dev.off()

png("2. Funnel plot - MM.png", width = W_fr, height = H_fr)
par(mar = c(6, 6, 5, 4), mgp = c(3.75, 1, 0)) 
funnelmc0re_mm<-funnel(mc0re_mm, pch = 20,cex = 1, contour = c(0.9,0.95,0.99), col.contour = c("darkgray", "gray", "lightgray"), xlim = c(min(round(mc0re_mm$TE,3))-0.1, max(round(mc0re_mm$TE,3))+0.8), xlab = "", ylab = "", axes = FALSE)
coords <- par("usr")
rect(coords[1], coords[3], coords[2], coords[4], border = "black", lwd = 1)
axis(2, cex.axis = 1.5, pos = coords[1])
axis(1, cex.axis = 1.5, pos = coords[3])
legend(min(round(mc0re_mm$TE,3))-0.1, 0.01, c("0.1 > p > 0.05", "0.05 > p > 0.01", "p < 0.01"), fill = c("darkgray", "gray", "lightgray"), bty = "n", cex = 1.5)
title("Funnel plot with contours of statistical significance", cex.main = 2, xlab = "Mean Difference (MD)", ylab = "Standard Error (SE)", cex.lab = 1.5)
dev.off()

png("3. Radial plot with Egger's test - MM.png", width = W_fr, height = H_fr)
radialmc0re_mm<-radial(mc0re_mm)
eggersmc0re_mm<-metabias(mc0re_mm, method.bias = "linreg", plotit = TRUE, k.min = 5)
eggersmc0re_mm
regmc0re_mm<-lm(I(mc0re_mm$TE/mc0re_mm$seTE) ~ I(1/mc0re_mm$seTE))
par(mar = c(6, 6, 5, 4), mgp = c(3.75, 1, 0))
radial(mc0re_mm, cex = 1.5, cex.lab = 1.5, axes = FALSE)
coords <- par("usr")
rect(coords[1], coords[3], coords[2], coords[4], border = "black", lwd = 1)
axis(2, cex.axis = 1.5, pos = coords[1])
axis(1, cex.axis = 1.5, pos = coords[3])
abline(regmc0re_mm)
title("Radial plot with a solid regression line for Egger's test", cex.main = 2)
dev.off()

# Small study Effects
png("4. Funnel plot with small study effects - MM.png", width = W_fr, height = H_fr)
par(mar = c(6, 6, 5, 4), mgp = c(3.75, 1, 0)) 
l2_mm<-limitmeta(mc0re_mm)
funnel(l2_mm, cex = 1.5, xlim = c(min(round(mc0re_mm$TE,3))-0.1, max(round(mc0re_mm$TE,3))+0.1), xlab = "", ylab = "", axes = FALSE)
coords <- par("usr")
rect(coords[1], coords[3], coords[2], coords[4], border = "black", lwd = 1)
axis(2, cex.axis = 1.5, pos = coords[1])
axis(1, cex.axis = 1.5, pos = coords[3])
print(l2_mm, digits = 3)
title("Funnel plot with a curved regression line for small study effects", cex.main = 2, xlab = "Mean Difference (MD)", ylab = "Standard Error (SE)", cex.lab = 1.5)
dev.off()


# Now we will proceed with the MRA on the new dataset: dset_mm

# Meta Regression Analysis for publication year in all studies (using metafor)
png("5. MRA plot according to publication year - MM.png", width = W_pool, height = H_pool)
par(mar = c(5, 5.5, 4, 2) + 0.1)
mc0remtf_mm<-escalc(measure = "MD",n1i = nexp, m1i = mexp, sd1i = sexp, n2i = nctrl,m2i = mctrl, sd2i = sctrl, data = dset_mm, append = TRUE)
mc0remtfma_mm<-rma(yi, vi, data = mc0remtf_mm, method = "REML", mods = ~year)
yearvec_mm<-seq(min(dset_mm$year), max(dset_mm$year), 1)
preds2_mm<-predict(mc0remtfma_mm, newmods = yearvec_mm)
wi2_mm<-1/sqrt(mc0remtf_mm$vi+mc0remtfma_mm$tau2)
size2_mm<-1+2*(wi2_mm - min(wi2_mm)) / (max(wi2_mm) - min(wi2_mm))
plot(mc0remtf_mm$year, mc0remtf_mm$yi, pch = 1, cex = size2_mm, xlim = c(min(mc0remtf_mm$year)-1, max(mc0remtf_mm$year)+1), ylim = c(min(mc0remtf_mm$yi)-1, max(mc0remtf_mm$yi)+1), ylab = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlab = "Year of publication", main = "Meta-regression of the comparative effect in all studies", cex.main = 2, cex.lab = 1.5, xaxt = "n", yaxt = "n")
lines(yearvec_mm, preds2_mm$pred, lwd = 2, col = "red")
lines(yearvec_mm, preds2_mm$ci.lb, lty = "dashed")
lines(yearvec_mm, preds2_mm$ci.ub, lty = "dashed")
abline(h=0, lwd = 2, lty = "dotted")
axis(side = 2, at = seq(round(min(mc0remtf_mm$yi)-1,0), round(max(mc0remtf_mm$yi)+1,0), 1), cex.axis = 1.5)
axis(side = 1, cex.axis = 1.5)
dev.off()


png("6. MRA plot according to NOS quality stars - MM.png", width = W_pool, height = H_pool)
par(mar = c(5, 5.5, 4, 2) + 0.1)
# Meta Regression Analysis for Newcastle-Ottawa Scale (NOS) quality stars in all studies (using metafor)
mc0remtfstars_mm<-escalc(measure = "MD",n1i = nexp, m1i = mexp, sd1i = sexp, n2i = nctrl,m2i = mctrl, sd2i = sctrl, data = dset_mm, append = TRUE)
mc0remtfmastars_mm<-rma(yi, vi, data = mc0remtfstars_mm, method = "REML", mods = ~stars)
starsvec_mm<-seq(min(dset_mm$stars), max(dset_mm$stars), 1)
predsstars2_mm<-predict(mc0remtfmastars_mm, newmods = starsvec_mm)
wistars2_mm<-1/sqrt(mc0remtfstars_mm$vi+mc0remtfmastars_mm$tau2)
sizestars2_mm<-1+2*(wistars2_mm - min(wistars2_mm)) / (max(wistars2_mm) - min(wistars2_mm))
plot(mc0remtfstars_mm$stars, mc0remtfstars_mm$yi, pch = 1, cex = sizestars2_mm, xlim = c(min(mc0remtfstars_mm$stars)-1, max(mc0remtfstars_mm$stars)+1), ylim = c(min(mc0remtfstars_mm$yi)-1, max(mc0remtfstars_mm$yi)+1), ylab = "Mean difference in per minute blood loss (RPN / RAPN vs. OPN) in ml/min", xlab = "NOS quality stars", main = "Meta-regression of the comparative effect in all studies", cex.main = 2, cex.lab = 1.5, xaxt = "n", yaxt = "n")
lines(starsvec_mm, predsstars2_mm$pred, lwd = 2, col = "red")
lines(starsvec_mm, predsstars2_mm$ci.lb, lty = "dashed")
lines(starsvec_mm, predsstars2_mm$ci.ub, lty = "dashed")
abline(h=0, lwd = 2, lty = "dotted")
axis(side = 2, at = seq(round(min(mc0remtfstars_mm$yi)-1,0), round(max(mc0remtfstars_mm$yi)+1,0), 1), cex.axis = 1.5)
axis(side = 1, cex.axis = 1.5)
dev.off()


