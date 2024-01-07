############################################### MAIN FUNCTION ##########################################################

# First we create a function to estimate Q from EBL & OT according to: 
# van Kempen et al. 2000: "Mean and Variance of Ratio Estimators Used in Fluorescence Ratio Imaging"

van_Kempen <- function(mEBL, mOT, sEBL, sOT, rho) {
  n <- length(mEBL)  # Length of the vectors
  
  mQ <- numeric(n)  # Initialize mQ vector
  sQ <- numeric(n)  # Initialize sQ vector
  
  for (i in 1:n) {
    mQ[i] <- (mEBL[i] / mOT[i]) + ((mEBL[i] * sOT[i]^2) / mOT[i]^3) - ((rho[i] * sEBL[i] * sOT[i]) / mOT[i]^2)
    
    sQ[i] <- sqrt((sEBL[i]^2 / mOT[i]^2) + ((mEBL[i]^2 * sOT[i]^2) / mOT[i]^4) - ((2 * rho[i] * mEBL[i] * sEBL[i] * sOT[i]) / mOT[i]^3))
  }
  
  return(list(mQ = mQ, sQ = sQ))
}

############################################### DATASETS ############################################################

# Then we define the datasets of interest
EBL<-read.csv("C:/Users/sotbi/Desktop/Data/EBL.csv", sep=";", header=TRUE)
OT<-read.csv("C:/Users/sotbi/Desktop/Data/OT.csv", sep=";", header=TRUE)
RS<-read.csv("C:/Users/sotbi/Desktop/Data/RS.csv", sep=";", header=TRUE)

####################################### CORRELATION COEFFICIENT #####################################################

# Theoretical values of the correlation coefficient (r)
rhos <- c(-0.99, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0 ,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99)

####################################### Q FOR MONTE CARLO (Q_MC) ###################################################

# We obtain the main structure of Q_MC
Q_MC <- EBL
Q_MC[,c("mexp","sexp","mctrl","sctrl")] <- NA

# We create a new column "r_exp" for the correlation coefficients that we obtained for the experimental arm (RPN / RAPN) with the Monte Carlo simulation
r_exp <- RS$Estimated.r.in.exp
# Now we determine the position to insert the new column
insert_after_column <- "sexp"
# Then we get the index of the column after which to insert the new column
insert_index <- match(insert_after_column, names(Q_MC))
# Finally we insert the new column after the specified column
Q_MC <- cbind(Q_MC[, 1:insert_index], r_exp, Q_MC[, (insert_index + 1):ncol(Q_MC)])

# Now we create a new column "r_ctrl" for the correlation coefficients that we obtained for the control arm (OPN) with the Monte Carlo simulation
Q_MC$r_ctrl <- RS$Estimated.r.in.ctrl

# We attach the EBL & OT data to Q_MC
Q_MC$mEBL_exp <- EBL$mexp
Q_MC$mOT_exp <- OT$mexp
Q_MC$sEBL_exp <- EBL$sexp
Q_MC$sOT_exp <- OT$sexp
Q_MC$mEBL_ctrl <- EBL$mctrl
Q_MC$mOT_ctrl <- OT$mctrl
Q_MC$sEBL_ctrl <- EBL$sctrl
Q_MC$sOT_ctrl <- OT$sctrl

# For the Experimental arm (RPN / RAPN)
Q_MC_exp <- van_Kempen(Q_MC$mEBL_exp, Q_MC$mOT_exp, Q_MC$sEBL_exp, Q_MC$sOT_exp, Q_MC$r_exp)
Q_MC$mexp <- Q_MC_exp$mQ 
Q_MC$sexp <- Q_MC_exp$sQ 

# For the Experimental arm (RPN / RAPN)
Q_MC_ctrl <- van_Kempen(Q_MC$mEBL_ctrl, Q_MC$mOT_ctrl, Q_MC$sEBL_ctrl, Q_MC$sOT_ctrl, Q_MC$r_ctrl)
Q_MC$mctrl <- Q_MC_ctrl$mQ 
Q_MC$sctrl <- Q_MC_ctrl$sQ 

# Now we will remove the redundant columns and refine Q_MC with precision of 5 decimal places
columns_to_remove <- c("r_exp", "r_ctrl", "mEBL_exp", "mOT_exp", "sEBL_exp", "sOT_exp", "mEBL_ctrl", "mOT_ctrl", "sEBL_ctrl", "sOT_ctrl")
Q_MC <- Q_MC[, !(names(Q_MC) %in% columns_to_remove)]
columns_to_round <- c("mexp", "sexp", "mctrl", "sctrl")
Q_MC[, columns_to_round] <- round(Q_MC[, columns_to_round], 5)


################################################ EXCEL FOLDER ######################################################

# We create a new folder named "Excel" on desktop 
folder_name <- "Excel"
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


################################################ SAVE Q_MC ####################################################

# Export the dataset as an MS Excel file
install.packages("openxlsx")
library(openxlsx)

# Specify the file path for the Excel file
file_path <- "C:/Users/sotbi/Desktop/Excel/Q_MC.xlsx"

# Create a new workbook
wb <- createWorkbook()

# Add the dataset to the workbook
addWorksheet(wb, "Sheet1")
writeData(wb, "Sheet1", Q_MC)

# Specify the number of decimal places for the Value columns
style <- createStyle(numFmt = "#,##0.000000")
addStyle(wb, "Sheet1", style, rows = NULL, cols = "B")

# Save the workbook as an Excel file
saveWorkbook(wb, file_path)

# Print a confirmation message
cat("Data exported successfully to", file_path)



############################### Q FOR SENSITIVITY ANALYSIS (Q_SA) ####################################################

# We obtain the main structure of Q_SA
Q_SA0 <- Q_MC
Q_SA0$mEBL_exp <- EBL$mexp
Q_SA0$mOT_exp <- OT$mexp
Q_SA0$sEBL_exp <- EBL$sexp
Q_SA0$sOT_exp <- OT$sexp
Q_SA0$mEBL_ctrl <- EBL$mctrl
Q_SA0$mOT_ctrl <- OT$mctrl
Q_SA0$sEBL_ctrl <- EBL$sctrl
Q_SA0$sOT_ctrl <- OT$sctrl

# Copy lines and append 21 times (length(rhos))
Q_SA <- Q_SA0
k <-length(rhos)-1
for (i in 1:k) {
  Q_SA <- rbind(Q_SA, Q_SA0)
}

# Reset row names of dataset "Q_SA"
rownames(Q_SA) <- NULL

# Assign the values of r (correlation coefficient)
Q_SA$r <- rep(rhos, each = length(EBL$study))

# For the Experimental arm (RPN / RAPN)
Q_SA_exp <- van_Kempen(Q_SA$mEBL_exp, Q_SA$mOT_exp, Q_SA$sEBL_exp, Q_SA$sOT_exp, Q_SA$r)
Q_SA$mexp <- Q_SA_exp$mQ 
Q_SA$sexp <- Q_SA_exp$sQ 

# For the Control arm (OPN)
Q_SA_ctrl <- van_Kempen(Q_SA$mEBL_ctrl, Q_SA$mOT_ctrl, Q_SA$sEBL_ctrl, Q_SA$sOT_ctrl, Q_SA$r)
Q_SA$mctrl <- Q_SA_ctrl$mQ 
Q_SA$sctrl <- Q_SA_ctrl$sQ 

# Now we will remove the redundant columns and refine Q_SA with precision of 3 decimal places
columns_to_remove <- c("mEBL_exp", "mOT_exp", "sEBL_exp", "sOT_exp", "mEBL_ctrl", "mOT_ctrl", "sEBL_ctrl", "sOT_ctrl")
Q_SA <- Q_SA[, !(names(Q_SA) %in% columns_to_remove)]
columns_to_round <- c("mexp", "sexp", "mctrl", "sctrl")
Q_SA[, columns_to_round] <- round(Q_SA[, columns_to_round], 3)



############################################### SAVE Q_SA #########################################################

# Specify the file path for the Excel file
file_path <- "C:/Users/sotbi/Desktop/Excel/Q_SA.xlsx"

# Create a new workbook
wb <- createWorkbook()

# Add the dataset to the workbook
addWorksheet(wb, "Sheet1")
writeData(wb, "Sheet1", Q_SA)

# Specify the number of decimal places for the Value columns
style <- createStyle(numFmt = "#,##0.000000")
addStyle(wb, "Sheet1", style, rows = NULL, cols = "B")

# Save the workbook as an Excel file
saveWorkbook(wb, file_path)

# Print a confirmation message
cat("Data exported successfully to", file_path)

