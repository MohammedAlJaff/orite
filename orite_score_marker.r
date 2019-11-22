### To run the follwoing code, the package filesstrings is needed. 
### To run the script: set the working directory to the one containing the files to
### be scored. Mark origin and terminus with one click each. The first click represents the
### origin and the second click represents the terminus. 
### After this you will be able to provide a score of how certain you were of the prediction. 
### All the information will be stored in the file "Orite_Score.txt" which is comma separated. 
### The original files that have been processed gets move from the current directory to 
### a "Processed" directory.
### If you accidentaly mark a position you will have to manually remove the position
### from "Orite_Score.txt". 


rm(list = ls(all.names = TRUE))
library(filesstrings)



# Creates the file "Orite_Score" to save the results in.

if(!file.exists("Orite_Score.txt")){
  y <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(y) <- c("ID","Ori","Ter","Score")
  write.table(y, file = "Orite_Score.txt", sep = ",",row.names = FALSE, col.names = TRUE)
}


# Extracts all the file names in current directory
current_dir <- getwd()
file_names <- list.files(path= current_dir, pattern = "[:digit:]")


for(name in file_names){
# Decide position for Origin,Terminus and assign a Score where 1 is cluess and 5 is very certain.
  df <- read.table(name, h=T)
  plot(df$pos, df$cGCsk, t="l")
  x <- locator(2)
  Ori <- floor(x$x[1])
  Ter <- floor(x$x[2])
  Score <- readline("How certain are you of this prediction? (1-5, where 1 is clueless and 5 is very certain): ")
  Score <- strtoi(Score, base = 0L)
 
# Removes the file extension.
   name = unlist(strsplit(name, split='.', fixed=TRUE))[1]
   
   test_df <- data.frame(name,Ori,Ter,Score )
# Add the data to the "Orite_Score" text file. 
     write.table(test_df,file = "Orite_Score.txt", append = TRUE, row.names = FALSE, col.names = FALSE, sep = ",")

# Removes the files that have been scored and move them from current directory to processed directory     
  processed_dir = paste(current_dir, '/processed', sep = '')
  processed_name = paste(current_dir, '/', name, '.txt', sep = '')
  file.move(processed_name, processed_dir)
  setwd(current_dir)
}