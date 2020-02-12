# Code to convert MSDIAL's default output of 
# .txt files to .csv files with nice names

filenames <- list.files("data_raw", pattern = "\\.txt", full.names = TRUE)

filedata <- lapply(X = filenames, FUN = read.csv, 
                   sep = "\t")
new_filenames <- gsub(pattern = "txt", replacement = "csv", filenames)
new_filenames <- gsub(pattern = "_0_[[:digit:]]*", replacement = "_HILICPos", new_filenames)
new_filenames <- gsub(pattern = "_1_[[:digit:]]*", replacement = "_HILICNeg", new_filenames)
names(filedata) <- new_filenames

sapply(seq_along(filedata), function(x){
  write.table(filedata[[x]], file = names(filedata)[x], 
            sep=",", col.names = FALSE, row.names = FALSE)
})

file.remove(filenames)
