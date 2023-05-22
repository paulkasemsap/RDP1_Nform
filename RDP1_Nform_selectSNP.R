+++++++++++++++++++++++++++++++++
#### Mission Accomplished: v 1.3 (May 2023).
#### Goals: 
  # 1) Select significant SNPs from GWAS results of each trait for further analyses and report in the publication.
  # 2) Filter SNPs efficiently across populations, treatments and models.
  # 3) Visualize significant SNPs.
#### Last update: Based on the source file "gwas-selectSNP.r" last committed on 5/22/2023.
+++++++++++++++++++++++++++++++++
#### Prerequisite ####
{#*---- Install & Load Required Packages ----
workstation <- "C:/Users/paulk/Box/RiceGWAS2023" #Change prefix of work station here first!
  
packagelist <- c(
    "rMVP", #manhattan/qq plots
    "tidyverse", #process data
    "data.table", #process data
    "viridis", #figure color
    "patchwork",#to arrange figures
    "UpSetR", #upsetR for visualizing sets
    "ComplexUpset"
) #Add packages into this list for ease of installation and loading
#install.packages(packagelist) #Only use once
lapply(packagelist, library, character.only = TRUE) #Load the required packages
}
getDTthreads() #check number of threads; should be 4 for higher efficiency
writeLines(capture.output(sessionInfo()), "sessionInfo.txt") #save session info and check for compatibility later, if needed
+++++++++++++++++++++++++++++++++
#### Working Code - Compile and Filter P-values from GAPIT GWAS results for EACH population and EACH model ####
#GAPIT periodically change file name formatting! Must check file name format before running codes!
#Depending on how GAPIT was run, need to compile/select the result files before starting this step.
#In this workflow, results for each population were saved in individual folders, but the file name contains strings that indicate "model", "nitrogen treatment", and "organ". Hence, the codes below select the result files based on the file names.

{#OK*---- Intermediate result "Sig[Exp or Marker]WiseSNP.list" ====
    {#* RDP1 - Specify work directory and population ####
      setwd(paste0(workstation,"/rdp1"))
      pop = "rdp1-pc3"
      {#|--------Step 1 Merge p-value across traits--------| ####
        start.time <- Sys.time()
        writeLines(capture.output(sessionInfo()), "sessionInfo_mergesnp.txt") #save session info and check for compatibility later, if needed
        writeLines(capture.output(getwd()), "wd_mergesnp.txt")#wd
        {#*1. BLINK root - Merge p-value across traits for Manhattan Plot ----
          {#*---- a. Specify organ to automate file selection ----
            method = "Blink"
            organ = "root"
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
        }
        {#*2. BLINK shoot - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "Blink"
            organ = "shoot"
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
        }
        {#*3. BLINK whole - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "Blink"
            organ = "whole"
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*4.  Farm root - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "FarmCPU"
            
            organ = "root"
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*5. Farm shoot - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "FarmCPU"
            
            
            organ = "shoot"
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*6.  Farm whole - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "FarmCPU"
            
            organ = "whole"
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*7. MLM root - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "MLM"
            
            organ = "root"
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*8.  MLM shoot - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "MLM"
            
            
            organ = "shoot"
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*9. MLM whole - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            
            method = "MLM"
            
            organ = "whole"
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*10.  GLM root - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            
            method = "GLM"
            organ = "root"
            
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*11.  GLM shoot - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            
            method = "GLM"
            
            
            organ = "shoot"
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*12.  GLM whole - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            
            method = "GLM"
            
            organ = "whole"
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
        }
        
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        time.taken
        write(time.taken,"timetaken_mergesnp.txt")
      }
      {#|--------Step 2 Filter significant SNPs--------| ####
        start.time <- Sys.time()
        writeLines(capture.output(sessionInfo()), "sessionInfo_filtersnp.txt") #save session info and check for compatibility later, if needed
        writeLines(capture.output(getwd()), "wd_filtersnp.txt")#wd
        {#*---- a. Prepare data files ----
          # Read file that contains values for manhattan plots
          EditTextHere <- "*.ManhattanPlot" #Edit file name here
          dir(pattern=EditTextHere) #Confirm selection
          file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
          
          # Read csv files for manhattan plots
          allpvalue <- lapply(file.list,function(x) {
            file <- fread(x) #read the data into R
            file <- rename_with(file, ~ gsub("_whole", "_", .x, fixed = T))
            file <- rename_with(file, ~ gsub("_root", "_", .x, fixed = T))
            file <- rename_with(file, ~ gsub("_shoot", "_", .x, fixed = T))
            x1 <- strsplit(x,"\\.")
            organ = x1[[1]][3]
            GWASmodel = x1[[1]][2]
            population = x1[[1]][1]
            file <- cbind (population, GWASmodel, organ, file) # add file name as a new column
            return (file)
          })
          
        }
        {#*---- a extra. For GAPIT results with different file name format; Rename column by removing prefix and sffix ----
          allpvalue<-allpvalue %>%
            map(~ .x %>%
                  rename_with(~ gsub("GAPIT.", "", .x, fixed = TRUE))%>%
                  rename_with(~ gsub(".GWAS.Results", "", .x, fixed = TRUE))%>%
                  rename_with(~ gsub("Blink.", "", .x, fixed = TRUE))%>%
                  rename_with(~ gsub("FarmCPU.", "", .x, fixed = TRUE))%>%
                  rename_with(~ gsub("MLM.", "", .x, fixed = TRUE))%>%
                  rename_with(~ gsub("GLM.", "", .x, fixed = TRUE))
            )
        }
        {#*---- b. Select significant SNPs across N treatments for each organ ----
          # Filter only SNPs with p-value below threshold (Bonferoni Correction = 0.1 for 700,000k SNP)
          SignificantMarkerWiseSNP<- allpvalue %>%
            map(~ .x %>%
                  filter_at(vars(matches('_')), any_vars(.< (0.001))) #p-value < 0.001
            )
          SignificantMarkerWise10e4SNP<- allpvalue %>%
            map(~ .x %>%
                  filter_at(vars(matches('_')), any_vars(.< (0.0001))) #p-value < 0.001
            )
          SignificantExpWiseSNP<- allpvalue %>%
            map(~ .x %>%
                  filter_at(vars(matches('_')), any_vars(.< (0.1/700000))) # Bonferoni correction alpha = 0.1/number of SNPs used.
            )
          
          # Merge all significant SNPs into one data table
          AllSigMarkerWiseSNP = rbindlist(SignificantMarkerWiseSNP,fill=F, use.names=T) #match column names by name, not fill with NAs if missing
          AllSigMarkerWise10e4SNP = rbindlist(SignificantMarkerWise10e4SNP,fill=F, use.names=T ) #match column names by name, not fill with NAs if missing
          AllSigExpWiseSNP = rbindlist(SignificantExpWiseSNP,fill=F, use.names=T ) #match column names by name, not fill with NAs if missing
          #https://www.rdocumentation.org/packages/data.table/versions/1.14.8/topics/rbindlist
          
          write.csv(AllSigMarkerWiseSNP, file=paste(pop,"AllModels","SigMarkerWiseSNP","list.csv",sep="."), row.names=FALSE)
          write.csv(AllSigMarkerWise10e4SNP, file=paste(pop,"AllModels","SigMarkerWise10e4SNP","list.csv",sep="."), row.names=FALSE)
          write.csv(AllSigExpWiseSNP, file=paste(pop,"AllModels","SigExpWiseSNP","list.csv",sep="."), row.names=FALSE)
          
          # Count number of significant SNPs
          AllSigMarkerWiseSNP%>%
            #mutate(newgroup=paste(organ,GWASmodel,sep= "."))%>%
            group_by(population,organ, GWASmodel)%>%
            tally()%>%
            write.csv(., file=paste(pop,"SigMarkerWiseSNP","count.csv",sep="."),quote=F,row.names=F)
          
          AllSigMarkerWise10e4SNP%>%
            #mutate(newgroup=paste(organ,GWASmodel,sep= "."))%>%
            group_by(population,organ, GWASmodel)%>%
            tally()%>%
            write.csv(., file=paste(pop,"SigMarkerWise10e4SNP","count.csv",sep="."),quote=F,row.names=F)
          
          AllSigExpWiseSNP%>%
            #mutate(newgroup=paste(organ,GWASmodel,sep= "."))%>%
            group_by(population,organ, GWASmodel)%>%
            tally()%>%
            write.csv(., file=paste(pop,"SigExpWiseSNP","count.csv",sep="."),quote=F,row.names=F)
          
          # Save a separate csv file for each GWAS model
          for(i in unique(AllSigMarkerWiseSNP$GWASmodel)){write.csv(AllSigMarkerWiseSNP[AllSigMarkerWiseSNP$GWASmodel==i,],paste(pop,i,"model","SigMarkerWiseSNP","list.csv",sep="."),quote=F,row.names=F)} #write 4 new files based on GWAS models (all organs)
          for(i in unique(AllSigMarkerWise10e4SNP$GWASmodel)){write.csv(AllSigMarkerWise10e4SNP[AllSigMarkerWise10e4SNP$GWASmodel==i,],paste(pop,i,"model","SigMarkerWise10e4SNP","list.csv",sep="."),quote=F,row.names=F)} #write 4 new files based on GWAS models (all organs)
          for(i in unique(AllSigExpWiseSNP$GWASmodel)){write.csv(AllSigExpWiseSNP[AllSigExpWiseSNP$GWASmodel==i,],paste(pop,i,"model","SigExpWiseSNP","list.csv",sep="."),quote=F,row.names=F)} #write 4 new files based on GWAS models (all organs)
          for(i in unique(AllSigMarkerWiseSNP$organ)){write.csv(AllSigMarkerWiseSNP[AllSigMarkerWiseSNP$organ==i,],paste(pop,i,"organ","SigMarkerWiseSNP","list.csv",sep="."),quote=F,row.names=F)} #write 3 new files based on organs (all models)
          
        }
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        time.taken
        write(time.taken,"timetaken_filterSNP.txt")
        # The SNP will be browsed against the genome browser
        #--http://rs-bt-mccouch4.biotech.cornell.edu:443/cgi-bin/hgGateway?clade=poaceae&org=O.+sativa&db=orySat2
        #--https://snp-seek.irri.org/
      }
      {#|--------Step 3 Plotting (multiple) Manhattan Plots *with only significant SNPs --------| ####
        #Currently I set it to plot manhattan plots using only sig SNPS with p-values <0.001 using the files I saved from Step 2 for EACH GWAS model. If I want to use all SNPs, that can be done, too. Just need to use the full GWAS result files, or the "ManhattanPlot" files I saved from Step 1.
        start.time <- Sys.time()
        writeLines(capture.output(sessionInfo()), "sessionInfo_manhattan.txt") #save session info and check for compatibility later, if needed
        writeLines(capture.output(getwd()), "wd_manhattan.txt")#wd
        {#*---- a.  Load data ----
          # Read file that contains values for manhattan plots; Use ONLY significant SNPs
          EditTextHere <- paste("SigMarkerWiseSNP.list")
          dir(pattern=EditTextHere) #Confirm selection
          file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
        }
        {#*---- b. Making multiple manhattan plots ----
          lapply(file.list,function(x) {
            file <- fread(x) #read the data into R
            x1 <- strsplit(x,"\\.")
            id = x1[[1]][2]
            confirm =  x1[[1]][3]
            population = x1[[1]][1]
            file <-MVP.Report(file[,4:10],
                              plot.type=c("m"),
                              multracks=T, #If this one is true, will get 2 plots:multitracks and multitraits, instead of n plot (for n traits)
                              threshold=c(0.1/700000,1e-3),
                              threshold.lty=c(2,1),
                              threshold.lwd=c(2,2),
                              threshold.col=c("black","black"),
                              amplify=F,
                              ylim=c(3,round(-log10(min(file[,7:10], na.rm=F)))+3),
                              band=3,
                              file.type="pdf",
                              memo=paste0(population,".",id,".",confirm),
                              dpi=5,
                              verbose = T)
            return (file)
          })
        }
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        time.taken
        write(time.taken,"timetaken_manhattan.txt")
      }
      {#|-------- Step 4 Plotting (multiple) QQ plots with all GWAS p-values--------| ####
        #This step can take a long time; only run when needed
        start.time <- Sys.time()
        writeLines(capture.output(sessionInfo()), "sessionInfo_QQ.txt") #save session info and check for compatibility later, if needed
        writeLines(capture.output(getwd()), "wd_QQ.txt")#wd
        {#*---- a.  Load data ----
          # Read file that contains values for manhattan plots
          EditTextHere <- paste(".ManhattanPlot")
          dir(pattern=EditTextHere) #Confirm selection
          file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
        }
        {#*---- b. Making multiple QQ plots from full GWAS results ----
          lapply(file.list,function(x) {
            file <- fread(x) #read the data into R
            file <- rename_with(file, ~ gsub("_whole_", "_", .x, fixed = T))
            file <- rename_with(file, ~ gsub("_root_", "_", .x, fixed = T))
            file <- rename_with(file, ~ gsub("_shoot_", "_", .x, fixed = T))
            file <- rename_with(file,~ gsub("GAPIT.", "", .x, fixed = TRUE))
            file <- rename_with(file,~ gsub(".GWAS.Results", "", .x, fixed = TRUE))
            file <- rename_with(file,~ gsub("Blink.", "", .x, fixed = TRUE))
            file <- rename_with(file,~ gsub("FarmCPU.", "", .x, fixed = TRUE))
            file <- rename_with(file,~ gsub("MLM.", "", .x, fixed = TRUE))
            file <- rename_with(file,~ gsub("GLM.", "", .x, fixed = TRUE))
            
            x1 <- strsplit(x,"\\.")
            GWASmodel = x1[[1]][2]
            population = x1[[1]][1]
            organ = x1[[1]][3]
            file <-MVP.Report(file[,1:7],
                              plot.type=c("q"),
                              multracks=T, #If this one is true, will get 2 plots: overlayed QQ plots and multiple QQ plots on the same plane
                              file.type="pdf",
                              memo=paste0(population,".",GWASmodel,".",organ),
                              dpi=5,
                              verbose = T)
            return (file)
          })
        }
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        time.taken
        write(time.taken,"timetaken_QQ.txt")
      }
    }
    {#* Aus - Specify work directory and population ####
      setwd(paste0(workstation,"/aus"))
      pop = "aus"
      {#|--------Step 1 Merge p-value across traits--------| ####
        start.time <- Sys.time()
        writeLines(capture.output(sessionInfo()), "sessionInfo_mergesnp.txt") #save session info and check for compatibility later, if needed
        writeLines(capture.output(getwd()), "wd_mergesnp.txt")#wd
        {#*1. BLINK root - Merge p-value across traits for Manhattan Plot ----
          {#*---- a. Specify organ to automate file selection ----
            method = "Blink"
            organ = "root"
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
        }
        {#*2. BLINK shoot - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "Blink"
            organ = "shoot"
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
        }
        {#*3. BLINK whole - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "Blink"
            organ = "whole"
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*4.  Farm root - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "FarmCPU"
            
            organ = "root"
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*5. Farm shoot - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "FarmCPU"
            
            
            organ = "shoot"
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*6.  Farm whole - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "FarmCPU"
            
            organ = "whole"
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*7. MLM root - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "MLM"
            
            organ = "root"
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*8.  MLM shoot - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "MLM"
            
            
            organ = "shoot"
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*9. MLM whole - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            
            method = "MLM"
            
            organ = "whole"
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*10.  GLM root - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            
            method = "GLM"
            organ = "root"
            
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*11.  GLM shoot - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            
            method = "GLM"
            
            
            organ = "shoot"
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*12.  GLM whole - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            
            method = "GLM"
            
            organ = "whole"
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
        }
        
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        time.taken
        write(time.taken,"timetaken_mergesnp.txt")
      }
      {#|--------Step 2 Filter significant SNPs--------| ####
        start.time <- Sys.time()
        writeLines(capture.output(sessionInfo()), "sessionInfo_filtersnp.txt") #save session info and check for compatibility later, if needed
        writeLines(capture.output(getwd()), "wd_filtersnp.txt")#wd
        {#*---- a. Prepare data files ----
          # Read file that contains values for manhattan plots
          EditTextHere <- "*.ManhattanPlot" #Edit file name here
          dir(pattern=EditTextHere) #Confirm selection
          file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
          
          # Read csv files for manhattan plots
          allpvalue <- lapply(file.list,function(x) {
            file <- fread(x) #read the data into R
            file <- rename_with(file, ~ gsub("_whole", "_", .x, fixed = T))
            file <- rename_with(file, ~ gsub("_root", "_", .x, fixed = T))
            file <- rename_with(file, ~ gsub("_shoot", "_", .x, fixed = T))
            x1 <- strsplit(x,"\\.")
            organ = x1[[1]][3]
            GWASmodel = x1[[1]][2]
            population = x1[[1]][1]
            file <- cbind (population, GWASmodel, organ, file) # add file name as a new column
            return (file)
          })
          
        }
        {#*---- a extra. For GAPIT results with different file name format; Rename column by removing prefix and sffix ----
          allpvalue<-allpvalue %>%
            map(~ .x %>%
                  rename_with(~ gsub("GAPIT.", "", .x, fixed = TRUE))%>%
                  rename_with(~ gsub(".GWAS.Results", "", .x, fixed = TRUE))%>%
                  rename_with(~ gsub("Blink.", "", .x, fixed = TRUE))%>%
                  rename_with(~ gsub("FarmCPU.", "", .x, fixed = TRUE))%>%
                  rename_with(~ gsub("MLM.", "", .x, fixed = TRUE))%>%
                  rename_with(~ gsub("GLM.", "", .x, fixed = TRUE))
            )
        }
        {#*---- b. Select significant SNPs across N treatments for each organ ----
          # Filter only SNPs with p-value below threshold (Bonferoni Correction = 0.1 for 700,000k SNP)
          SignificantMarkerWiseSNP<- allpvalue %>%
            map(~ .x %>%
                  filter_at(vars(matches('_')), any_vars(.< (0.001))) #p-value < 0.001
            )
          SignificantMarkerWise10e4SNP<- allpvalue %>%
            map(~ .x %>%
                  filter_at(vars(matches('_')), any_vars(.< (0.0001))) #p-value < 0.001
            )
          SignificantExpWiseSNP<- allpvalue %>%
            map(~ .x %>%
                  filter_at(vars(matches('_')), any_vars(.< (0.1/700000))) # Bonferoni correction alpha = 0.1/number of SNPs used.
            )
          
          # Merge all significant SNPs into one data table
          AllSigMarkerWiseSNP = rbindlist(SignificantMarkerWiseSNP,fill=F, use.names=T) #match column names by name, not fill with NAs if missing
          AllSigMarkerWise10e4SNP = rbindlist(SignificantMarkerWise10e4SNP,fill=F, use.names=T ) #match column names by name, not fill with NAs if missing
          AllSigExpWiseSNP = rbindlist(SignificantExpWiseSNP,fill=F, use.names=T ) #match column names by name, not fill with NAs if missing
          #https://www.rdocumentation.org/packages/data.table/versions/1.14.8/topics/rbindlist
          
          write.csv(AllSigMarkerWiseSNP, file=paste(pop,"AllModels","SigMarkerWiseSNP","list.csv",sep="."), row.names=FALSE)
          write.csv(AllSigMarkerWise10e4SNP, file=paste(pop,"AllModels","SigMarkerWise10e4SNP","list.csv",sep="."), row.names=FALSE)
          write.csv(AllSigExpWiseSNP, file=paste(pop,"AllModels","SigExpWiseSNP","list.csv",sep="."), row.names=FALSE)
          
          # Count number of significant SNPs
          AllSigMarkerWiseSNP%>%
            #mutate(newgroup=paste(organ,GWASmodel,sep= "."))%>%
            group_by(population,organ, GWASmodel)%>%
            tally()%>%
            write.csv(., file=paste(pop,"SigMarkerWiseSNP","count.csv",sep="."),quote=F,row.names=F)
          
          AllSigMarkerWise10e4SNP%>%
            #mutate(newgroup=paste(organ,GWASmodel,sep= "."))%>%
            group_by(population,organ, GWASmodel)%>%
            tally()%>%
            write.csv(., file=paste(pop,"SigMarkerWise10e4SNP","count.csv",sep="."),quote=F,row.names=F)
          
          AllSigExpWiseSNP%>%
            #mutate(newgroup=paste(organ,GWASmodel,sep= "."))%>%
            group_by(population,organ, GWASmodel)%>%
            tally()%>%
            write.csv(., file=paste(pop,"SigExpWiseSNP","count.csv",sep="."),quote=F,row.names=F)
          
          # Save a separate csv file for each GWAS model
          for(i in unique(AllSigMarkerWiseSNP$GWASmodel)){write.csv(AllSigMarkerWiseSNP[AllSigMarkerWiseSNP$GWASmodel==i,],paste(pop,i,"model","SigMarkerWiseSNP","list.csv",sep="."),quote=F,row.names=F)} #write 4 new files based on GWAS models (all organs)
          for(i in unique(AllSigMarkerWise10e4SNP$GWASmodel)){write.csv(AllSigMarkerWise10e4SNP[AllSigMarkerWise10e4SNP$GWASmodel==i,],paste(pop,i,"model","SigMarkerWise10e4SNP","list.csv",sep="."),quote=F,row.names=F)} #write 4 new files based on GWAS models (all organs)
          for(i in unique(AllSigExpWiseSNP$GWASmodel)){write.csv(AllSigExpWiseSNP[AllSigExpWiseSNP$GWASmodel==i,],paste(pop,i,"model","SigExpWiseSNP","list.csv",sep="."),quote=F,row.names=F)} #write 4 new files based on GWAS models (all organs)
          for(i in unique(AllSigMarkerWiseSNP$organ)){write.csv(AllSigMarkerWiseSNP[AllSigMarkerWiseSNP$organ==i,],paste(pop,i,"organ","SigMarkerWiseSNP","list.csv",sep="."),quote=F,row.names=F)} #write 3 new files based on organs (all models)
          
        }
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        time.taken
        write(time.taken,"timetaken_filterSNP.txt")
        # The SNP will be browsed against the genome browser
        #--http://rs-bt-mccouch4.biotech.cornell.edu:443/cgi-bin/hgGateway?clade=poaceae&org=O.+sativa&db=orySat2
        #--https://snp-seek.irri.org/
      }
      {#|--------Step 3 Plotting (multiple) Manhattan Plots *with only significant SNPs --------| ####
        #Currently I set it to plot manhattan plots using only sig SNPS with p-values <0.001 using the files I saved from Step 2 for EACH GWAS model. If I want to use all SNPs, that can be done, too. Just need to use the full GWAS result files, or the "ManhattanPlot" files I saved from Step 1.
        start.time <- Sys.time()
        writeLines(capture.output(sessionInfo()), "sessionInfo_manhattan.txt") #save session info and check for compatibility later, if needed
        writeLines(capture.output(getwd()), "wd_manhattan.txt")#wd
        {#*---- a.  Load data ----
          # Read file that contains values for manhattan plots; Use ONLY significant SNPs
          EditTextHere <- paste("SigMarkerWiseSNP.list")
          dir(pattern=EditTextHere) #Confirm selection
          file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
        }
        {#*---- b. Making multiple manhattan plots ----
          lapply(file.list,function(x) {
            file <- fread(x) #read the data into R
            x1 <- strsplit(x,"\\.")
            id = x1[[1]][2]
            confirm =  x1[[1]][3]
            population = x1[[1]][1]
            file <-MVP.Report(file[,4:10],
                              plot.type=c("m"),
                              multracks=T, #If this one is true, will get 2 plots:multitracks and multitraits, instead of n plot (for n traits)
                              threshold=c(0.1/700000,1e-3),
                              threshold.lty=c(2,1),
                              threshold.lwd=c(2,2),
                              threshold.col=c("black","black"),
                              amplify=F,
                              ylim=c(3,round(-log10(min(file[,7:10], na.rm=F)))+3),
                              band=3,
                              file.type="pdf",
                              memo=paste0(population,".",id,".",confirm),
                              dpi=5,
                              verbose = T)
            return (file)
          })
        }
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        time.taken
        write(time.taken,"timetaken_manhattan.txt")
      }
      {#|-------- Step 4 Plotting (multiple) QQ plots with all GWAS p-values--------| ####
        #This step can take a long time; only run when needed
        start.time <- Sys.time()
        writeLines(capture.output(sessionInfo()), "sessionInfo_QQ.txt") #save session info and check for compatibility later, if needed
        writeLines(capture.output(getwd()), "wd_QQ.txt")#wd
        {#*---- a.  Load data ----
          # Read file that contains values for manhattan plots
          EditTextHere <- paste(".ManhattanPlot")
          dir(pattern=EditTextHere) #Confirm selection
          file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
        }
        {#*---- b. Making multiple QQ plots from full GWAS results ----
          lapply(file.list,function(x) {
            file <- fread(x) #read the data into R
            file <- rename_with(file, ~ gsub("_whole_", "_", .x, fixed = T))
            file <- rename_with(file, ~ gsub("_root_", "_", .x, fixed = T))
            file <- rename_with(file, ~ gsub("_shoot_", "_", .x, fixed = T))
            file <- rename_with(file,~ gsub("GAPIT.", "", .x, fixed = TRUE))
            file <- rename_with(file,~ gsub(".GWAS.Results", "", .x, fixed = TRUE))
            file <- rename_with(file,~ gsub("Blink.", "", .x, fixed = TRUE))
            file <- rename_with(file,~ gsub("FarmCPU.", "", .x, fixed = TRUE))
            file <- rename_with(file,~ gsub("MLM.", "", .x, fixed = TRUE))
            file <- rename_with(file,~ gsub("GLM.", "", .x, fixed = TRUE))
            
            x1 <- strsplit(x,"\\.")
            GWASmodel = x1[[1]][2]
            population = x1[[1]][1]
            organ = x1[[1]][3]
            file <-MVP.Report(file[,1:7],
                              plot.type=c("q"),
                              multracks=T, #If this one is true, will get 2 plots: overlayed QQ plots and multiple QQ plots on the same plane
                              file.type="pdf",
                              memo=paste0(population,".",GWASmodel,".",organ),
                              dpi=5,
                              verbose = T)
            return (file)
          })
        }
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        time.taken
        write(time.taken,"timetaken_QQ.txt")
      }
    }
    {#* Indica - Specify work directory and population ####
      setwd(paste0(workstation,"/ind"))
      pop = "ind"
      {#|--------Step 1 Merge p-value across traits--------| ####
        start.time <- Sys.time()
        writeLines(capture.output(sessionInfo()), "sessionInfo_mergesnp.txt") #save session info and check for compatibility later, if needed
        writeLines(capture.output(getwd()), "wd_mergesnp.txt")#wd
        {#*1. BLINK root - Merge p-value across traits for Manhattan Plot ----
          {#*---- a. Specify organ to automate file selection ----
            method = "Blink"
            organ = "root"
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
        }
        {#*2. BLINK shoot - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "Blink"
            organ = "shoot"
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
        }
        {#*3. BLINK whole - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "Blink"
            organ = "whole"
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*4.  Farm root - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "FarmCPU"
            
            organ = "root"
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*5. Farm shoot - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "FarmCPU"
            
            
            organ = "shoot"
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*6.  Farm whole - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "FarmCPU"
            
            organ = "whole"
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*7. MLM root - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "MLM"
            
            organ = "root"
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*8.  MLM shoot - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "MLM"
            
            
            organ = "shoot"
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*9. MLM whole - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            
            method = "MLM"
            
            organ = "whole"
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*10.  GLM root - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            
            method = "GLM"
            organ = "root"
            
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*11.  GLM shoot - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            
            method = "GLM"
            
            
            organ = "shoot"
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*12.  GLM whole - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            
            method = "GLM"
            
            organ = "whole"
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
        }
        
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        time.taken
        write(time.taken,"timetaken_mergesnp.txt")
      }
      {#|--------Step 2 Filter significant SNPs--------| ####
        start.time <- Sys.time()
        writeLines(capture.output(sessionInfo()), "sessionInfo_filtersnp.txt") #save session info and check for compatibility later, if needed
        writeLines(capture.output(getwd()), "wd_filtersnp.txt")#wd
        {#*---- a. Prepare data files ----
          # Read file that contains values for manhattan plots
          EditTextHere <- "*.ManhattanPlot" #Edit file name here
          dir(pattern=EditTextHere) #Confirm selection
          file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
          
          # Read csv files for manhattan plots
          allpvalue <- lapply(file.list,function(x) {
            file <- fread(x) #read the data into R
            file <- rename_with(file, ~ gsub("_whole", "_", .x, fixed = T))
            file <- rename_with(file, ~ gsub("_root", "_", .x, fixed = T))
            file <- rename_with(file, ~ gsub("_shoot", "_", .x, fixed = T))
            x1 <- strsplit(x,"\\.")
            organ = x1[[1]][3]
            GWASmodel = x1[[1]][2]
            population = x1[[1]][1]
            file <- cbind (population, GWASmodel, organ, file) # add file name as a new column
            return (file)
          })
          
        }
        {#*---- a extra. For GAPIT results with different file name format; Rename column by removing prefix and sffix ----
          allpvalue<-allpvalue %>%
            map(~ .x %>%
                  rename_with(~ gsub("GAPIT.", "", .x, fixed = TRUE))%>%
                  rename_with(~ gsub(".GWAS.Results", "", .x, fixed = TRUE))%>%
                  rename_with(~ gsub("Blink.", "", .x, fixed = TRUE))%>%
                  rename_with(~ gsub("FarmCPU.", "", .x, fixed = TRUE))%>%
                  rename_with(~ gsub("MLM.", "", .x, fixed = TRUE))%>%
                  rename_with(~ gsub("GLM.", "", .x, fixed = TRUE))
            )
        }
        {#*---- b. Select significant SNPs across N treatments for each organ ----
          # Filter only SNPs with p-value below threshold (Bonferoni Correction = 0.1 for 700,000k SNP)
          SignificantMarkerWiseSNP<- allpvalue %>%
            map(~ .x %>%
                  filter_at(vars(matches('_')), any_vars(.< (0.001))) #p-value < 0.001
            )
          SignificantMarkerWise10e4SNP<- allpvalue %>%
            map(~ .x %>%
                  filter_at(vars(matches('_')), any_vars(.< (0.0001))) #p-value < 0.001
            )
          SignificantExpWiseSNP<- allpvalue %>%
            map(~ .x %>%
                  filter_at(vars(matches('_')), any_vars(.< (0.1/700000))) # Bonferoni correction alpha = 0.1/number of SNPs used.
            )
          
          # Merge all significant SNPs into one data table
          AllSigMarkerWiseSNP = rbindlist(SignificantMarkerWiseSNP,fill=F, use.names=T) #match column names by name, not fill with NAs if missing
          AllSigMarkerWise10e4SNP = rbindlist(SignificantMarkerWise10e4SNP,fill=F, use.names=T ) #match column names by name, not fill with NAs if missing
          AllSigExpWiseSNP = rbindlist(SignificantExpWiseSNP,fill=F, use.names=T ) #match column names by name, not fill with NAs if missing
          #https://www.rdocumentation.org/packages/data.table/versions/1.14.8/topics/rbindlist
          
          write.csv(AllSigMarkerWiseSNP, file=paste(pop,"AllModels","SigMarkerWiseSNP","list.csv",sep="."), row.names=FALSE)
          write.csv(AllSigMarkerWise10e4SNP, file=paste(pop,"AllModels","SigMarkerWise10e4SNP","list.csv",sep="."), row.names=FALSE)
          write.csv(AllSigExpWiseSNP, file=paste(pop,"AllModels","SigExpWiseSNP","list.csv",sep="."), row.names=FALSE)
          
          # Count number of significant SNPs
          AllSigMarkerWiseSNP%>%
            #mutate(newgroup=paste(organ,GWASmodel,sep= "."))%>%
            group_by(population,organ, GWASmodel)%>%
            tally()%>%
            write.csv(., file=paste(pop,"SigMarkerWiseSNP","count.csv",sep="."),quote=F,row.names=F)
          
          AllSigMarkerWise10e4SNP%>%
            #mutate(newgroup=paste(organ,GWASmodel,sep= "."))%>%
            group_by(population,organ, GWASmodel)%>%
            tally()%>%
            write.csv(., file=paste(pop,"SigMarkerWise10e4SNP","count.csv",sep="."),quote=F,row.names=F)
          
          AllSigExpWiseSNP%>%
            #mutate(newgroup=paste(organ,GWASmodel,sep= "."))%>%
            group_by(population,organ, GWASmodel)%>%
            tally()%>%
            write.csv(., file=paste(pop,"SigExpWiseSNP","count.csv",sep="."),quote=F,row.names=F)
          
          # Save a separate csv file for each GWAS model
          for(i in unique(AllSigMarkerWiseSNP$GWASmodel)){write.csv(AllSigMarkerWiseSNP[AllSigMarkerWiseSNP$GWASmodel==i,],paste(pop,i,"model","SigMarkerWiseSNP","list.csv",sep="."),quote=F,row.names=F)} #write 4 new files based on GWAS models (all organs)
          for(i in unique(AllSigMarkerWise10e4SNP$GWASmodel)){write.csv(AllSigMarkerWise10e4SNP[AllSigMarkerWise10e4SNP$GWASmodel==i,],paste(pop,i,"model","SigMarkerWise10e4SNP","list.csv",sep="."),quote=F,row.names=F)} #write 4 new files based on GWAS models (all organs)
          for(i in unique(AllSigExpWiseSNP$GWASmodel)){write.csv(AllSigExpWiseSNP[AllSigExpWiseSNP$GWASmodel==i,],paste(pop,i,"model","SigExpWiseSNP","list.csv",sep="."),quote=F,row.names=F)} #write 4 new files based on GWAS models (all organs)
          for(i in unique(AllSigMarkerWiseSNP$organ)){write.csv(AllSigMarkerWiseSNP[AllSigMarkerWiseSNP$organ==i,],paste(pop,i,"organ","SigMarkerWiseSNP","list.csv",sep="."),quote=F,row.names=F)} #write 3 new files based on organs (all models)
          
        }
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        time.taken
        write(time.taken,"timetaken_filterSNP.txt")
        # The SNP will be browsed against the genome browser
        #--http://rs-bt-mccouch4.biotech.cornell.edu:443/cgi-bin/hgGateway?clade=poaceae&org=O.+sativa&db=orySat2
        #--https://snp-seek.irri.org/
      }
      {#|--------Step 3 Plotting (multiple) Manhattan Plots *with only significant SNPs --------| ####
        #Currently I set it to plot manhattan plots using only sig SNPS with p-values <0.001 using the files I saved from Step 2 for EACH GWAS model. If I want to use all SNPs, that can be done, too. Just need to use the full GWAS result files, or the "ManhattanPlot" files I saved from Step 1.
        start.time <- Sys.time()
        writeLines(capture.output(sessionInfo()), "sessionInfo_manhattan.txt") #save session info and check for compatibility later, if needed
        writeLines(capture.output(getwd()), "wd_manhattan.txt")#wd
        {#*---- a.  Load data ----
          # Read file that contains values for manhattan plots; Use ONLY significant SNPs
          EditTextHere <- paste("SigMarkerWiseSNP.list")
          dir(pattern=EditTextHere) #Confirm selection
          file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
        }
        {#*---- b. Making multiple manhattan plots ----
          lapply(file.list,function(x) {
            file <- fread(x) #read the data into R
            x1 <- strsplit(x,"\\.")
            id = x1[[1]][2]
            confirm =  x1[[1]][3]
            population = x1[[1]][1]
            file <-MVP.Report(file[,4:10],
                              plot.type=c("m"),
                              multracks=T, #If this one is true, will get 2 plots:multitracks and multitraits, instead of n plot (for n traits)
                              threshold=c(0.1/700000,1e-3),
                              threshold.lty=c(2,1),
                              threshold.lwd=c(2,2),
                              threshold.col=c("black","black"),
                              amplify=F,
                              ylim=c(3,round(-log10(min(file[,7:10], na.rm=F)))+3),
                              band=3,
                              file.type="pdf",
                              memo=paste0(population,".",id,".",confirm),
                              dpi=5,
                              verbose = T)
            return (file)
          })
        }
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        time.taken
        write(time.taken,"timetaken_manhattan.txt")
      }
      {#|-------- Step 4 Plotting (multiple) QQ plots with all GWAS p-values--------| ####
        #This step can take a long time; only run when needed
        start.time <- Sys.time()
        writeLines(capture.output(sessionInfo()), "sessionInfo_QQ.txt") #save session info and check for compatibility later, if needed
        writeLines(capture.output(getwd()), "wd_QQ.txt")#wd
        {#*---- a.  Load data ----
          # Read file that contains values for manhattan plots
          EditTextHere <- paste(".ManhattanPlot")
          dir(pattern=EditTextHere) #Confirm selection
          file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
        }
        {#*---- b. Making multiple QQ plots from full GWAS results ----
          lapply(file.list,function(x) {
            file <- fread(x) #read the data into R
            file <- rename_with(file, ~ gsub("_whole_", "_", .x, fixed = T))
            file <- rename_with(file, ~ gsub("_root_", "_", .x, fixed = T))
            file <- rename_with(file, ~ gsub("_shoot_", "_", .x, fixed = T))
            file <- rename_with(file,~ gsub("GAPIT.", "", .x, fixed = TRUE))
            file <- rename_with(file,~ gsub(".GWAS.Results", "", .x, fixed = TRUE))
            file <- rename_with(file,~ gsub("Blink.", "", .x, fixed = TRUE))
            file <- rename_with(file,~ gsub("FarmCPU.", "", .x, fixed = TRUE))
            file <- rename_with(file,~ gsub("MLM.", "", .x, fixed = TRUE))
            file <- rename_with(file,~ gsub("GLM.", "", .x, fixed = TRUE))
            
            x1 <- strsplit(x,"\\.")
            GWASmodel = x1[[1]][2]
            population = x1[[1]][1]
            organ = x1[[1]][3]
            file <-MVP.Report(file[,1:7],
                              plot.type=c("q"),
                              multracks=T, #If this one is true, will get 2 plots: overlayed QQ plots and multiple QQ plots on the same plane
                              file.type="pdf",
                              memo=paste0(population,".",GWASmodel,".",organ),
                              dpi=5,
                              verbose = T)
            return (file)
          })
        }
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        time.taken
        write(time.taken,"timetaken_QQ.txt")
      }
    }
    {#* TEJ- Specify work directory and population ####
      setwd(paste0(workstation,"/tej"))
      pop = "tej"
      {#|--------Step 1 Merge p-value across traits--------| ####
        start.time <- Sys.time()
        writeLines(capture.output(sessionInfo()), "sessionInfo_mergesnp.txt") #save session info and check for compatibility later, if needed
        writeLines(capture.output(getwd()), "wd_mergesnp.txt")#wd
        {#*1. BLINK root - Merge p-value across traits for Manhattan Plot ----
          {#*---- a. Specify organ to automate file selection ----
            method = "Blink"
            organ = "root"
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
        }
        {#*2. BLINK shoot - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "Blink"
            organ = "shoot"
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
        }
        {#*3. BLINK whole - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "Blink"
            organ = "whole"
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*4.  Farm root - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "FarmCPU"
            
            organ = "root"
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*5. Farm shoot - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "FarmCPU"
            
            
            organ = "shoot"
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*6.  Farm whole - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "FarmCPU"
            
            organ = "whole"
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*7. MLM root - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "MLM"
            
            organ = "root"
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*8.  MLM shoot - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "MLM"
            
            
            organ = "shoot"
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*9. MLM whole - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            
            method = "MLM"
            
            organ = "whole"
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*10.  GLM root - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            
            method = "GLM"
            organ = "root"
            
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*11.  GLM shoot - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            
            method = "GLM"
            
            
            organ = "shoot"
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*12.  GLM whole - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            
            method = "GLM"
            
            organ = "whole"
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
        }
        
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        time.taken
        write(time.taken,"timetaken_mergesnp.txt")
      }
      {#|--------Step 2 Filter significant SNPs--------| ####
        start.time <- Sys.time()
        writeLines(capture.output(sessionInfo()), "sessionInfo_filtersnp.txt") #save session info and check for compatibility later, if needed
        writeLines(capture.output(getwd()), "wd_filtersnp.txt")#wd
        {#*---- a. Prepare data files ----
          # Read file that contains values for manhattan plots
          EditTextHere <- "*.ManhattanPlot" #Edit file name here
          dir(pattern=EditTextHere) #Confirm selection
          file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
          
          # Read csv files for manhattan plots
          allpvalue <- lapply(file.list,function(x) {
            file <- fread(x) #read the data into R
            file <- rename_with(file, ~ gsub("_whole", "_", .x, fixed = T))
            file <- rename_with(file, ~ gsub("_root", "_", .x, fixed = T))
            file <- rename_with(file, ~ gsub("_shoot", "_", .x, fixed = T))
            x1 <- strsplit(x,"\\.")
            organ = x1[[1]][3]
            GWASmodel = x1[[1]][2]
            population = x1[[1]][1]
            file <- cbind (population, GWASmodel, organ, file) # add file name as a new column
            return (file)
          })
          
        }
        {#*---- a extra. For GAPIT results with different file name format; Rename column by removing prefix and sffix ----
          allpvalue<-allpvalue %>%
            map(~ .x %>%
                  rename_with(~ gsub("GAPIT.", "", .x, fixed = TRUE))%>%
                  rename_with(~ gsub(".GWAS.Results", "", .x, fixed = TRUE))%>%
                  rename_with(~ gsub("Blink.", "", .x, fixed = TRUE))%>%
                  rename_with(~ gsub("FarmCPU.", "", .x, fixed = TRUE))%>%
                  rename_with(~ gsub("MLM.", "", .x, fixed = TRUE))%>%
                  rename_with(~ gsub("GLM.", "", .x, fixed = TRUE))
            )
        }
        {#*---- b. Select significant SNPs across N treatments for each organ ----
          # Filter only SNPs with p-value below threshold (Bonferoni Correction = 0.1 for 700,000k SNP)
          SignificantMarkerWiseSNP<- allpvalue %>%
            map(~ .x %>%
                  filter_at(vars(matches('_')), any_vars(.< (0.001))) #p-value < 0.001
            )
          SignificantMarkerWise10e4SNP<- allpvalue %>%
            map(~ .x %>%
                  filter_at(vars(matches('_')), any_vars(.< (0.0001))) #p-value < 0.001
            )
          SignificantExpWiseSNP<- allpvalue %>%
            map(~ .x %>%
                  filter_at(vars(matches('_')), any_vars(.< (0.1/700000))) # Bonferoni correction alpha = 0.1/number of SNPs used.
            )
          
          # Merge all significant SNPs into one data table
          AllSigMarkerWiseSNP = rbindlist(SignificantMarkerWiseSNP,fill=F, use.names=T) #match column names by name, not fill with NAs if missing
          AllSigMarkerWise10e4SNP = rbindlist(SignificantMarkerWise10e4SNP,fill=F, use.names=T ) #match column names by name, not fill with NAs if missing
          AllSigExpWiseSNP = rbindlist(SignificantExpWiseSNP,fill=F, use.names=T ) #match column names by name, not fill with NAs if missing
          #https://www.rdocumentation.org/packages/data.table/versions/1.14.8/topics/rbindlist
          
          write.csv(AllSigMarkerWiseSNP, file=paste(pop,"AllModels","SigMarkerWiseSNP","list.csv",sep="."), row.names=FALSE)
          write.csv(AllSigMarkerWise10e4SNP, file=paste(pop,"AllModels","SigMarkerWise10e4SNP","list.csv",sep="."), row.names=FALSE)
          write.csv(AllSigExpWiseSNP, file=paste(pop,"AllModels","SigExpWiseSNP","list.csv",sep="."), row.names=FALSE)
          
          # Count number of significant SNPs
          AllSigMarkerWiseSNP%>%
            #mutate(newgroup=paste(organ,GWASmodel,sep= "."))%>%
            group_by(population,organ, GWASmodel)%>%
            tally()%>%
            write.csv(., file=paste(pop,"SigMarkerWiseSNP","count.csv",sep="."),quote=F,row.names=F)
          
          AllSigMarkerWise10e4SNP%>%
            #mutate(newgroup=paste(organ,GWASmodel,sep= "."))%>%
            group_by(population,organ, GWASmodel)%>%
            tally()%>%
            write.csv(., file=paste(pop,"SigMarkerWise10e4SNP","count.csv",sep="."),quote=F,row.names=F)
          
          AllSigExpWiseSNP%>%
            #mutate(newgroup=paste(organ,GWASmodel,sep= "."))%>%
            group_by(population,organ, GWASmodel)%>%
            tally()%>%
            write.csv(., file=paste(pop,"SigExpWiseSNP","count.csv",sep="."),quote=F,row.names=F)
          
          # Save a separate csv file for each GWAS model
          for(i in unique(AllSigMarkerWiseSNP$GWASmodel)){write.csv(AllSigMarkerWiseSNP[AllSigMarkerWiseSNP$GWASmodel==i,],paste(pop,i,"model","SigMarkerWiseSNP","list.csv",sep="."),quote=F,row.names=F)} #write 4 new files based on GWAS models (all organs)
          for(i in unique(AllSigMarkerWise10e4SNP$GWASmodel)){write.csv(AllSigMarkerWise10e4SNP[AllSigMarkerWise10e4SNP$GWASmodel==i,],paste(pop,i,"model","SigMarkerWise10e4SNP","list.csv",sep="."),quote=F,row.names=F)} #write 4 new files based on GWAS models (all organs)
          for(i in unique(AllSigExpWiseSNP$GWASmodel)){write.csv(AllSigExpWiseSNP[AllSigExpWiseSNP$GWASmodel==i,],paste(pop,i,"model","SigExpWiseSNP","list.csv",sep="."),quote=F,row.names=F)} #write 4 new files based on GWAS models (all organs)
          for(i in unique(AllSigMarkerWiseSNP$organ)){write.csv(AllSigMarkerWiseSNP[AllSigMarkerWiseSNP$organ==i,],paste(pop,i,"organ","SigMarkerWiseSNP","list.csv",sep="."),quote=F,row.names=F)} #write 3 new files based on organs (all models)
          
        }
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        time.taken
        write(time.taken,"timetaken_filterSNP.txt")
        # The SNP will be browsed against the genome browser
        #--http://rs-bt-mccouch4.biotech.cornell.edu:443/cgi-bin/hgGateway?clade=poaceae&org=O.+sativa&db=orySat2
        #--https://snp-seek.irri.org/
      }
      {#|--------Step 3 Plotting (multiple) Manhattan Plots *with only significant SNPs --------| ####
        #Currently I set it to plot manhattan plots using only sig SNPS with p-values <0.001 using the files I saved from Step 2 for EACH GWAS model. If I want to use all SNPs, that can be done, too. Just need to use the full GWAS result files, or the "ManhattanPlot" files I saved from Step 1.
        start.time <- Sys.time()
        writeLines(capture.output(sessionInfo()), "sessionInfo_manhattan.txt") #save session info and check for compatibility later, if needed
        writeLines(capture.output(getwd()), "wd_manhattan.txt")#wd
        {#*---- a.  Load data ----
          # Read file that contains values for manhattan plots; Use ONLY significant SNPs
          EditTextHere <- paste("SigMarkerWiseSNP.list")
          dir(pattern=EditTextHere) #Confirm selection
          file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
        }
        {#*---- b. Making multiple manhattan plots ----
          lapply(file.list,function(x) {
            file <- fread(x) #read the data into R
            x1 <- strsplit(x,"\\.")
            id = x1[[1]][2]
            confirm =  x1[[1]][3]
            population = x1[[1]][1]
            file <-MVP.Report(file[,4:10],
                              plot.type=c("m"),
                              multracks=T, #If this one is true, will get 2 plots:multitracks and multitraits, instead of n plot (for n traits)
                              threshold=c(0.1/700000,1e-3),
                              threshold.lty=c(2,1),
                              threshold.lwd=c(2,2),
                              threshold.col=c("black","black"),
                              amplify=F,
                              ylim=c(3,round(-log10(min(file[,7:10], na.rm=F)))+3),
                              band=3,
                              file.type="pdf",
                              memo=paste0(population,".",id,".",confirm),
                              dpi=5,
                              verbose = T)
            return (file)
          })
        }
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        time.taken
        write(time.taken,"timetaken_manhattan.txt")
      }
      {#|-------- Step 4 Plotting (multiple) QQ plots with all GWAS p-values--------| ####
        #This step can take a long time; only run when needed
        start.time <- Sys.time()
        writeLines(capture.output(sessionInfo()), "sessionInfo_QQ.txt") #save session info and check for compatibility later, if needed
        writeLines(capture.output(getwd()), "wd_QQ.txt")#wd
        {#*---- a.  Load data ----
          # Read file that contains values for manhattan plots
          EditTextHere <- paste(".ManhattanPlot")
          dir(pattern=EditTextHere) #Confirm selection
          file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
        }
        {#*---- b. Making multiple QQ plots from full GWAS results ----
          lapply(file.list,function(x) {
            file <- fread(x) #read the data into R
            file <- rename_with(file, ~ gsub("_whole_", "_", .x, fixed = T))
            file <- rename_with(file, ~ gsub("_root_", "_", .x, fixed = T))
            file <- rename_with(file, ~ gsub("_shoot_", "_", .x, fixed = T))
            file <- rename_with(file,~ gsub("GAPIT.", "", .x, fixed = TRUE))
            file <- rename_with(file,~ gsub(".GWAS.Results", "", .x, fixed = TRUE))
            file <- rename_with(file,~ gsub("Blink.", "", .x, fixed = TRUE))
            file <- rename_with(file,~ gsub("FarmCPU.", "", .x, fixed = TRUE))
            file <- rename_with(file,~ gsub("MLM.", "", .x, fixed = TRUE))
            file <- rename_with(file,~ gsub("GLM.", "", .x, fixed = TRUE))
            
            x1 <- strsplit(x,"\\.")
            GWASmodel = x1[[1]][2]
            population = x1[[1]][1]
            organ = x1[[1]][3]
            file <-MVP.Report(file[,1:7],
                              plot.type=c("q"),
                              multracks=T, #If this one is true, will get 2 plots: overlayed QQ plots and multiple QQ plots on the same plane
                              file.type="pdf",
                              memo=paste0(population,".",GWASmodel,".",organ),
                              dpi=5,
                              verbose = T)
            return (file)
          })
        }
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        time.taken
        write(time.taken,"timetaken_QQ.txt")
      }
    }
    {#* TRJ- Specify work directory and population ####
      setwd(paste0(workstation,"/trj"))
      pop = "trj"
      {#|--------Step 1 Merge p-value across traits--------| ####
        start.time <- Sys.time()
        writeLines(capture.output(sessionInfo()), "sessionInfo_mergesnp.txt") #save session info and check for compatibility later, if needed
        writeLines(capture.output(getwd()), "wd_mergesnp.txt")#wd
        {#*1. BLINK root - Merge p-value across traits for Manhattan Plot ----
          {#*---- a. Specify organ to automate file selection ----
            method = "Blink"
            organ = "root"
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
        }
        {#*2. BLINK shoot - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "Blink"
            organ = "shoot"
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
        }
        {#*3. BLINK whole - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "Blink"
            organ = "whole"
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*4.  Farm root - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "FarmCPU"
            
            organ = "root"
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*5. Farm shoot - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "FarmCPU"
            
            
            organ = "shoot"
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*6.  Farm whole - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "FarmCPU"
            
            organ = "whole"
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*7. MLM root - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "MLM"
            
            organ = "root"
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*8.  MLM shoot - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "MLM"
            
            
            organ = "shoot"
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*9. MLM whole - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            
            method = "MLM"
            
            organ = "whole"
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*10.  GLM root - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            
            method = "GLM"
            organ = "root"
            
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*11.  GLM shoot - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            
            method = "GLM"
            
            
            organ = "shoot"
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*12.  GLM whole - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            
            method = "GLM"
            
            organ = "whole"
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
        }
        
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        time.taken
        write(time.taken,"timetaken_mergesnp.txt")
      }
      {#|--------Step 2 Filter significant SNPs--------| ####
        start.time <- Sys.time()
        writeLines(capture.output(sessionInfo()), "sessionInfo_filtersnp.txt") #save session info and check for compatibility later, if needed
        writeLines(capture.output(getwd()), "wd_filtersnp.txt")#wd
        {#*---- a. Prepare data files ----
          # Read file that contains values for manhattan plots
          EditTextHere <- "*.ManhattanPlot" #Edit file name here
          dir(pattern=EditTextHere) #Confirm selection
          file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
          
          # Read csv files for manhattan plots
          allpvalue <- lapply(file.list,function(x) {
            file <- fread(x) #read the data into R
            file <- rename_with(file, ~ gsub("_whole", "_", .x, fixed = T))
            file <- rename_with(file, ~ gsub("_root", "_", .x, fixed = T))
            file <- rename_with(file, ~ gsub("_shoot", "_", .x, fixed = T))
            x1 <- strsplit(x,"\\.")
            organ = x1[[1]][3]
            GWASmodel = x1[[1]][2]
            population = x1[[1]][1]
            file <- cbind (population, GWASmodel, organ, file) # add file name as a new column
            return (file)
          })
          
        }
        {#*---- a extra. For GAPIT results with different file name format; Rename column by removing prefix and sffix ----
          allpvalue<-allpvalue %>%
            map(~ .x %>%
                  rename_with(~ gsub("GAPIT.", "", .x, fixed = TRUE))%>%
                  rename_with(~ gsub(".GWAS.Results", "", .x, fixed = TRUE))%>%
                  rename_with(~ gsub("Blink.", "", .x, fixed = TRUE))%>%
                  rename_with(~ gsub("FarmCPU.", "", .x, fixed = TRUE))%>%
                  rename_with(~ gsub("MLM.", "", .x, fixed = TRUE))%>%
                  rename_with(~ gsub("GLM.", "", .x, fixed = TRUE))
            )
        }
        {#*---- b. Select significant SNPs across N treatments for each organ ----
          # Filter only SNPs with p-value below threshold (Bonferoni Correction = 0.1 for 700,000k SNP)
          SignificantMarkerWiseSNP<- allpvalue %>%
            map(~ .x %>%
                  filter_at(vars(matches('_')), any_vars(.< (0.001))) #p-value < 0.001
            )
          SignificantMarkerWise10e4SNP<- allpvalue %>%
            map(~ .x %>%
                  filter_at(vars(matches('_')), any_vars(.< (0.0001))) #p-value < 0.001
            )
          SignificantExpWiseSNP<- allpvalue %>%
            map(~ .x %>%
                  filter_at(vars(matches('_')), any_vars(.< (0.1/700000))) # Bonferoni correction alpha = 0.1/number of SNPs used.
            )
          
          # Merge all significant SNPs into one data table
          AllSigMarkerWiseSNP = rbindlist(SignificantMarkerWiseSNP,fill=F, use.names=T) #match column names by name, not fill with NAs if missing
          AllSigMarkerWise10e4SNP = rbindlist(SignificantMarkerWise10e4SNP,fill=F, use.names=T ) #match column names by name, not fill with NAs if missing
          AllSigExpWiseSNP = rbindlist(SignificantExpWiseSNP,fill=F, use.names=T ) #match column names by name, not fill with NAs if missing
          #https://www.rdocumentation.org/packages/data.table/versions/1.14.8/topics/rbindlist
          
          write.csv(AllSigMarkerWiseSNP, file=paste(pop,"AllModels","SigMarkerWiseSNP","list.csv",sep="."), row.names=FALSE)
          write.csv(AllSigMarkerWise10e4SNP, file=paste(pop,"AllModels","SigMarkerWise10e4SNP","list.csv",sep="."), row.names=FALSE)
          write.csv(AllSigExpWiseSNP, file=paste(pop,"AllModels","SigExpWiseSNP","list.csv",sep="."), row.names=FALSE)
          
          # Count number of significant SNPs
          AllSigMarkerWiseSNP%>%
            #mutate(newgroup=paste(organ,GWASmodel,sep= "."))%>%
            group_by(population,organ, GWASmodel)%>%
            tally()%>%
            write.csv(., file=paste(pop,"SigMarkerWiseSNP","count.csv",sep="."),quote=F,row.names=F)
          
          AllSigMarkerWise10e4SNP%>%
            #mutate(newgroup=paste(organ,GWASmodel,sep= "."))%>%
            group_by(population,organ, GWASmodel)%>%
            tally()%>%
            write.csv(., file=paste(pop,"SigMarkerWise10e4SNP","count.csv",sep="."),quote=F,row.names=F)
          
          AllSigExpWiseSNP%>%
            #mutate(newgroup=paste(organ,GWASmodel,sep= "."))%>%
            group_by(population,organ, GWASmodel)%>%
            tally()%>%
            write.csv(., file=paste(pop,"SigExpWiseSNP","count.csv",sep="."),quote=F,row.names=F)
          
          # Save a separate csv file for each GWAS model
          for(i in unique(AllSigMarkerWiseSNP$GWASmodel)){write.csv(AllSigMarkerWiseSNP[AllSigMarkerWiseSNP$GWASmodel==i,],paste(pop,i,"model","SigMarkerWiseSNP","list.csv",sep="."),quote=F,row.names=F)} #write 4 new files based on GWAS models (all organs)
          for(i in unique(AllSigMarkerWise10e4SNP$GWASmodel)){write.csv(AllSigMarkerWise10e4SNP[AllSigMarkerWise10e4SNP$GWASmodel==i,],paste(pop,i,"model","SigMarkerWise10e4SNP","list.csv",sep="."),quote=F,row.names=F)} #write 4 new files based on GWAS models (all organs)
          for(i in unique(AllSigExpWiseSNP$GWASmodel)){write.csv(AllSigExpWiseSNP[AllSigExpWiseSNP$GWASmodel==i,],paste(pop,i,"model","SigExpWiseSNP","list.csv",sep="."),quote=F,row.names=F)} #write 4 new files based on GWAS models (all organs)
          for(i in unique(AllSigMarkerWiseSNP$organ)){write.csv(AllSigMarkerWiseSNP[AllSigMarkerWiseSNP$organ==i,],paste(pop,i,"organ","SigMarkerWiseSNP","list.csv",sep="."),quote=F,row.names=F)} #write 3 new files based on organs (all models)
          
        }
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        time.taken
        write(time.taken,"timetaken_filterSNP.txt")
        # The SNP will be browsed against the genome browser
        #--http://rs-bt-mccouch4.biotech.cornell.edu:443/cgi-bin/hgGateway?clade=poaceae&org=O.+sativa&db=orySat2
        #--https://snp-seek.irri.org/
      }
      {#|--------Step 3 Plotting (multiple) Manhattan Plots *with only significant SNPs --------| ####
        #Currently I set it to plot manhattan plots using only sig SNPS with p-values <0.001 using the files I saved from Step 2 for EACH GWAS model. If I want to use all SNPs, that can be done, too. Just need to use the full GWAS result files, or the "ManhattanPlot" files I saved from Step 1.
        start.time <- Sys.time()
        writeLines(capture.output(sessionInfo()), "sessionInfo_manhattan.txt") #save session info and check for compatibility later, if needed
        writeLines(capture.output(getwd()), "wd_manhattan.txt")#wd
        {#*---- a.  Load data ----
          # Read file that contains values for manhattan plots; Use ONLY significant SNPs
          EditTextHere <- paste("SigMarkerWiseSNP.list")
          dir(pattern=EditTextHere) #Confirm selection
          file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
        }
        {#*---- b. Making multiple manhattan plots ----
          lapply(file.list,function(x) {
            file <- fread(x) #read the data into R
            x1 <- strsplit(x,"\\.")
            id = x1[[1]][2]
            confirm =  x1[[1]][3]
            population = x1[[1]][1]
            file <-MVP.Report(file[,4:10],
                              plot.type=c("m"),
                              multracks=T, #If this one is true, will get 2 plots:multitracks and multitraits, instead of n plot (for n traits)
                              threshold=c(0.1/700000,1e-3),
                              threshold.lty=c(2,1),
                              threshold.lwd=c(2,2),
                              threshold.col=c("black","black"),
                              amplify=F,
                              ylim=c(3,round(-log10(min(file[,7:10], na.rm=F)))+3),
                              band=3,
                              file.type="pdf",
                              memo=paste0(population,".",id,".",confirm),
                              dpi=5,
                              verbose = T)
            return (file)
          })
        }
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        time.taken
        write(time.taken,"timetaken_manhattan.txt")
      }
      {#|-------- Step 4 Plotting (multiple) QQ plots with all GWAS p-values--------| ####
        #This step can take a long time; only run when needed
        start.time <- Sys.time()
        writeLines(capture.output(sessionInfo()), "sessionInfo_QQ.txt") #save session info and check for compatibility later, if needed
        writeLines(capture.output(getwd()), "wd_QQ.txt")#wd
        {#*---- a.  Load data ----
          # Read file that contains values for manhattan plots
          EditTextHere <- paste(".ManhattanPlot")
          dir(pattern=EditTextHere) #Confirm selection
          file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
        }
        {#*---- b. Making multiple QQ plots from full GWAS results ----
          lapply(file.list,function(x) {
            file <- fread(x) #read the data into R
            file <- rename_with(file, ~ gsub("_whole_", "_", .x, fixed = T))
            file <- rename_with(file, ~ gsub("_root_", "_", .x, fixed = T))
            file <- rename_with(file, ~ gsub("_shoot_", "_", .x, fixed = T))
            file <- rename_with(file,~ gsub("GAPIT.", "", .x, fixed = TRUE))
            file <- rename_with(file,~ gsub(".GWAS.Results", "", .x, fixed = TRUE))
            file <- rename_with(file,~ gsub("Blink.", "", .x, fixed = TRUE))
            file <- rename_with(file,~ gsub("FarmCPU.", "", .x, fixed = TRUE))
            file <- rename_with(file,~ gsub("MLM.", "", .x, fixed = TRUE))
            file <- rename_with(file,~ gsub("GLM.", "", .x, fixed = TRUE))
            
            x1 <- strsplit(x,"\\.")
            GWASmodel = x1[[1]][2]
            population = x1[[1]][1]
            organ = x1[[1]][3]
            file <-MVP.Report(file[,1:7],
                              plot.type=c("q"),
                              multracks=T, #If this one is true, will get 2 plots: overlayed QQ plots and multiple QQ plots on the same plane
                              file.type="pdf",
                              memo=paste0(population,".",GWASmodel,".",organ),
                              dpi=5,
                              verbose = T)
            return (file)
          })
        }
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        time.taken
        write(time.taken,"timetaken_QQ.txt")
      }
    }
    {#* Var Ind, pc0 - Specify work directory and population ####
      setwd(paste0(workstation,"/varind"))
      pop = "varind-pc0"
      {#|--------Step 1 Merge p-value across traits--------| ####
        start.time <- Sys.time()
        writeLines(capture.output(sessionInfo()), "sessionInfo_mergesnp.txt") #save session info and check for compatibility later, if needed
        writeLines(capture.output(getwd()), "wd_mergesnp.txt")#wd
        {#*1. BLINK root - Merge p-value across traits for Manhattan Plot ----
          {#*---- a. Specify organ to automate file selection ----
            method = "Blink"
            organ = "root"
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
        }
        {#*2. BLINK shoot - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "Blink"
            organ = "shoot"
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
        }
        {#*3. BLINK whole - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "Blink"
            organ = "whole"
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*4.  Farm root - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "FarmCPU"
            
            organ = "root"
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*5. Farm shoot - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "FarmCPU"
            
            
            organ = "shoot"
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*6.  Farm whole - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "FarmCPU"
            
            organ = "whole"
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*7. MLM root - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "MLM"
            
            organ = "root"
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*8.  MLM shoot - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "MLM"
            
            
            organ = "shoot"
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*9. MLM whole - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            
            method = "MLM"
            
            organ = "whole"
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*10.  GLM root - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            
            method = "GLM"
            organ = "root"
            
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*11.  GLM shoot - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            
            method = "GLM"
            
            
            organ = "shoot"
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*12.  GLM whole - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            
            method = "GLM"
            
            organ = "whole"
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
        }
        
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        time.taken
        write(time.taken,"timetaken_mergesnp.txt")
      }
      {#|--------Step 2 Filter significant SNPs--------| ####
        start.time <- Sys.time()
        writeLines(capture.output(sessionInfo()), "sessionInfo_filtersnp.txt") #save session info and check for compatibility later, if needed
        writeLines(capture.output(getwd()), "wd_filtersnp.txt")#wd
        {#*---- a. Prepare data files ----
          # Read file that contains values for manhattan plots
          EditTextHere <- "*.ManhattanPlot" #Edit file name here
          dir(pattern=EditTextHere) #Confirm selection
          file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
          
          # Read csv files for manhattan plots
          allpvalue <- lapply(file.list,function(x) {
            file <- fread(x) #read the data into R
            file <- rename_with(file, ~ gsub("_whole", "_", .x, fixed = T))
            file <- rename_with(file, ~ gsub("_root", "_", .x, fixed = T))
            file <- rename_with(file, ~ gsub("_shoot", "_", .x, fixed = T))
            x1 <- strsplit(x,"\\.")
            organ = x1[[1]][3]
            GWASmodel = x1[[1]][2]
            population = x1[[1]][1]
            file <- cbind (population, GWASmodel, organ, file) # add file name as a new column
            return (file)
          })
          
        }
        {#*---- a extra. For GAPIT results with different file name format; Rename column by removing prefix and sffix ----
          allpvalue<-allpvalue %>%
            map(~ .x %>%
                  rename_with(~ gsub("GAPIT.", "", .x, fixed = TRUE))%>%
                  rename_with(~ gsub(".GWAS.Results", "", .x, fixed = TRUE))%>%
                  rename_with(~ gsub("Blink.", "", .x, fixed = TRUE))%>%
                  rename_with(~ gsub("FarmCPU.", "", .x, fixed = TRUE))%>%
                  rename_with(~ gsub("MLM.", "", .x, fixed = TRUE))%>%
                  rename_with(~ gsub("GLM.", "", .x, fixed = TRUE))
            )
        }
        {#*---- b. Select significant SNPs across N treatments for each organ ----
          # Filter only SNPs with p-value below threshold (Bonferoni Correction = 0.1 for 700,000k SNP)
          SignificantMarkerWiseSNP<- allpvalue %>%
            map(~ .x %>%
                  filter_at(vars(matches('_')), any_vars(.< (0.001))) #p-value < 0.001
            )
          SignificantMarkerWise10e4SNP<- allpvalue %>%
            map(~ .x %>%
                  filter_at(vars(matches('_')), any_vars(.< (0.0001))) #p-value < 0.001
            )
          SignificantExpWiseSNP<- allpvalue %>%
            map(~ .x %>%
                  filter_at(vars(matches('_')), any_vars(.< (0.1/700000))) # Bonferoni correction alpha = 0.1/number of SNPs used.
            )
          
          # Merge all significant SNPs into one data table
          AllSigMarkerWiseSNP = rbindlist(SignificantMarkerWiseSNP,fill=F, use.names=T) #match column names by name, not fill with NAs if missing
          AllSigMarkerWise10e4SNP = rbindlist(SignificantMarkerWise10e4SNP,fill=F, use.names=T ) #match column names by name, not fill with NAs if missing
          AllSigExpWiseSNP = rbindlist(SignificantExpWiseSNP,fill=F, use.names=T ) #match column names by name, not fill with NAs if missing
          #https://www.rdocumentation.org/packages/data.table/versions/1.14.8/topics/rbindlist
          
          write.csv(AllSigMarkerWiseSNP, file=paste(pop,"AllModels","SigMarkerWiseSNP","list.csv",sep="."), row.names=FALSE)
          write.csv(AllSigMarkerWise10e4SNP, file=paste(pop,"AllModels","SigMarkerWise10e4SNP","list.csv",sep="."), row.names=FALSE)
          write.csv(AllSigExpWiseSNP, file=paste(pop,"AllModels","SigExpWiseSNP","list.csv",sep="."), row.names=FALSE)
          
          # Count number of significant SNPs
          AllSigMarkerWiseSNP%>%
            #mutate(newgroup=paste(organ,GWASmodel,sep= "."))%>%
            group_by(population,organ, GWASmodel)%>%
            tally()%>%
            write.csv(., file=paste(pop,"SigMarkerWiseSNP","count.csv",sep="."),quote=F,row.names=F)
          
          AllSigMarkerWise10e4SNP%>%
            #mutate(newgroup=paste(organ,GWASmodel,sep= "."))%>%
            group_by(population,organ, GWASmodel)%>%
            tally()%>%
            write.csv(., file=paste(pop,"SigMarkerWise10e4SNP","count.csv",sep="."),quote=F,row.names=F)
          
          AllSigExpWiseSNP%>%
            #mutate(newgroup=paste(organ,GWASmodel,sep= "."))%>%
            group_by(population,organ, GWASmodel)%>%
            tally()%>%
            write.csv(., file=paste(pop,"SigExpWiseSNP","count.csv",sep="."),quote=F,row.names=F)
          
          # Save a separate csv file for each GWAS model
          for(i in unique(AllSigMarkerWiseSNP$GWASmodel)){write.csv(AllSigMarkerWiseSNP[AllSigMarkerWiseSNP$GWASmodel==i,],paste(pop,i,"model","SigMarkerWiseSNP","list.csv",sep="."),quote=F,row.names=F)} #write 4 new files based on GWAS models (all organs)
          for(i in unique(AllSigMarkerWise10e4SNP$GWASmodel)){write.csv(AllSigMarkerWise10e4SNP[AllSigMarkerWise10e4SNP$GWASmodel==i,],paste(pop,i,"model","SigMarkerWise10e4SNP","list.csv",sep="."),quote=F,row.names=F)} #write 4 new files based on GWAS models (all organs)
          for(i in unique(AllSigExpWiseSNP$GWASmodel)){write.csv(AllSigExpWiseSNP[AllSigExpWiseSNP$GWASmodel==i,],paste(pop,i,"model","SigExpWiseSNP","list.csv",sep="."),quote=F,row.names=F)} #write 4 new files based on GWAS models (all organs)
          for(i in unique(AllSigMarkerWiseSNP$organ)){write.csv(AllSigMarkerWiseSNP[AllSigMarkerWiseSNP$organ==i,],paste(pop,i,"organ","SigMarkerWiseSNP","list.csv",sep="."),quote=F,row.names=F)} #write 3 new files based on organs (all models)
          
        }
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        time.taken
        write(time.taken,"timetaken_filterSNP.txt")
        # The SNP will be browsed against the genome browser
        #--http://rs-bt-mccouch4.biotech.cornell.edu:443/cgi-bin/hgGateway?clade=poaceae&org=O.+sativa&db=orySat2
        #--https://snp-seek.irri.org/
      }
      {#|--------Step 3 Plotting (multiple) Manhattan Plots *with only significant SNPs --------| ####
        #Currently I set it to plot manhattan plots using only sig SNPS with p-values <0.001 using the files I saved from Step 2 for EACH GWAS model. If I want to use all SNPs, that can be done, too. Just need to use the full GWAS result files, or the "ManhattanPlot" files I saved from Step 1.
        start.time <- Sys.time()
        writeLines(capture.output(sessionInfo()), "sessionInfo_manhattan.txt") #save session info and check for compatibility later, if needed
        writeLines(capture.output(getwd()), "wd_manhattan.txt")#wd
        {#*---- a.  Load data ----
          # Read file that contains values for manhattan plots; Use ONLY significant SNPs
          EditTextHere <- paste("SigMarkerWiseSNP.list")
          dir(pattern=EditTextHere) #Confirm selection
          file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
        }
        {#*---- b. Making multiple manhattan plots ----
          lapply(file.list,function(x) {
            file <- fread(x) #read the data into R
            x1 <- strsplit(x,"\\.")
            id = x1[[1]][2]
            confirm =  x1[[1]][3]
            population = x1[[1]][1]
            file <-MVP.Report(file[,4:10],
                              plot.type=c("m"),
                              multracks=T, #If this one is true, will get 2 plots:multitracks and multitraits, instead of n plot (for n traits)
                              threshold=c(0.1/700000,1e-3),
                              threshold.lty=c(2,1),
                              threshold.lwd=c(2,2),
                              threshold.col=c("black","black"),
                              amplify=F,
                              ylim=c(3,round(-log10(min(file[,7:10], na.rm=F)))+3),
                              band=3,
                              file.type="pdf",
                              memo=paste0(population,".",id,".",confirm),
                              dpi=5,
                              verbose = T)
            return (file)
          })
        }
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        time.taken
        write(time.taken,"timetaken_manhattan.txt")
      }
      {#|-------- Step 4 Plotting (multiple) QQ plots with all GWAS p-values--------| ####
        #This step can take a long time; only run when needed
        start.time <- Sys.time()
        writeLines(capture.output(sessionInfo()), "sessionInfo_QQ.txt") #save session info and check for compatibility later, if needed
        writeLines(capture.output(getwd()), "wd_QQ.txt")#wd
        {#*---- a.  Load data ----
          # Read file that contains values for manhattan plots
          EditTextHere <- paste(".ManhattanPlot")
          dir(pattern=EditTextHere) #Confirm selection
          file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
        }
        {#*---- b. Making multiple QQ plots from full GWAS results ----
          lapply(file.list,function(x) {
            file <- fread(x) #read the data into R
            file <- rename_with(file, ~ gsub("_whole_", "_", .x, fixed = T))
            file <- rename_with(file, ~ gsub("_root_", "_", .x, fixed = T))
            file <- rename_with(file, ~ gsub("_shoot_", "_", .x, fixed = T))
            file <- rename_with(file,~ gsub("GAPIT.", "", .x, fixed = TRUE))
            file <- rename_with(file,~ gsub(".GWAS.Results", "", .x, fixed = TRUE))
            file <- rename_with(file,~ gsub("Blink.", "", .x, fixed = TRUE))
            file <- rename_with(file,~ gsub("FarmCPU.", "", .x, fixed = TRUE))
            file <- rename_with(file,~ gsub("MLM.", "", .x, fixed = TRUE))
            file <- rename_with(file,~ gsub("GLM.", "", .x, fixed = TRUE))
            
            x1 <- strsplit(x,"\\.")
            GWASmodel = x1[[1]][2]
            population = x1[[1]][1]
            organ = x1[[1]][3]
            file <-MVP.Report(file[,1:7],
                              plot.type=c("q"),
                              multracks=T, #If this one is true, will get 2 plots: overlayed QQ plots and multiple QQ plots on the same plane
                              file.type="pdf",
                              memo=paste0(population,".",GWASmodel,".",organ),
                              dpi=5,
                              verbose = T)
            return (file)
          })
        }
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        time.taken
        write(time.taken,"timetaken_QQ.txt")
      }
    }
    {#* Var Jap, pc6- Specify work directory and population ####
      setwd(paste0(workstation,"/varjap"))
      pop = "varjap-pc6"
      {#|--------Step 1 Merge p-value across traits--------| ####
        start.time <- Sys.time()
        writeLines(capture.output(sessionInfo()), "sessionInfo_mergesnp.txt") #save session info and check for compatibility later, if needed
        writeLines(capture.output(getwd()), "wd_mergesnp.txt")#wd
        {#*1. BLINK root - Merge p-value across traits for Manhattan Plot ----
          {#*---- a. Specify organ to automate file selection ----
            method = "Blink"
            organ = "root"
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
        }
        {#*2. BLINK shoot - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "Blink"
            organ = "shoot"
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
        }
        {#*3. BLINK whole - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "Blink"
            organ = "whole"
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*4.  Farm root - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "FarmCPU"
            
            organ = "root"
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*5. Farm shoot - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "FarmCPU"
            
            
            organ = "shoot"
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*6.  Farm whole - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "FarmCPU"
            
            organ = "whole"
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*7. MLM root - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "MLM"
            
            organ = "root"
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*8.  MLM shoot - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            method = "MLM"
            
            
            organ = "shoot"
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*9. MLM whole - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            
            method = "MLM"
            
            organ = "whole"
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*10.  GLM root - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            
            method = "GLM"
            organ = "root"
            
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*11.  GLM shoot - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            
            method = "GLM"
            
            
            organ = "shoot"
            
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
          
          
        }
        {#*12.  GLM whole - Merge p-value across traits for Manhattan Plot ####
          {#*---- a. Specify organ to automate file selection ----
            
            method = "GLM"
            
            organ = "whole"
            
          }
          {#*---- b. Select files from directory ----
            #--https://stackoverflow.com/questions/30829470/apply-a-function-from-a-specific-r-package-to-all-files-in-folder
            
            # File overview
            dir(pattern= "GAPIT.Association.GWAS_Results.*.csv") #Check how many files, file pattern
            
            # Select files
            EditTextHere <- paste0(method,".","*","_",organ) #Edit file name here
            dir(pattern=EditTextHere) #Confirm selection
            file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
            
            # Alternatives include:
            #file.list <- dir(pattern = "txt$")
            #file.list = dir(pattern="*01.csv")
            #file.list <- list.files(pattern="*.csv")
          }
          {#*---- c. Keep only SNP + P.value columns ----
            dataforManhattan <- lapply(file.list,function(x) {
              file <- fread(x) #read the data into R
              file <- subset(file, select=c("SNP","Chromosome", "Position", "P.value")) #subset only SNP + p.value columns (first 4 columns)
              trtname = gsub(paste0("GAPIT.Association.GWAS_Results.",method,".","|.csv"), "", x) # Get the start of filename prefix; make it as short and specific as possible
              file<- setnames(file, "P.value", trtname) #change each p.value column header to match trait name
              return (file)
            })
            # Select only SNP + P.value column from each file
            #--https://stackoverflow.com/questions/50513935/rename-every-column-of-df-with-a-part-of-their-filename-r
            
            #If I want to add index file name, I would follow these two instead:
            #--https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45522323#45522323
            #--https://stackoverflow.com/questions/58679646/is-there-a-way-to-set-filename-as-column-names-automatically
            #--https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
          }
          {#*---- d. Merge P.value from multi-traits into one file by SNP ----
            mergeCols <- c("SNP", "Chromosome", "Position")
            ManhattanReady = Reduce(function(...) merge(..., by = mergeCols, all=T),dataforManhattan)
          }
          {#*---- e. Export CSV file for record ----
            write.csv(ManhattanReady, file=paste(pop,method,organ,"ManhattanPlot.csv",sep="."), row.names=FALSE)
          }
        }
        
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        time.taken
        write(time.taken,"timetaken_mergesnp.txt")
      }
      {#|--------Step 2 Filter significant SNPs--------| ####
        start.time <- Sys.time()
        writeLines(capture.output(sessionInfo()), "sessionInfo_filtersnp.txt") #save session info and check for compatibility later, if needed
        writeLines(capture.output(getwd()), "wd_filtersnp.txt")#wd
        {#*---- a. Prepare data files ----
          # Read file that contains values for manhattan plots
          EditTextHere <- "*.ManhattanPlot" #Edit file name here
          dir(pattern=EditTextHere) #Confirm selection
          file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
          
          # Read csv files for manhattan plots
          allpvalue <- lapply(file.list,function(x) {
            file <- fread(x) #read the data into R
            file <- rename_with(file, ~ gsub("_whole", "_", .x, fixed = T))
            file <- rename_with(file, ~ gsub("_root", "_", .x, fixed = T))
            file <- rename_with(file, ~ gsub("_shoot", "_", .x, fixed = T))
            x1 <- strsplit(x,"\\.")
            organ = x1[[1]][3]
            GWASmodel = x1[[1]][2]
            population = x1[[1]][1]
            file <- cbind (population, GWASmodel, organ, file) # add file name as a new column
            return (file)
          })
          
        }
        {#*---- a extra. For GAPIT results with different file name format; Rename column by removing prefix and sffix ----
          allpvalue<-allpvalue %>%
            map(~ .x %>%
                  rename_with(~ gsub("GAPIT.", "", .x, fixed = TRUE))%>%
                  rename_with(~ gsub(".GWAS.Results", "", .x, fixed = TRUE))%>%
                  rename_with(~ gsub("Blink.", "", .x, fixed = TRUE))%>%
                  rename_with(~ gsub("FarmCPU.", "", .x, fixed = TRUE))%>%
                  rename_with(~ gsub("MLM.", "", .x, fixed = TRUE))%>%
                  rename_with(~ gsub("GLM.", "", .x, fixed = TRUE))
            )
        }
        {#*---- b. Select significant SNPs across N treatments for each organ ----
          # Filter only SNPs with p-value below threshold (Bonferoni Correction = 0.1 for 700,000k SNP)
          SignificantMarkerWiseSNP<- allpvalue %>%
            map(~ .x %>%
                  filter_at(vars(matches('_')), any_vars(.< (0.001))) #p-value < 0.001
            )
          SignificantMarkerWise10e4SNP<- allpvalue %>%
            map(~ .x %>%
                  filter_at(vars(matches('_')), any_vars(.< (0.0001))) #p-value < 0.001
            )
          SignificantExpWiseSNP<- allpvalue %>%
            map(~ .x %>%
                  filter_at(vars(matches('_')), any_vars(.< (0.1/700000))) # Bonferoni correction alpha = 0.1/number of SNPs used.
            )
          
          # Merge all significant SNPs into one data table
          AllSigMarkerWiseSNP = rbindlist(SignificantMarkerWiseSNP,fill=F, use.names=T) #match column names by name, not fill with NAs if missing
          AllSigMarkerWise10e4SNP = rbindlist(SignificantMarkerWise10e4SNP,fill=F, use.names=T ) #match column names by name, not fill with NAs if missing
          AllSigExpWiseSNP = rbindlist(SignificantExpWiseSNP,fill=F, use.names=T ) #match column names by name, not fill with NAs if missing
          #https://www.rdocumentation.org/packages/data.table/versions/1.14.8/topics/rbindlist
          
          write.csv(AllSigMarkerWiseSNP, file=paste(pop,"AllModels","SigMarkerWiseSNP","list.csv",sep="."), row.names=FALSE)
          write.csv(AllSigMarkerWise10e4SNP, file=paste(pop,"AllModels","SigMarkerWise10e4SNP","list.csv",sep="."), row.names=FALSE)
          write.csv(AllSigExpWiseSNP, file=paste(pop,"AllModels","SigExpWiseSNP","list.csv",sep="."), row.names=FALSE)
          
          # Count number of significant SNPs
          AllSigMarkerWiseSNP%>%
            #mutate(newgroup=paste(organ,GWASmodel,sep= "."))%>%
            group_by(population,organ, GWASmodel)%>%
            tally()%>%
            write.csv(., file=paste(pop,"SigMarkerWiseSNP","count.csv",sep="."),quote=F,row.names=F)
          
          AllSigMarkerWise10e4SNP%>%
            #mutate(newgroup=paste(organ,GWASmodel,sep= "."))%>%
            group_by(population,organ, GWASmodel)%>%
            tally()%>%
            write.csv(., file=paste(pop,"SigMarkerWise10e4SNP","count.csv",sep="."),quote=F,row.names=F)
          
          AllSigExpWiseSNP%>%
            #mutate(newgroup=paste(organ,GWASmodel,sep= "."))%>%
            group_by(population,organ, GWASmodel)%>%
            tally()%>%
            write.csv(., file=paste(pop,"SigExpWiseSNP","count.csv",sep="."),quote=F,row.names=F)
          
          # Save a separate csv file for each GWAS model
          for(i in unique(AllSigMarkerWiseSNP$GWASmodel)){write.csv(AllSigMarkerWiseSNP[AllSigMarkerWiseSNP$GWASmodel==i,],paste(pop,i,"model","SigMarkerWiseSNP","list.csv",sep="."),quote=F,row.names=F)} #write 4 new files based on GWAS models (all organs)
          for(i in unique(AllSigMarkerWise10e4SNP$GWASmodel)){write.csv(AllSigMarkerWise10e4SNP[AllSigMarkerWise10e4SNP$GWASmodel==i,],paste(pop,i,"model","SigMarkerWise10e4SNP","list.csv",sep="."),quote=F,row.names=F)} #write 4 new files based on GWAS models (all organs)
          for(i in unique(AllSigExpWiseSNP$GWASmodel)){write.csv(AllSigExpWiseSNP[AllSigExpWiseSNP$GWASmodel==i,],paste(pop,i,"model","SigExpWiseSNP","list.csv",sep="."),quote=F,row.names=F)} #write 4 new files based on GWAS models (all organs)
          for(i in unique(AllSigMarkerWiseSNP$organ)){write.csv(AllSigMarkerWiseSNP[AllSigMarkerWiseSNP$organ==i,],paste(pop,i,"organ","SigMarkerWiseSNP","list.csv",sep="."),quote=F,row.names=F)} #write 3 new files based on organs (all models)
          
        }
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        time.taken
        write(time.taken,"timetaken_filterSNP.txt")
        # The SNP will be browsed against the genome browser
        #--http://rs-bt-mccouch4.biotech.cornell.edu:443/cgi-bin/hgGateway?clade=poaceae&org=O.+sativa&db=orySat2
        #--https://snp-seek.irri.org/
      }
      {#|--------Step 3 Plotting (multiple) Manhattan Plots *with only significant SNPs --------| ####
        #Currently I set it to plot manhattan plots using only sig SNPS with p-values <0.001 using the files I saved from Step 2 for EACH GWAS model. If I want to use all SNPs, that can be done, too. Just need to use the full GWAS result files, or the "ManhattanPlot" files I saved from Step 1.
        start.time <- Sys.time()
        writeLines(capture.output(sessionInfo()), "sessionInfo_manhattan.txt") #save session info and check for compatibility later, if needed
        writeLines(capture.output(getwd()), "wd_manhattan.txt")#wd
        {#*---- a.  Load data ----
          # Read file that contains values for manhattan plots; Use ONLY significant SNPs
          EditTextHere <- paste("SigMarkerWiseSNP.list")
          dir(pattern=EditTextHere) #Confirm selection
          file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
        }
        {#*---- b. Making multiple manhattan plots ----
          lapply(file.list,function(x) {
            file <- fread(x) #read the data into R
            x1 <- strsplit(x,"\\.")
            id = x1[[1]][2]
            confirm =  x1[[1]][3]
            population = x1[[1]][1]
            file <-MVP.Report(file[,4:10],
                              plot.type=c("m"),
                              multracks=T, #If this one is true, will get 2 plots:multitracks and multitraits, instead of n plot (for n traits)
                              threshold=c(0.1/700000,1e-3),
                              threshold.lty=c(2,1),
                              threshold.lwd=c(2,2),
                              threshold.col=c("black","black"),
                              amplify=F,
                              ylim=c(3,round(-log10(min(file[,7:10], na.rm=F)))+3),
                              band=3,
                              file.type="pdf",
                              memo=paste0(population,".",id,".",confirm),
                              dpi=5,
                              verbose = T)
            return (file)
          })
        }
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        time.taken
        write(time.taken,"timetaken_manhattan.txt")
      }
      {#|-------- Step 4 Plotting (multiple) QQ plots with all GWAS p-values--------| ####
        #This step can take a long time; only run when needed
        start.time <- Sys.time()
        writeLines(capture.output(sessionInfo()), "sessionInfo_QQ.txt") #save session info and check for compatibility later, if needed
        writeLines(capture.output(getwd()), "wd_QQ.txt")#wd
        {#*---- a.  Load data ----
          # Read file that contains values for manhattan plots
          EditTextHere <- paste(".ManhattanPlot")
          dir(pattern=EditTextHere) #Confirm selection
          file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
        }
        {#*---- b. Making multiple QQ plots from full GWAS results ----
          lapply(file.list,function(x) {
            file <- fread(x) #read the data into R
            file <- rename_with(file, ~ gsub("_whole_", "_", .x, fixed = T))
            file <- rename_with(file, ~ gsub("_root_", "_", .x, fixed = T))
            file <- rename_with(file, ~ gsub("_shoot_", "_", .x, fixed = T))
            file <- rename_with(file,~ gsub("GAPIT.", "", .x, fixed = TRUE))
            file <- rename_with(file,~ gsub(".GWAS.Results", "", .x, fixed = TRUE))
            file <- rename_with(file,~ gsub("Blink.", "", .x, fixed = TRUE))
            file <- rename_with(file,~ gsub("FarmCPU.", "", .x, fixed = TRUE))
            file <- rename_with(file,~ gsub("MLM.", "", .x, fixed = TRUE))
            file <- rename_with(file,~ gsub("GLM.", "", .x, fixed = TRUE))
            
            x1 <- strsplit(x,"\\.")
            GWASmodel = x1[[1]][2]
            population = x1[[1]][1]
            organ = x1[[1]][3]
            file <-MVP.Report(file[,1:7],
                              plot.type=c("q"),
                              multracks=T, #If this one is true, will get 2 plots: overlayed QQ plots and multiple QQ plots on the same plane
                              file.type="pdf",
                              memo=paste0(population,".",GWASmodel,".",organ),
                              dpi=5,
                              verbose = T)
            return (file)
          })
        }
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        time.taken
        write(time.taken,"timetaken_QQ.txt")
      }
    }
}
#### Working Code - Compile P-values across ALL populations and models ####
#Evaluate Q-Q plots manually and select models with at least one significant markers at experiment-wise level for further analyses
#First, Put all SNP.list in the same folder (both experiment- and marker- wise list) before running the codes

setwd(paste0(workstation,"/Box/RiceGWAS2023/~")) #select folder where all list files are located

{#OK*---- Intermediate result "tableallSigMarkerWiseSNP" - Summarizing significant SNPs at 0.0001----
  EditTextHere <- "*.SigMarkerWise10e4SNP.list" #Edit file name here
  dir(pattern=EditTextHere) #Confirm selection
  file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files  
  
  # Read csv files
  allSigMarkerWiseSNP <- lapply(file.list,function(x) {
    file <- fread(x) #read the data into R
    return (file)
  })  
  {#*---- a extra. For GAPIT results with different file name format; Rename column by removing prefix and suffix ----
    allSigMarkerWiseSNP <-allSigMarkerWiseSNP  %>%
      map(~ .x %>% 
            rename_with(~ gsub("Chromosome", "Chr", .x, fixed = TRUE))%>%
            rename_with(~ gsub("Position", "Pos", .x, fixed = TRUE))
      ) %>%
      
      #Correct name from the new file name formatting 4/4/2023
      map(~ .x %>% mutate(GWASmodel = str_replace_all(GWASmodel, "Blink", "BLINK")))
    
  }
  tableallSigMarkerWiseSNP = rbindlist(allSigMarkerWiseSNP,fill=F, use.names=T ) #match column names by name, not fill with NAs if missing
  length(unique(tableallSigMarkerWiseSNP$SNP)) #9741 as of 5/14/2023
}
{#OK*---- Intermediate result "tableallSigExpWiseSNP" & TABLE S6 - Summarize number of shared conditions; SNPs that were significant at experiment-wise level across 336 conditions (7 pop x 4 models x 4 trt x 3 organs)----
  {#OK*----Read values ----
    # Read file that contain sig SNPs
    EditTextHere <- ".SigExpWiseSNP.list" #Edit file name here
    # EditTextHere <- "AllModels.SigExpWiseSNP.list" #Edit file name here
    dir(pattern=EditTextHere) #Confirm selection
    file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files  
    
    # Read csv files
    allSigExpWiseSNP <- lapply(file.list,function(x) {
      file <- fread(x) #read the data into R
      return (file)
    })  
    
    #Correct column names for some result files  
    {#*---- a extra. For GAPIT results with different file name format; Rename column by removing prefix and suffix ----
      allSigExpWiseSNP <-allSigExpWiseSNP  %>%
        map(~ .x %>% 
              rename_with(~ gsub("Chromosome", "Chr", .x, fixed = TRUE))%>%
              rename_with(~ gsub("Position", "Pos", .x, fixed = TRUE))
        )
    }
    
    #Bind all files together
    tableallSigExpWiseSNP = rbindlist(allSigExpWiseSNP,fill=F, use.names=T ) #216 SNP in total. match column names by names; This file has trt and exp rep combinations as variable name; OK but not the final file I should save. The actual number if I pivot longer = 236 SNP 
    
    #write.csv(tableallSigExpWiseSNP, file=paste("AllModels.SigExpWiseSNP.summary.list.csv",sep="."), row.names=FALSE)
    length(unique(tableallSigExpWiseSNP$SNP)) #176 for only selected models for rdp1-pc3, aus, ind, tej, trj, var ind pc0, var jap pc 6
    
    #The number of associations in each population here will be LESS than the actual one because this was selected on a row basis, based on trt x experiment replicate combinations. The actual numbers will be HIGHER when I put both trt and exp rep into individual variables.
  }
  {#OK*----Count number of significant associations  ---- 
    # The sum of each count will be >236 when arranged by individual condition; BUT, when pivot wider for individual factor (the number of rows will vary), the sum of count will be 236!
    {#OK*----Count number of significant association across ALL conditions (336 conditions across 4 factors) ---- 
      # The next step expand all variables into columns and then count the number of significant SNPs across conditions.
      
      # I must count cond for each factor separately and merge all the counts together by the SNP, Chr, Pos, pop, GWASmodel, organ, exp, trt later. I can't keep the count_xxx column because that would get different numbers and hence will be treated as distinct entries!!
      
      AllSigExpWiseSNPbyCond_all<- #Keep all p-values, but add count number of conditions with significant P (0.1/700000 markers)
        tableallSigExpWiseSNP %>%
        #Make trt and exp one of the variables from the column of trt x exp rep combinations
        pivot_longer( 
          cols = contains("_"),
          names_to = c("trt","exp"),
          names_sep = "_",
          values_to = "Pvalue",
          values_drop_na = F)%>%
        
        #Keep only entries with p-values
        filter(!is.na(Pvalue))%>%
        # filter(Pvalue<0.1/700000)%>% #get the same result regardless of this line because I counted based on this threshold anyway
        
        #Expand (widen) all the 3 factors (GWASmodel, organ, trt) and population into new factor combinations (AKA conditions), so that we can see how many conditions individual unique SNPs were significant across the whole experiment
        pivot_wider(
          id_cols = c(SNP,Chr,Pos),
          id_expand = FALSE,
          names_from = c(population, GWASmodel, organ, trt , exp),
          names_prefix = "FACTOR",
          names_sep = "_",
          names_glue = NULL,
          names_sort = FALSE,
          names_vary = "fastest",
          names_expand = FALSE,
          names_repair = "check_unique",
          values_from = Pvalue,
          values_fill = NULL,
          values_fn = NULL,
          unused_fn = NULL
        )%>%
        
        #Count the number of P-value across all factor combinations (336 columns) in each uniques SNP row
        mutate(count_all=rowSums(across(!c(SNP, Chr, Pos),~as.numeric(.x<(0.1/700000))),na.rm=TRUE))%>%
        arrange(desc(count_all))%>% #now we can use this "countcond" to filter; == "1" for unique condition, > 1 for shared across conds
        
        #Separate all factors into individual variable again    
        pivot_longer(
          cols = contains("FACTOR"),
          names_to = c("population","GWASmodel","organ","trt","exp"),
          names_prefix = "FACTOR",
          names_sep = "_",
          values_to = "Pvalue",
          values_drop_na = T #This is equivalent to filter(!is.na(P)) after pivot_longer; This values_drop_na leaves me with a tibble of 864 rows instead of 19,008 rows (with all NA entries if I choose values_drop_na = F)
          
        )%>%
        arrange(Pvalue)#this is the file for manhattan plot! (Contains 236 sig associations among all 864 p-values we selected from at least one n trt being significant )
    }
    {#OK*----Count number of significant association in each condition: Population ---- 
      
      AllSigExpWiseSNPbyCond_pop<- #Count for each factor - Population
        AllSigExpWiseSNPbyCond_all  %>%   
        pivot_wider(
          id_cols = c(SNP,Chr,Pos, GWASmodel, organ, trt , exp),
          id_expand = FALSE,
          names_from = c(population),
          names_prefix = "pop_",
          names_sep = "_",
          names_glue = NULL,
          names_sort = FALSE,
          names_vary = "fastest",
          names_expand = FALSE,
          names_repair = "check_unique",
          values_from = Pvalue,
          values_fill = NULL,
          values_fn = NULL,
          unused_fn = NULL
        )%>%
        
        mutate(count_pop=rowSums(across(starts_with("pop_"),~as.numeric(.x<(0.1/700000))),na.rm=TRUE))%>%
        
        #Separate all factors into individual variable again    
        pivot_longer(
          cols = contains("pop_"),
          names_to = c("population"),
          names_prefix = "pop_",
          #names_sep = "_",
          values_to = "Pvalue",
          values_drop_na = T 
        )%>%
        arrange(Pvalue)
    }
    {#OK*----Count number of significant association in each condition: GWAS Model ---- 
      
      AllSigExpWiseSNPbyCond_mod<- #Count for each factor - GWAS model
        AllSigExpWiseSNPbyCond_all  %>%   
        
        pivot_wider(
          id_cols = c(SNP,Chr,Pos,  population, organ, trt , exp),
          id_expand = FALSE,
          names_from = c(GWASmodel),
          names_prefix = "mod_",
          names_sep = "_",
          names_glue = NULL,
          names_sort = FALSE,
          names_vary = "fastest",
          names_expand = FALSE,
          names_repair = "check_unique",
          values_from = Pvalue,
          values_fill = NULL,
          values_fn = NULL,
          unused_fn = NULL
        )%>%
        
        
        mutate(count_mod=rowSums(across(starts_with("mod_"),~as.numeric(.x<(0.1/700000))),na.rm=TRUE))%>%
        
        #Separate all factors into individual variable again    
        pivot_longer(
          cols = contains("mod_"),
          names_to = c("GWASmodel"),
          names_prefix = "mod_",
          #names_sep = "_",
          values_to = "Pvalue",
          values_drop_na = T 
        )  %>%
        arrange(Pvalue)
    }
    {#OK*----Count number of significant association in each condition: Organ ---- 
      
      AllSigExpWiseSNPbyCond_org<- #Count for each factor - Organ
        AllSigExpWiseSNPbyCond_all  %>%   
        
        pivot_wider(
          id_cols = c(SNP,Chr,Pos, population, GWASmodel, trt , exp),
          id_expand = FALSE,
          names_from = c(organ),
          names_prefix = "organ_",
          names_sep = "_",
          names_glue = NULL,
          names_sort = FALSE,
          names_vary = "fastest",
          names_expand = FALSE,
          names_repair = "check_unique",
          values_from = Pvalue,
          values_fill = NULL,
          values_fn = NULL,
          unused_fn = NULL
        )%>%
        
        mutate(count_organ=rowSums(across(starts_with("organ_"),~as.numeric(.x<(0.1/700000))),na.rm=TRUE))%>%
        
        #Separate all factors into individual variable again    
        pivot_longer(
          cols = contains("organ_"),
          names_to = c("organ"),
          names_prefix = "organ_",
          #names_sep = "_",
          values_to = "Pvalue",
          values_drop_na = T
        )    %>%
        arrange(Pvalue)
    }
    {#OK*----Count number of significant association in each condition: N treatment ---- 
      
      AllSigExpWiseSNPbyCond_trt<- #Count for each factor - N treatment
        AllSigExpWiseSNPbyCond_all  %>%   
        
        pivot_wider(
          id_cols = c(SNP,Chr,Pos, population, GWASmodel, organ , exp),
          id_expand = FALSE,
          names_from = c(trt),
          names_prefix = "trt_",
          names_sep = "_",
          names_glue = NULL,
          names_sort = FALSE,
          names_vary = "fastest",
          names_expand = FALSE,
          names_repair = "check_unique",
          values_from = Pvalue,
          values_fill = NULL,
          values_fn = NULL,
          unused_fn = NULL
        )%>%
        
        mutate(count_trt=rowSums(across(starts_with("trt_"),~as.numeric(.x<(0.1/700000))),na.rm=TRUE))%>%
        
        #Separate all factors into individual variable again    
        pivot_longer(
          cols = contains("trt_"),
          names_to = c("trt"),
          names_prefix = "trt_",
          #names_sep = "_",
          values_to = "Pvalue",
          values_drop_na = T
        )      %>%
        arrange(Pvalue)
    }
   
  }
  {#OK*----TABLE S6 - Merge all counts together into one df and write csv ----    
    
    AllSigExpWiseSNPbyCond<- #864 intotal - but this includes p<0.1/700000
      AllSigExpWiseSNPbyCond_all%>%
      inner_join(AllSigExpWiseSNPbyCond_pop)%>%
      inner_join(AllSigExpWiseSNPbyCond_mod)%>%
      inner_join(AllSigExpWiseSNPbyCond_org)%>%
      inner_join(AllSigExpWiseSNPbyCond_trt)

    # Write CSV for only associations significant P-value 
    AllSigExpWiseSNPbyCond%>% 
      filter(Pvalue<0.1/700000)%>% #filter only 236 first row of AllSigExpWiseSNPbyCond
      arrange(Pvalue)%>% 
      dplyr::select (! exp)%>%
      write.csv(., file=paste("All","SigExWiseSNP","byCond","count.csv",sep="."),quote=F,row.names=F) #this is the list of all 931 exp wise significant association (updated using rdp1-pc3)
    length(unique(AllSigExpWiseSNPbyCond$SNP)) #176 for only selected models for rdp1-pc3, aus, ind, tej, trj, var ind pc0, var jap pc 6
    
  }
  #SUMMARY 5/14/23 - 236 Associations, 176 unique SNPs.
  #Still when looking at multiple conditions, I will need to look at count_all only. If only look at count_factor (any), I will only be looking at intersect, not union!
}
#### Working Code - Designate peak SNP with at least 3 adjacent significant SNP as QTL  ####
#Narrow down the candidate loci by filter only SNP with at least 3 co-located SNP being significant in the same population and nitrogen conditions

{#OK*---- Intermediate result "[pop]_[trt]_[Chr]_[SNP]_adjacentSNP_list" - Summarizing significant SNPs at 0.0001 adjacent to the highly significant SNP (exp-wise threshold)----
  #If run previous steps prior to this one, the dataframes containing ExpWise and MarkerWise SNP are in the global environment already
  #Specify directory (new folder) to save intermediate results
  setwd(paste0(workstation,"/USETHIS_MAF_QN/QTLbyOrgan/adjacentSNP_popxtrt0.0001")) 
  
  {#*----Step 1 - Match SNP at 0.001 within 150 kb of each exp wise SNP----
    require(data.table)
    #https://stackoverflow.com/questions/24480031/overlap-join-with-start-and-end-positions/25655497#25655497
    #https://stackoverflow.com/questions/27619381/data-frame-lookup-value-in-range-and-return-different-column
    #https://stackoverflow.com/questions/19748535/finding-overlapping-ranges-between-two-interval-data/19752720#19752720
    
    keySNP<-AllSigExpWiseSNPbyCond%>%
      filter(Pvalue<0.1/700000)%>% #filter only significant SNPs
      # filter(count_all>1)%>% #use this line if want to filter number of shared condition
      group_by(population, trt, Chr,SNP,Pos)%>%
      tally()%>%
      mutate(start=Pos-75000)%>%
      mutate(end=Pos+75000)%>%
      # arrange(desc(n))%>%
      as.data.frame()%>%
      setDT()%>%
      setkey(SNP)%>%
      select(c(population, trt, Chr,SNP,start,end))
    
    {#OK*----Make new markerwise SNP file to include trt and pop---- 
      AllSigMarkerWiseSNP_trtxpop_Long <- tableallSigMarkerWiseSNP%>%
        #Make trt and exp one of the variables from the column of trt x exp rep combinations
        pivot_longer( 
          cols = contains("_"),
          names_to = c("trt","exp"),
          names_sep = "_",
          values_to = "P",
          values_drop_na = F)%>%
        
        #Keep only entries with p-values
        filter(!is.na(P)) %>%
        
        filter(P<0.0001) %>%
        
        arrange(desc(P))%>% #The lowest P must be >= 0.0001
        
        group_by(population, trt, Chr,SNP,Pos)%>% 
        tally()%>%
        arrange(desc(n))#nrow = number of unique SNP
      }
    selectSNP<-AllSigMarkerWiseSNP_trtxpop_Long%>%
      select(c(population, trt, Chr,SNP,Pos))%>%
      setDT()%>%
      setkey(SNP)
    
    #This one works! 
    for(l in unique(keySNP$population)){
      for(k in unique(keySNP$trt)){
        for(j in unique(keySNP$Chr)){
          for(i in unique(keySNP$SNP)){
            write.csv(
              selectSNP[selectSNP$population==l&selectSNP$trt==k&selectSNP$Chr==j&inrange(Pos, keySNP[keySNP$population==l&keySNP$trt==k&keySNP$Chr==j&keySNP$SNP==i,]$start, keySNP[keySNP$population==l&keySNP$trt==k&keySNP$Chr==j&keySNP$SNP==i,]$end)],
              paste(l,k,j,i,"adjacentSNP","list.csv",sep="_"),quote=F,row.names=F)}
        }
      }
    }
  }
  {#OK*----Step 2 - Merge individual file of each peak SNP together ----
      {#OK*---- Select files from directory at 0.0001 ----
        # Select files; should get 197904 files!!!!
        EditTextHere <- paste0("*adjacentSNP_list") #Edit file name here
        dir(pattern=EditTextHere) #Confirm selection
        file.list = list.files(pattern=EditTextHere) #I choose this one to select all files I need; edit text before * to filter files
      }
      {#OK*---- Read file at 0.0001----
        alladjacentSNP_popxtrt <- lapply(file.list,function(x) {
          file <- fread(x) #read the data into R
          
          x1 <- strsplit(x,"\\_")
          # population = x1[[1]][1]
          # trt = x1[[1]][2]
          peakSNP = x1[[1]][4]
          
          file <- cbind (peakSNP, file) # add file name as a new column
          file<- setnames(file, "SNP", "adjacentSNP") 
          return (file)
        })
        
        tablealladjacentSNP_popxtrt = rbindlist(alladjacentSNP_popxtrt,fill=F, use.names=T ) #match column names by names
        length(unique(tablealladjacentSNP_popxtrt$peakSNP))
        length(unique(tablealladjacentSNP_popxtrt$Chr)) 
        tablealladjacentSNP_popxtrt%>%
          write.csv(., file=paste("Table","adjacentSNP_popxtrt","csv",sep="."),quote=F,row.names=F)
      }
    
  }
  {#OK*----Step 3 -  Count number of SNP----
    
    BySNP_tablealladjacentSNP_popxtrt <-
      tablealladjacentSNP_popxtrt%>% #keep only non-duplicate columns
      filter(!is.na(adjacentSNP))%>%
      group_by(Chr,peakSNP)%>%
      data.table()
    
    BySNP_tablealladjacentSNP_popxtrt_count <-
      tablealladjacentSNP_popxtrt%>% #keep only non-duplicate columns
      filter(!is.na(adjacentSNP))%>%
      group_by(Chr,peakSNP)%>%
      tally()%>%
      arrange(desc(n))%>%
      rename(n0.0001= "n")%>%
      data.table()
    
    BySNP_tablealladjacentSNP_popxtrt%>%
      write.csv(., file=paste("All","adjacentSNP_popxtrt","NoNA","csv",sep="."),quote=F,row.names=F) #this is the list of all unique peak SNP at exp threshold
    
    BySNP_tablealladjacentSNP_popxtrt_count%>%
      write.csv(., file=paste("All","adjacentSNP_popxtrt","count.csv",sep="."),quote=F,row.names=F) #this is the count of adjacent SNPs of all  unique peak SNP at exp threshold 
  }
}
{#OK*---- TABLE S7 - Select peak SNP (QTL) and prepare interval to query candidate genes within LD----
  {#OK*----Find range of each peak SNP----
  setwd(paste0(workstation,"/USETHIS_MAF_QN/ByOrgan")) 
  
  # For x unique SNPs, there are x SNPs within +-75kb of these peak SNPs
  # BySNP_tablealladjacentSNP_popxtrt<- #list all QTL with adjacent SNPs
    
    #load a file, if didn't run on the same workstation
    BySNP_tablealladjacentSNP_popxtrt<-read.csv("All.adjacentSNP_popxtrt.0.0001.NoNA.csv", header=T) # I went fo 10e-4
  
 min_pos<- BySNP_tablealladjacentSNP_popxtrt%>%
  group_by(peakSNP, population, trt)%>%
  slice(which.min(Pos))%>%
   rename(start = Pos)
 
 max_pos<- BySNP_tablealladjacentSNP_popxtrt%>%
   group_by(peakSNP, population, trt)%>%
   slice(which.max(Pos))%>%
   rename(end = Pos)
 
 BySNP_tablealladjacentSNP_popxtrt_range <- #This file can be used to find candidate genes within the range of each peak SNP; keep all the adjacent SNP for reference
 BySNP_tablealladjacentSNP_popxtrt%>%
   group_by(peakSNP, population, trt)%>%
   add_tally()%>%
   rename(n_adjacentSNP_popxtrt = n)%>%
   group_by(peakSNP)%>%
   add_tally()%>%
   rename(n_adjacentSNP = n)%>%
   rename(SNP=adjacentSNP)%>%
   
   left_join(min_pos)%>%
   rename(start_SNP = adjacentSNP)%>%
   left_join(max_pos)%>%
   rename(end_SNP = adjacentSNP)%>%
   group_by(Chr,peakSNP, population, trt)%>%
   arrange(Chr, peakSNP,population,trt)
  }
  {#OK*----Pre-DATAFILE Make a list with only 196 rows (keep only unique trt x pop combination and find candidate genes using this list----
    ByPeakSNP_tablealladjacentSNP_popxtrt_range<-#this one only filter out unique peakSNP for all population x trt combination (196 rather than 176 SNPs; some SNPs are significant in multiple pop x trt combinations)
      BySNP_tablealladjacentSNP_popxtrt_range%>% 
   distinct(population,trt, n_adjacentSNP,n_adjacentSNP_popxtrt,  Chr,peakSNP,start_SNP, end_SNP, start,end) %>%
      mutate(range_kb=(end-start)/1000)%>%
      mutate(search_term = paste0("chr",Chr, ":", start,"-",end))
    ByPeakSNP_tablealladjacentSNP_popxtrt_range%>%
      write.csv(., file=paste("PeakSNP","range","csv",sep="."),quote=F,row.names=F) 
      
  }
  {#OK*----Count number of peak SNP for each pop x trt (AKA QTL) BEFORE filtering out <3 adjacent SNP----   
    #OK Number of peak SNP (AKA QTL) for each population x trt combination (15 combinations; total = 196)
    BySNP_tablealladjacentSNP_popxtrt_range%>%
      dplyr::select(population,trt,peakSNP)%>%
      group_by(population,trt,peakSNP)%>%
      tally()%>%
      group_by(population,trt)%>%
      tally()%>%
      write.csv(., file=paste("PeakSNP","Bypopxtrt","count","csv",sep="."),quote=F,row.names=F)
    
    #OK Number of All SNP at 0.0001 for each population x trt combination (total = 2109 adjacent SNPs)
  BySNP_tablealladjacentSNP_popxtrt_range%>%
    dplyr::select(population,trt,peakSNP)%>%
    group_by(population,trt,peakSNP)%>% 
    tally()%>%
    write.csv(., file=paste("AdjacentSNP","Bypopxtrt","count","csv",sep="."),quote=F,row.names=F)  # This is similar to column "n_adjacentSNP" in "BySNP_tablealladjacentSNP_popxtrt_range"
  
  }
 {#OK*----DATAFILE TABLE S7 - Merge QTL with number of conditions; use n adjacentSNP for pop x trt >3 as a cut point and use the number of conditions to prioritize QTL---- 
 
   #Use this as main to join the range with, because within each pop x trt combination, each individual model x organ x exp rep does have varying number of count_xxx for each factor!
   
   QTLlistbypopxtrt_ALL<-
    AllSigExpWiseSNPbyCond%>% #864 SNP (including non-sig Pvalue)
  filter(Pvalue<0.1/700000)%>% #list of all 236 significant associations
   rename(peakSNP=SNP)%>%
     group_by(peakSNP,population,trt)%>% #196 unique SNP, from 236 sig associations
     add_tally()%>%
     rename(n_popxtrt=n) %>%#at this step, there are 236 rows of associations, grouped into 196 unique combinations of pop x trt
inner_join(ByPeakSNP_tablealladjacentSNP_popxtrt_range)%>% #We can use this list to search for the next step!
     filter(n_adjacentSNP_popxtrt>2)%>% #or filter to have at least 3 significant SNPs (0.001) at the peak SNP. Total of 138 peakSNP/QTL left (112 combination of peakSNP x pop x trt) from 196 SNPs previously. (i.e.58 SNP has <3 adjacent SNPs)
     arrange(desc(count_all))%>%#then sort by count_all
     
     #update trt name back to the full trt name
     mutate(trt = (replace(trt,trt == "H", "NH4+~10mM")))%>%
     mutate(trt = (replace(trt,trt == "L", "NH4+~0.3mM")))%>%
     mutate(trt = (replace(trt,trt == "M", "NH4+~3mM")))%>%
     mutate(trt = (replace(trt,trt == "N", "NO3-~3mM")))%>%
     group_by(peakSNP,population,trt) #group by pop x trt (138 QTL left)
   
      QTLlistbypopxtrt<-QTLlistbypopxtrt_ALL%>%
     slice(which.min(Pvalue))%>% #keep only the lowest p-value when the same QTL have multiple significant P-values within the same peakSNPxtrtxpop combination. Total of 112 QTL left (same as the number of peakSNP x pop x trt combinations)
     rename(Lowest_Pvalue=Pvalue)%>%
     dplyr::select (!c(exp))
   
   QTLlistbypopxtrt%>% #This file should guide plotting up manhattan plots!
     write.csv(., file=paste("RDP1_Nform_Table_S7_QTL_list","csv",sep="."),quote=F,row.names=F) #112 QTL
   
   AllSigExpWiseSNPbyCond%>% 
     filter(Pvalue<0.1/700000)%>% #list of all significant associations
     rename(peakSNP=SNP)%>%
     group_by(peakSNP,population,trt)%>%
     add_tally()%>%
     rename(n_popxtrt=n) %>%
     inner_join(ByPeakSNP_tablealladjacentSNP_popxtrt_range)%>% #We can use this list to search for the next step!
     filter(n_adjacentSNP_popxtrt>2)%>% #or filter to have at least 3 significant SNPs (0.0001) at the peak SNP
     arrange(desc(count_all))%>%#then sort by count_all
     
     #update trt name back to the full trt name
     mutate(trt = (replace(trt,trt == "H", "NH4+~10mM")))%>%
     mutate(trt = (replace(trt,trt == "L", "NH4+~0.3mM")))%>%
     mutate(trt = (replace(trt,trt == "M", "NH4+~3mM")))%>%
     mutate(trt = (replace(trt,trt == "N", "NO3-~3mM")))%>%
     group_by(peakSNP,population,trt)%>% #group by pop x trt
     slice(which.min(Pvalue))%>% #keep only the lowest p-value
     rename(Lowest_Pvalue=Pvalue)%>%
     dplyr::select (!c(exp))%>%
     group_by(peakSNP,population)%>% #We want to query candidate genes based on population and QTL range, so we can ignore from which trt they were significant. Next, I kept only the entry with the largest range to look for candidate genes. 
     slice(which.max(range_kb))%>% #keep only the largest range of each unique SNP in each population. So, we get 98 searchable ranges for individual QTL from 112 QTL. This means than among 98 QTL, some has multiple N conditions (e.g. 14 repeats; 112-98)
     
     write.csv(., file=paste("RDP1_Nform_Table_S7_unique_QTL_list","csv",sep="."),quote=F,row.names=F) #98 QTL
 }
}
+++++++++++++++++++++++++++++++++
#### Working Code - Data Reporting  ####
{#OK*---- FIGURE 5 & FIGURE S3 - Summarizing number of significant SNPs ----
  #https://upset.app/about/#about #https://github.com/hms-dbmi/UpSetR
  #https://cran.r-project.org/web/packages/pivottabler/vignettes/v00-vignettes.html
  {#OK-----Fig 5C UpsetR for shared SNP across N conditions-----
    #Recommended readings
    #https://cran.r-project.org/web/packages/ComplexUpset/index.html 
    #https://github.com/cxli233/customized_upset_plots
    #https://krassowski.github.io/complex-upset/articles/Examples_R.html#further-adjustments-using-ggplot2-functions
    #https://krassowski.github.io/complex-upset/index.html
    
    UpsetR_trt<-tableallSigExpWiseSNP %>%
      #Make trt and exp one of the variables from the column of trt x exp rep combinations
      pivot_longer( 
        cols = contains("_"),
        names_to = c("trt","exp"),
        names_sep = "_",
        values_to = "P",
        values_drop_na = F)%>%
      
      #Keep only entries with p-values
      filter(P<0.1/700000)%>%
      filter(!is.na(P))%>%
      
      #Count for each factor - Treatment
      pivot_wider(
        id_cols = c(SNP,Chr,Pos, population, GWASmodel, organ , exp),
        id_expand = FALSE,
        names_from = c(trt),
        names_prefix = "trt_",
        names_sep = "_",
        names_glue = NULL,
        names_sort = FALSE,
        names_vary = "fastest",
        names_expand = FALSE,
        names_repair = "check_unique",
        values_from = P,
        values_fill = 0,
        values_fn = NULL,
        unused_fn = NULL
      )%>%
      # mutate(across(starts_with("trt_"), replace(across(starts_with("trt_")),.x>0, "1")))
      mutate(trt_H = as.numeric(replace(trt_H,trt_H>0, "1")))%>%
      mutate(trt_L = as.numeric(replace(trt_L,trt_L>0, "1")))%>%
      mutate(trt_M = as.numeric(replace(trt_M,trt_M>0, "1")))%>%
      mutate(trt_N = as.numeric(replace(trt_N,trt_N>0, "1")))%>%
      # select(starts_with("trt_"))%>%
      rename("NH4+ 10 mM"=trt_H)%>%
      rename("NH4+ 0.3 mM"=trt_L)%>%
      rename("NH4+ 3 mM"=trt_M)%>%
      rename("NO3- 3 mM"=trt_N)%>%
      
      mutate(population = str_replace_all(population, "rdp1", "Whole Panel"))%>%
      mutate(population = str_replace_all(population, "varind-pc0", "INDICA"))%>%
      mutate(population = str_replace_all(population, "varjap-pc6", "JAPONICA"))%>%
      mutate(population = str_replace_all(population, "tej", "temperate-japonica"))%>%
      mutate(population = str_replace_all(population, "trj", "tropical-japonica"))%>%
      mutate(population = str_replace_all(population, "ind", "indica"))%>%
      mutate(population=factor(population, levels=c("Whole Panel","INDICA","aus","indica","JAPONICA","temperate-japonica", "tropical-japonica")))%>%
      as.data.frame()
    
    # UpsetR_trt%>%
    #   write.csv(., file=paste("UpsetR_trt","csv",sep="."),quote=F,row.names=F) 
    
    
    Fig_Upset_trt<-ComplexUpset::upset(
      UpsetR_trt, c("NH4+ 10 mM", "NH4+ 0.3 mM", "NH4+ 3 mM","NO3- 3 mM"),
      width_ratio=0.2,
      height_ratio=0.2,
      name="",
      min_size=1,
      # group_by='sets'
      # stripes=c("#20A387FF","#FDE725FF", "#39568CFF","#440154FF"),#https://www.thinkingondata.com/something-about-viridis-library/
      stripes='white',

      base_annotations=list(
        'Intersection size'=intersection_size(
          counts=TRUE,
          # mapping=aes(fill=population),
          text=list(
            vjust=0.5,
            hjust=0,
            angle=90,
            size=3
          )
        )+ ylab('Number of\n association')
        +theme(panel.background = element_rect(fill = 'transparent', color = "black",linewidth = .5, linetype = 1))
        +theme(panel.border = element_rect(fill = 'transparent', color = "black",linewidth = .5, linetype = 1))
        +theme(panel.grid.major.y =element_blank(),panel.grid.major.x =element_blank(),panel.grid.minor.y =element_blank())
        +theme(plot.tag = element_text(face = 'bold', size = 10))
        + theme(axis.text.y=element_text(angle=0, size=8))
        + theme(axis.title.y=element_text(angle=90, size=8))
      ),
      set_sizes=FALSE,
      matrix=(
        intersection_matrix(
          geom=geom_point(
            shape='square filled',
            # color="white",
            # fill="red",
            stroke = 0.1,
            size=2.5)
        )
        +  scale_color_manual(
          values=c('NH4+ 10 mM'="#440154FF", 'NH4+ 0.3 mM'="#39568CFF", 'NH4+ 3 mM'="#20A387FF", 'NO3- 3 mM'="#FDE725FF"))),
      
      queries=list(
        upset_query(set='NH4+ 10 mM', fill="#440154FF"),
        upset_query(set='NH4+ 0.3 mM', fill="#39568CFF"),
        upset_query(set='NH4+ 3 mM', fill="#20A387FF"),
        upset_query(set='NO3- 3 mM', fill="#FDE725FF"))
    )
    Fig_Upset_trt
  }
  {#OK-----Fig 5AB + Fig supp Number of association by each component-----
    
    AllSigExpWiseSNPby_Pop<- 
      AllSigExpWiseSNPbyCond%>%       
      filter(Pvalue<0.1/700000)%>% 
      mutate(population = str_replace_all(population, "rdp1-pc3", "Whole Panel"))%>% 
      mutate(population = str_replace_all(population, "varind-pc0", "INDICA"))%>% 
      mutate(population = str_replace_all(population, "varjap-pc6", "JAPONICA"))%>% 
      mutate(population = str_replace_all(population, "tej", "temperate\njaponica"))%>% 
      mutate(population = str_replace_all(population, "trj", "tropical\njaponica"))%>% 
      mutate(population = str_replace_all(population, "ind", "indica"))%>% 
      
      group_by(population)%>%
      tally()%>%
      
      ggplot(aes(x=factor(population, levels=c("Whole Panel","INDICA","aus","indica","JAPONICA","temperate\njaponica", "tropical\njaponica")), 
                 y=n))+
      
      geom_bar(stat = "identity")+
      
      xlab("")+
      ylab("Number of \nassociation")+
      theme(panel.background = element_rect(fill = "white", color = "black",linewidth = .5, linetype = 1),
            panel.border = element_rect(fill = "transparent", color = "black",linewidth = 0.75, linetype = 1),
            # panel.grid.major.y =element_line(color = "gray",
            #                                  linewidth = 0.25,
            #                                  linetype = 2),
            # 
            legend.text = element_text(size=6),
            legend.title = element_text(
              size=7,
              face = "bold"
            ),
            # legend.key = element_rect(fill = "white", color = "black"),
            
            legend.position="none",
            # legend.box.background = element_rect(fill = "white", colour = "black"),
            # legend.position = c(.93, .83),
            # legend.justification = c("right", "top"),
            # legend.box.just = "right",
            # legend.margin = margin(3, 3, 3, 3),
            # axis.ticks.x =  element_blank(),
            axis.title.y=element_text(size=7),   
            axis.title.x=element_blank(),
            axis.text.x = element_text(angle = 30, vjust = 0, hjust=0,size = 5),
            axis.text.y = element_text(size = 8),
            axis.text = element_text(color = "black"),
            
            # strip.text.y = element_text(size=9,face="bold",angle = 270),
            # strip.text.x = element_blank(),
            # strip.text.x = element_text(size=9,face="bold"),
            
            # strip.background = element_rect(colour="NA", fill="white")
            # plot.tag = element_text(face = 'bold', size = 10)
            plot.tag = element_text(size = 12)
            
            
      )+
      # scale_fill_viridis(discrete = TRUE)+
      scale_x_discrete(position = "top",
                       breaks=c("Whole Panel","INDICA","aus","indica","JAPONICA","temperate\njaponica", "tropical\njaponica"),
                       labels=c("Whole Panel","INDICA","aus","indica","JAPONICA","temperate\njaponica", "tropical\njaponica")) +
      
      geom_text(aes(label=n,vjust =-0.2, hjust=0.5),size=2, color="black")+
      scale_y_continuous(breaks = seq(0, 150, by=50),limits=c(0,160))
    
    AllSigExpWiseSNPby_Pop
    
    AllSigExpWiseSNPby_PopxModel<- 
      AllSigExpWiseSNPbyCond%>%       
      filter(Pvalue<0.1/700000)%>% 
      mutate(population = str_replace_all(population, "rdp1-pc3", "Whole Panel"))%>% 
      mutate(population = str_replace_all(population, "varind-pc0", "INDICA"))%>% 
      mutate(population = str_replace_all(population, "varjap-pc6", "JAPONICA"))%>% 
      mutate(population = str_replace_all(population, "tej", "temperate-japonica"))%>% 
      mutate(population = str_replace_all(population, "trj", "tropical-japonica"))%>% 
      mutate(population = str_replace_all(population, "ind", "indica"))%>% 
      
      mutate(GWASmodel = str_replace_all(GWASmodel, "Blink", "BLINK"))%>% 
      
      
      group_by(population,GWASmodel)%>%
      tally()%>%
      
      ggplot(aes(x=factor(population, levels=c("Whole Panel","INDICA","aus","indica","JAPONICA","temperate-japonica", "tropical-japonica")), 
                 y=n,
                 fill=factor(GWASmodel, levels=c("BLINK","FarmCPU","MLM","GLM")))
      )+
      
      geom_bar(stat = "identity",position = "fill")+
      
      xlab("")+
      ylab("Proportion")+
      theme(panel.background = element_rect(fill = "white", color = "black",linewidth = .5, linetype = 1),
            panel.border = element_rect(fill = "transparent", color = "black",linewidth = 0.75, linetype = 1),
            # panel.grid.major.y =element_line(color = "gray",
            #                                  linewidth = 0.25,
            #                                  linetype = 2),
            
            legend.text = element_text(size=6),
            legend.title = element_text(
              size=7,
              face = "bold"
            ),
            # legend.key = element_rect(fill = "white", color = "black"),
            
            legend.position="right",
            # legend.box.background = element_rect(fill = "white", colour = "black"),
            # legend.position = c(.93, .83),
            # legend.justification = c("right", "top"),
            # legend.box.just = "right",
            # legend.margin = margin(3, 3, 3, 3),
            axis.ticks.x =  element_blank(),
            axis.title.y=element_text(size=7),   
            axis.title.x=element_blank(),
            axis.text.x = element_blank(),
            # axis.text.x = element_text(size=7),
            axis.text.y = element_text(size = 8),
            axis.text = element_text(color = "black"),
            
            # strip.text.y = element_text(size=9,face="bold",angle = 270),
            # strip.text.x = element_blank(),
            # strip.text.x = element_text(size=9,face="bold"),
            
            # strip.background = element_rect(colour="NA", fill="white")
            plot.tag = element_text( size = 12)
            
      )+
      scale_fill_viridis(discrete = TRUE,name="GWAS Model",
                         # breaks=c("BLINK","FarmCPU","GLM", "MLM"),
                         # labels=c("BLINK","FarmCPU","GLM", "MLM"))
                         #New analyses did not include any GLM results due to deviation from normality in QQplots
                         breaks=c("BLINK","FarmCPU", "MLM", "GLM"),
                         labels=c("BLINK","FarmCPU", "MLM","GLM"))+
      scale_y_continuous(breaks = seq(0, 1, by=0.2))
    
    AllSigExpWiseSNPby_PopxModel
    
    
    
    
    AllSigExpWiseSNPby_PopxOrgan<- 
      AllSigExpWiseSNPbyCond%>%       
      filter(Pvalue<0.1/700000)%>% 
      mutate(population = str_replace_all(population, "rdp1-pc3", "Whole Panel"))%>% 
      mutate(population = str_replace_all(population, "varind-pc0", "INDICA"))%>% 
      mutate(population = str_replace_all(population, "varjap-pc6", "JAPONICA"))%>% 
      mutate(population = str_replace_all(population, "tej", "temperate-japonica"))%>% 
      mutate(population = str_replace_all(population, "trj", "tropical-japonica"))%>% 
      mutate(population = str_replace_all(population, "ind", "indica"))%>% 
      
      group_by(population,organ)%>%
      tally()%>%
      ggplot(aes(x=factor(population, levels=c("Whole Panel","INDICA","aus","indica","JAPONICA","temperate-japonica", "tropical-japonica")), 
                 y=n,
                 fill=factor(organ, levels=c("whole","shoot","root")))
      )+
      
      geom_bar(stat = "identity",position = "fill")+
      
      xlab("")+
      ylab("Proportion")+
      theme(panel.background = element_rect(fill = "white", color = "black",linewidth = .5, linetype = 1),
            panel.border = element_rect(fill = "transparent", color = "black",linewidth = 0.75, linetype = 1),
            # panel.grid.major.y =element_line(color = "gray",
            #                                  linewidth = 0.25,
            #                                  linetype = 2),
            
            legend.text = element_text(size=6),
            legend.title = element_text(
              size=7,
              face = "bold"
            ),
            # legend.key = element_rect(fill = "white", color = "black"),
            
            legend.position="right",
            # legend.box.background = element_rect(fill = "white", colour = "black"),
            # legend.position = c(.93, .83),
            # legend.justification = c("right", "top"),
            # legend.box.just = "right",
            # legend.margin = margin(3, 3, 3, 3),
            axis.ticks.x =  element_blank(),
            axis.title.y=element_text(size=7),   
            axis.title.x=element_blank(),
            axis.text.x = element_blank(),
             # axis.text.x = element_text(size=7),
            axis.text.y = element_text(size = 8),
            axis.text = element_text(color = "black"),
            
            # strip.text.y = element_text(size=9,face="bold",angle = 270),
            # strip.text.x = element_blank(),
            # strip.text.x = element_text(size=9,face="bold"),
            
            # strip.background = element_rect(colour="NA", fill="white")
            plot.tag = element_text( size = 12)
            
      )+
      scale_fill_viridis(discrete = TRUE,name="Plant Organ",
                         breaks=c("whole","shoot","root"),
                         labels=c("Whole","Shoot","Root"))+
      scale_y_continuous(breaks = seq(0, 1, by=0.2))
    AllSigExpWiseSNPby_PopxOrgan
    
    
    AllSigExpWiseSNPby_PopxTrt<- 
      AllSigExpWiseSNPbyCond%>%       
      filter(Pvalue<0.1/700000)%>% 
      mutate(population = str_replace_all(population, "rdp1-pc3", "Whole Panel"))%>% 
      mutate(population = str_replace_all(population, "varind-pc0", "INDICA"))%>% 
      mutate(population = str_replace_all(population, "varjap-pc6", "JAPONICA"))%>% 
      mutate(population = str_replace_all(population, "tej", "temperate-japonica"))%>% 
      mutate(population = str_replace_all(population, "trj", "tropical-japonica"))%>% 
      mutate(population = str_replace_all(population, "ind", "indica"))%>% 
      group_by(population,trt)%>%
      tally()%>%
      ggplot(aes(x=factor(population, levels=c("Whole Panel","INDICA","aus","indica","JAPONICA","temperate-japonica", "tropical-japonica")), 
                 y=n,
                 fill=trt)
      )+
      
      geom_bar(stat = "identity",position = "fill")+
      
      xlab("")+
      ylab("Proportion")+
      theme(panel.background = element_rect(fill = "white", color = "black",linewidth = .5, linetype = 1),
            panel.border = element_rect(fill = "transparent", color = "black",linewidth = 0.75, linetype = 1),
            # panel.grid.major.y =element_line(color = "gray",
            #                                  linewidth = 0.25,
            #                                  linetype = 2),
            
            legend.text = element_text(size=6),
            legend.title = element_text(
              size=7,
              face = "bold"
            ),
            # legend.key = element_rect(fill = "white", color = "black"),
            
            # legend.position="right",
            legend.position="none",
            # legend.box.background = element_rect(fill = "white", colour = "black"),
            # legend.position = c(.93, .83),
            # legend.justification = c("right", "top"),
            # legend.box.just = "right",
            # legend.margin = margin(3, 3, 3, 3),
            axis.ticks.x =  element_blank(),
            axis.title.y=element_text(size=7),   
            axis.title.x=element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 8),
            axis.text = element_text(color = "black"),
            
            # strip.text.y = element_text(size=9,face="bold",angle = 270),
            # strip.text.x = element_blank(),
            # strip.text.x = element_text(size=9,face="bold"),
            
            # strip.background = element_rect(colour="NA", fill="white")
            # plot.tag = element_text(face = 'bold', size = 10)
            plot.tag = element_text(size = 12)            
      )+
      scale_fill_viridis(discrete = TRUE,name="Nitrogen Treatment",
                         breaks=c("H","L","M","N"),
                         labels=c("NH4+ 10 mM","NH4+ 0.3 mM","NH4+ 3 mM","NO-3 3 mM"))+
      scale_y_continuous(breaks = seq(0, 1, by=0.2))
    AllSigExpWiseSNPby_PopxTrt
    }
  {#OK*----FIG 5D QTL by pop x trt----   
    Fig5D<-AllSigExpWiseSNPbyCond%>% 
      filter(Pvalue<0.1/700000)%>% #list of all significant associations
      rename(peakSNP=SNP)%>%
      group_by(peakSNP,population,trt)%>%
      add_tally()%>%
      rename(n_popxtrt=n) %>%#at this step, there are 236 rows of associations, grouped into 196 unique combinations of pop x trt
      inner_join(ByPeakSNP_tablealladjacentSNP_popxtrt_range)%>% #We can use this list to search for the next step!
      filter(n_adjacentSNP_popxtrt>2)%>% #or filter to have at least 3 significant SNPs (0.0001) at the peak SNP
      arrange(desc(count_all))%>%#then sort by count_all
      
      
      #update trt name back to the full trt name
      mutate(trt = (replace(trt,trt == "H", "NH4+ 10 mM")))%>%
      mutate(trt = (replace(trt,trt == "L", "NH4+ 0.3 mM")))%>%
      mutate(trt = (replace(trt,trt == "M", "NH4+ 3 mM")))%>%
      mutate(trt = (replace(trt,trt == "N", "NO3- 3 mM")))%>%
      mutate(trt=factor(trt,levels=c("NH4+ 10 mM", "NH4+ 0.3 mM", "NH4+ 3 mM","NO3- 3 mM")))%>%
      
      slice(which.min(Pvalue)) %>%
      group_by(population,trt)%>% #group by pop x trt (138 rows of 112 unique peakSNP x pop x trt) into 9 remaining group
      #If I group by group_by(peakSNP, population,trt), I would get 138 rows -- if only filtered one from each QTL (lowest p-value only), I will get 112.
      
       tally()%>%   
      #update population name
      mutate(population = str_replace_all(population, "rdp1-pc3", "Whole Panel"))%>%
      mutate(population = str_replace_all(population, "varind-pc0", "INDICA"))%>%
      mutate(population = str_replace_all(population, "varjap-pc6", "JAPONICA"))%>%
      mutate(population = str_replace_all(population, "tej", "temperate\njaponica"))%>%
      mutate(population = str_replace_all(population, "trj", "tropical\njaponica"))%>%
      mutate(population = str_replace_all(population, "ind", "indica"))%>%
      mutate(population=factor(population, levels=c("Whole Panel","INDICA","aus","indica","JAPONICA","temperate\njaponica", "tropical\njaponica")))%>%
      
      ggplot(aes(x=trt, 
                 y=n, fill=trt)
             
      )+
      
      geom_bar(stat="identity",position="stack")+
      facet_grid(.~population)+ 
      
      xlab("Nitrogen Treatment")+
      ylab("Number of peak SNP (QTL)")+
      theme(panel.background = element_rect(fill = "white", color = "black",linewidth = 0.75, linetype = 1),
            panel.border = element_rect(fill = "transparent", color = "black",linewidth = 0.75, linetype = 1),
            panel.grid.major.y =element_blank(),
            panel.grid.major.x =element_blank(),
            panel.grid.minor.y =element_blank(),
            panel.grid.minor.x =element_blank(),
            
            legend.text = element_text(size=7),
            legend.title = element_text(
              size=7,
              face = "bold"
            ),
            # legend.key = element_rect(fill = "white", color = "black"),
            
            legend.position="right",
            #legend.box.background = element_rect(fill = "white", colour = "black"),
            # legend.position = c(.93, .83),
            # legend.justification = c("right", "top"),
            # legend.box.just = "right",
            # legend.margin = margin(3, 3, 3, 3),
            axis.ticks.x =  element_blank(),
            axis.title.y=element_text(size=7),   
            axis.title.x=element_text(size=7),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 8),
            axis.text = element_text(color = "black"),
            
            #strip.text.y = element_text(size=8,face="bold",angle = 270),
            # strip.text.x = element_blank(),
            strip.text.x = element_text(size=6),
            
            strip.background = element_rect(colour="NA", fill="white"),
            plot.tag = element_text(face = 'bold', size = 8)
      )+ 
      guides(fill=guide_legend(title=""))+
      scale_fill_viridis(discrete = TRUE)+
       geom_text(aes(label = n), position = position_dodge(width = 1), vjust = -0.5,size=2)+
      scale_y_continuous(breaks = seq(0, 60, by=20),limits=c(0,70))
    
    
    Fig5D
  }
  {#OK-----FIGURE 5 & FIGURE S3 - Save as figures-----
    #Variations I could tag
    # + plot_annotation("A")
    # + plot_annotation(tag_levels = "A")
    # +labs(tag = c("C","B"))+ plot_annotation(tag_levels = "A")
    
    AllSigExpWiseSNPby_Pop/AllSigExpWiseSNPby_PopxTrt +plot_layout(ncol = 1)+plot_annotation(tag_levels = list(c("A", "B")))
    ggsave("RDP1_Nform_Figure_5AB.png", dpi=300, height=3.5, width=3, units="in")
    ggsave("RDP1_Nform_Figure_5AB.pdf", dpi=300, height=3.5, width=3, units="in")     

    Fig_Upset_trt +plot_layout(ncol = 1) + plot_annotation("C")#I can't make it bold...
    ggsave("RDP1_Nform_Figure_5C.png", dpi=300, height=3.5, width=3, units="in", bg="transparent")
    ggsave("RDP1_Nform_Figure_5C.pdf", dpi=300, height=3.5, width=3, units="in", bg="transparent")
    
    Fig5D +plot_layout(ncol = 1)+ plot_annotation("D") #I can't make it bold...
    ggsave("RDP1_Nform_Figure_5D.png", dpi=300, height=2, width=6, units="in", bg="transparent")
    ggsave("RDP1_Nform_Figure_5D.pdf", dpi=300, height=2, width=6, units="in", bg="transparent")
    
    #Keep only popxtrt; The rest go to supplementary
    (AllSigExpWiseSNPby_Pop/AllSigExpWiseSNPby_PopxOrgan/AllSigExpWiseSNPby_PopxModel +plot_layout(ncol = 1)+ plot_annotation(tag_levels = "A"))
     ggsave("RDP1_Nform_Figure_S3.png", dpi=300, height=6, width=4, units="in")
     ggsave("RDP1_Nform_Figure_S3.pdf", dpi=300, height=6, width=4, units="in")
  }
}
{#OK*---- FIGURE 6 - Visualize QTL as manhattan plots by N treatments----
  #Data I can use to plot manhattan.
  #AllSigExpWiseSNPbyCond # this one has only SNP with at least one N condition > exp-wise threshold
  tableallSigMarkerWiseSNP # this one has only SNP with at least one N condition > marker-wise threshold
  QTLlistbypopxtrt_ALL # this one has only QTL with >= 3 adjacent SNP
  
  #Make two sets of plots, each with different N. 
  #6A - all marker-wise significant
  #6B - QTL only (138 points)
{#*----JPG----
{#OK*---- 6A ----
    MarkerWiseSig<-tableallSigMarkerWiseSNP%>%
      #Make trt and exp one of the variables from the column of trt x exp rep combinations
      pivot_longer( 
        cols = contains("_"),
        names_to = c("trt","exp"),
        names_sep = "_",
        values_to = "P",
        values_drop_na = F)%>% #409744
      
      #Keep only entries with p-values
      filter(!is.na(P)) %>% #409744
      
      #MUST FILTER ONLY P<0.0001 because now that I pivot long, there will be many entries below the threshold!
      filter(P<0.0001)%>%
      #filter(P<0.1/700000) #If I add this step I could just get the same as pivotlonger "tableallSigExpWiseSNP"
      mutate(trt = (replace(trt,trt == "H", "NH4+ 10 mM")))%>%
      mutate(trt = (replace(trt,trt == "L", "NH4+ 0.3 mM")))%>%
      mutate(trt = (replace(trt,trt == "M", "NH4+ 3 mM")))%>%
      mutate(trt = (replace(trt,trt == "N", "NO3- 3 mM")))%>%
      pivot_wider(
        id_cols = c(population,GWASmodel,organ,SNP,Chr,Pos,exp),
        id_expand = FALSE,
        names_from = trt,
        names_prefix = "",
        names_sep = "_",
        names_glue = NULL,
        names_sort = FALSE,
        names_vary = "fastest",
        names_expand = FALSE,
        names_repair = "check_unique",
        values_from = P,
        values_fill = NULL,
        values_fn = NULL,
        unused_fn = NULL
      )%>% #16836
      dplyr::select (!c(population,GWASmodel,organ,exp))

    MVP.Report(MarkerWiseSig,
               plot.type=c("m"),
               multracks=T, #If this one is true, will get 2 plots:multitracks and multitraits, instead of n plot (for n traits)
             
               threshold=c(0.1/700000,0.1/700000),
               threshold.lty=c(2,2),
               threshold.lwd=c(2,2),
               threshold.col=c("black","black"),
               
               col=c("black","gray"),
               
               signal.col = "blue",
               signal.pch = 19,
               signal.cex=c(0.5,1),
               signal.line=0,
               
               amplify=T,
               
               ylim=c(4,20),
               band=3,
               
               file.type="jpg",
               
               memo=paste0("6A"),
               dpi=300,
               verbose = T)
    
    
  
  
  }
{#OK*---- 6B ----
   QTLforManhattan<-
    QTLlistbypopxtrt_ALL%>%
    ungroup()%>%
    dplyr::select (c(peakSNP, Chr,Pos, trt, Pvalue))%>%
    dplyr::rename(SNP=peakSNP)%>%
    dplyr::rename(P=Pvalue)%>%
      #filter(P<0.1/700000) #If I add this step I could just get the same as pivotlonger "tableallSigExpWiseSNP"
      mutate(trt = (replace(trt,trt == "NH4+~10mM", "NH4+ 10 mM")))%>%
      mutate(trt = (replace(trt,trt == "NH4+~0.3mM", "NH4+ 0.3 mM")))%>%
      mutate(trt = (replace(trt,trt == "NH4+~3mM", "NH4+ 3 mM")))%>%
      mutate(trt = (replace(trt,trt == "NO3-~3mM", "NO3- 3 mM")))%>%
    mutate(row = row_number()) %>%
      pivot_wider(
        id_cols = c(row,SNP,Chr,Pos),
        id_expand = FALSE,
        names_from = trt,
        names_prefix = "",
        names_sep = "",
        names_glue = NULL,
        names_sort = FALSE,
        names_vary = "fastest",
        names_expand = FALSE,
        names_repair = "check_unique",
        values_from = P,
        values_fill = NULL
      )%>%
    dplyr::  select(-row) %>% 
    relocate('NO3- 3 mM',.before ='NH4+ 0.3 mM')%>% 
    relocate('NH4+ 0.3 mM',.before ='NH4+ 10 mM')%>% 
    relocate('NH4+ 10 mM',.before ='NH4+ 3 mM')

    MVP.Report(QTLforManhattan,
               plot.type=c("m"),
               multracks=T, #If this one is true, will get 2 plots:multitracks and multitraits, instead of n plot (for n traits)
             
               threshold=c(0.1/700000,0.1/700000),
               threshold.lty=c(2,2),
               threshold.lwd=c(2,2),
               threshold.col=c("black","black"),
               
               col=c("black","gray"),
               
               signal.col = "blue",
               signal.pch = 19,
               signal.cex=c(0.5,1),
               signal.line=0,
               
               amplify=T,
               
               highlight = c('SNP-4.23818648.','SNP-4.23792662.','SNP-12.25853862.','SNP-12.25890218.','SNP-12.25901203.'),
               highlight.col="red",
               highlight.pch = 17,
               highlight.cex =5, 
               
               ylim=c(4,20),
               band=3,
               
               file.type="jpg",
               
               memo=paste0("6B"),
               dpi=300,
               verbose = T)
  
  }
}
{#*----pdf----
  {#OK*---- 6A ----
    MarkerWiseSig<-tableallSigMarkerWiseSNP%>%
      #Make trt and exp one of the variables from the column of trt x exp rep combinations
      pivot_longer( 
        cols = contains("_"),
        names_to = c("trt","exp"),
        names_sep = "_",
        values_to = "P",
        values_drop_na = F)%>% #409744
      
      #Keep only entries with p-values
      filter(!is.na(P)) %>% #409744
      
      #MUST FILTER ONLY P<0.0001 because now that I pivot long, there will be many entries below the threshold!
      filter(P<0.0001)%>%
      #filter(P<0.1/700000) #If I add this step I could just get the same as pivotlonger "tableallSigExpWiseSNP"
      mutate(trt = (replace(trt,trt == "H", "NH4+ 10 mM")))%>%
      mutate(trt = (replace(trt,trt == "L", "NH4+ 0.3 mM")))%>%
      mutate(trt = (replace(trt,trt == "M", "NH4+ 3 mM")))%>%
      mutate(trt = (replace(trt,trt == "N", "NO3- 3 mM")))%>%
      pivot_wider(
        id_cols = c(population,GWASmodel,organ,SNP,Chr,Pos,exp),
        id_expand = FALSE,
        names_from = trt,
        names_prefix = "",
        names_sep = "_",
        names_glue = NULL,
        names_sort = FALSE,
        names_vary = "fastest",
        names_expand = FALSE,
        names_repair = "check_unique",
        values_from = P,
        values_fill = NULL,
        values_fn = NULL,
        unused_fn = NULL
      )%>% #16836
      dplyr::select (!c(population,GWASmodel,organ,exp))
    
    MVP.Report(MarkerWiseSig,
               plot.type=c("m"),
               multracks=T, #If this one is true, will get 2 plots:multitracks and multitraits, instead of n plot (for n traits)
               
               threshold=c(0.1/700000,0.1/700000),
               threshold.lty=c(2,2),
               threshold.lwd=c(2,2),
               threshold.col=c("black","black"),
               
               col=c("black","gray"),
               
               signal.col = "blue",
               signal.pch = 19,
               signal.cex=c(0.5,1),
               signal.line=0,
               
               amplify=T,
               
               ylim=c(4,20),
               band=3,
               
               file.type="pdf",
               
               memo=paste0("6A"),
               dpi=300,
               verbose = T)
  }
  {#OK*---- 6B ----
    QTLforManhattan<-
      QTLlistbypopxtrt_ALL%>%
      ungroup()%>%
      dplyr::select (c(peakSNP, Chr,Pos, trt, Pvalue))%>%
      dplyr::rename(SNP=peakSNP)%>%
      dplyr::rename(P=Pvalue)%>%
      #filter(P<0.1/700000) #If I add this step I could just get the same as pivotlonger "tableallSigExpWiseSNP"
      mutate(trt = (replace(trt,trt == "NH4+~10mM", "NH4+ 10 mM")))%>%
      mutate(trt = (replace(trt,trt == "NH4+~0.3mM", "NH4+ 0.3 mM")))%>%
      mutate(trt = (replace(trt,trt == "NH4+~3mM", "NH4+ 3 mM")))%>%
      mutate(trt = (replace(trt,trt == "NO3-~3mM", "NO3- 3 mM")))%>%
      mutate(row = row_number()) %>%
      pivot_wider(
        id_cols = c(row,SNP,Chr,Pos),
        id_expand = FALSE,
        names_from = trt,
        names_prefix = "",
        names_sep = "",
        names_glue = NULL,
        names_sort = FALSE,
        names_vary = "fastest",
        names_expand = FALSE,
        names_repair = "check_unique",
        values_from = P,
        values_fill = NULL
      )%>%
      dplyr::  select(-row) %>% 
      relocate('NO3- 3 mM',.before ='NH4+ 0.3 mM')%>% 
      relocate('NH4+ 0.3 mM',.before ='NH4+ 10 mM')%>% 
      relocate('NH4+ 10 mM',.before ='NH4+ 3 mM')
    
    MVP.Report(QTLforManhattan,
               plot.type=c("m"),
               multracks=T, #If this one is true, will get 2 plots:multitracks and multitraits, instead of n plot (for n traits)
               
               threshold=c(0.1/700000,0.1/700000),
               threshold.lty=c(2,2),
               threshold.lwd=c(2,2),
               threshold.col=c("black","black"),
               
               col=c("black","gray"),
               
               signal.col = "blue",
               signal.pch = 19,
               signal.cex=c(0.5,1),
               signal.line=0,
               
               amplify=T,
               
               highlight = c('SNP-4.23818648.','SNP-4.23792662.','SNP-12.25853862.','SNP-12.25890218.','SNP-12.25901203.'),
               highlight.col="red",
               highlight.pch = 17,
               highlight.cex =5, 
               
               ylim=c(4,20),
               band=3,
               
               file.type="pdf",
               
               memo=paste0("6B"),
               dpi=300,
               verbose = T)
    
  }  

}
}
+++++++++++++++++++++++++++++++++












+++++++++++++++++++++++++++++++++  