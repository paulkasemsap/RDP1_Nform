+++++++++++++++++++++++++++++++++
#### Mission Accomplished:  v 1.2 (May 2023).
#### Goal: Perform GWAS on rice growth phenotypes in response to different N forms
#### Last update: Based on the source file "gwas-GAPIT.r" last committed on 5/22/2023. Note the r/GAPIT verion used for data in the submitted manuscript is "GAPIT3_3.1.0".
#### Note ####
{#OK*---- Resources ----
# Genotype data from http://www.ricediversity.org/data/
## 700k - McCouch S, Wright M, Tung C-W, Maron L, McNally K, Fitzgerald M, Singh N, DeClerck G, Agosto Perez F, Korniliev P, Greenberg A, Nareda ME, Mercado SM, Harrington S, Shi Y, Branchini D, Kuser-Fal?ao, Leung H, Ebana K, Yano M, Eizenga G, McClung A, Mezey J. (2016) Open Access Resources for Genome Wide Association Mapping in Rice. Nature Comm, 7: 10532 doi 10.1038/ncomms10532 (http://rs-bt-mccouch4.biotech.cornell.edu/staged_data/HDRA-G6-4-RDP1-RDP2-NIAS.AGCT.hmp.txt.gz)
  
## 44k - Keyan Zhao, Chih-Wei Tung, Georgia C. Eizenga, Mark H. Wright, M. Liakat Ali, Adam H. Price, Gareth J. Norton, M. Rafiqul Islam, Andy Reynolds, Jason Mezey, Anna M. McClung, Carlos D. Bustamante & Susan R. McCouch (2011). Genome-wide association mapping reveals a rich genetic architecture of complex traits in Oryza sativa. Nat Comm 2:467 | DOI: 10.1038/ncomms1467, Published Online 13 Sep 2011. (http://www.ricediversity.org/data/sets/44kgwas)
}
+++++++++++++++++++++++++++++++++
#### Prerequisite ####
# GAPIT - Genomic Association and Prediction Integrated Tool
# Designed by Zhiwu Zhang
# GAPIT Source code written by Zhiwu Zhang, Alex Lipka, Feng Tian and You Tang
# http://zzlab.net/GAPIT/gapit_functions.txt
{#*---- Step 0: Import library and GAPIT functions; run this section each time to start R ----
# (Do this section only for new installation of R)
# As of GAPIT 3.3, seems like this step can be skipped?
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")
BiocManager::install("multtest")
  
install.packages("yaml")
install.packages("gplots")
install.packages("scatterplot3d")#The downloaded link at: http://cran.r-project.org/package=scatterplot3d
}
{#*---- Step 1.1: Load GAPIT ----
# Pick one option out of three to source functions from below; Choose github for latest IMO.
# Must note version used! Re-source if running again after a long time. Last source using GAPIT.Version="2022.4.16, GAPIT 3.1", Manual dated March 12, 2022; I will try to use the latest release from Github (2/21/23)
  
#From website
#source("http://zzlab.net/GAPIT/GAPIT.library.R"); doesn't need anymore from version 3.3
source("http://zzlab.net/GAPIT/gapit_functions.txt")
#source("http://www.zzlab.net/GAPIT/emma.txt")
  
#From Github
# install.packages("devtools")
devtools::install_github("jiabowang/GAPIT3",force=TRUE)
#("/Users/Zhiwu/Dropbox/Current/revolutionr/gapit/gapit_functions.txt")
  
#From Github archive 
# This one uses the latest release from Github at [https://github.com/jiabowang/GAPIT3/releases/tag/GAPIT3]?
install.packages("GAPIT3-GAPIT3.tar.gz", repos = NULL, type="source") #This file is dated 9/21/22 but the code is dated version 3.1.
# Examining differences between this one and previous version [https://github.com/jiabowang/GAPIT3/compare/GAPIT3.2...GAPIT3], seems like 3.1 is still the latest?!? -- Must check version description in the text, in session info, and confirm with the file name formatting.
}
{#*---- Step 1.2: Install & Load Other Required Packages ----
  packagelist <- c(
    "MASS",  # required for ginv
    "multtest",
    "gplots",
    "compiler", #required for cmpfun
    "scatterplot3d",
    "data.table",
    "GAPIT3",
    "rgl",
    "xfun"
  ) #Add packages into this list for ease of installation and loading
  #install.packages(packagelist) #Only use once
  lapply(packagelist, library, character.only = TRUE) #Load the required packages
  
  #Other packages that may fail to load automatically
  #install.packages("rgl_1.0.1.tar.gz", repos = NULL, type="source")
  #install.packages("xfun_0.37.tar.gz", repos = NULL, type="source")
  #install.Rtools(check = TRUE, check_r_update = TRUE, GUI = TRUE, ...)
}
getDTthreads() #check number of threads; should be 4 for higher efficiency
writeLines(capture.output(sessionInfo()), "sessionInfo.txt") #save session info and check for compatibility later, if needed
+++++++++++++++++++++++++++++++++
#### Working Code - GWAS ####
{#*---- Step 2: Set data directory and import files ----
workstation <- "C:/Users/kasemsap/Box/RiceGWAS2023/GAPIToutput/214Asmundson/MAF_selectionF_QN/" #Specify folder where all genotype and phenotype data are located for each population
setwd(workstation)
}
{#*---- Step 3: WHOLE PANEL - Load phenotypes and genotypes ----
  rdp1 <- read.csv("RDP1_Nform_Table_S2_QNbyOrganTrtGen_GWAS.csv", head = TRUE) #QN 
# See Step 3.1.2 (How we get the smaller genotype files) below in the archived section
# new files are <100 MB (Previously the whole file was > 2GB)
}
{#*---- Step 3.1: SUBPOP - Load phenotypes and genotypes ----
  # Prepare data files
  subpop <- read.csv("RDP1_Nform_accession.csv")
  rdp1 <- read.csv("RDP1_Nform_Table_S2_QNbyOrganTrtGen_GWAS.csv", head = TRUE) #QN

  # Merge df together. Keep only gen + phenotypes + subpop
  subpopdw  <- merge.data.table(rdp1,subpop, by = "gen")

  drops1 <- c("GSOR_ID","NSFTV_ID", "IRGC_ID" , "HDRA_ID", "varietal", "accession", "sampleset")
  subpopdata <- subpopdw [ , !(names(subpopdw) %in% drops1)]

  drops2 <- c("GSOR_ID","NSFTV_ID", "IRGC_ID" , "HDRA_ID", "subpop","accession", "sampleset")
  varietaldata <- subpopdw [ , !(names(subpopdw) %in% drops2)]

  # Filter and keep only each subpop together
  subpop_tej<- subset(subpopdata, subpop == "temperate-japonica" )
  subpop_trj<- subset(subpopdata, subpop == "tropical-japonica" )
  subpop_ind<- subset(subpopdata, subpop == "indica" )
  subpop_aus<- subset(subpopdata, subpop == "aus" )

  subpop_varind<- subset(varietaldata,  varietal == "indica" )
  subpop_varjap<- subset(varietaldata,  varietal == "japonica" )

  # Drop subpop column
  subpop_tej$subpop <- NULL
  subpop_trj$subpop <- NULL
  subpop_ind$subpop <- NULL
  subpop_aus$subpop <- NULL

  subpop_varind$varietal<- NULL
  subpop_varjap$varietal<- NULL
  #updated phenotype and genotype files, so that the genotype lists are matched.
  
  # See Step 3.1.2 (How we get the smaller genotype files) below in the archived section
  # new files are <100 MB (Previously the whole file was > 2GB)
}
{#*---- Step 4: Run GWAS with GAPIT #model selection F (SNP.MAF=0.05) ----
  {#*---- rdp1 ----
    pop= "rdp1"
    myY=  rdp1# change myY here
    
    {#*---- BLINK, FarmCPU, MLM - Random.model = T, mod-selection = F ----
      setwd(paste0(workstation,pop,"/Other"))
      
      ##Run with 700k SNP array
      writeLines(capture.output(sessionInfo()), paste0(pop,"_sessionInfo.txt")) #save session info and check for compatibility later, if needed
      
      start.time <- Sys.time()
      myGAPIT <- GAPIT(
        Y=myY,
        #G=myG,
        file.G =paste0("700kSNP_",pop,"_chr_"),
        file.Ext.G="txt",
        file.from=1,
        file.to=12,
        #file.path="", 
        
        PCA.total=3,
        #CV=myCV 
        SNP.MAF = 0.05,
        
        model= c("BLINK","FarmCPU","MLM"),
        Multiple_analysis= TRUE,
        Model.selection = FALSE,
        
        Random.model = TRUE
        
        
      )
      end.time <- Sys.time()
      time.taken <- end.time - start.time
      time.taken
      write(time.taken,paste0(pop,"_timetaken.txt"))
    }
    {#*---- GLM, Random.model = F, mod.selection = F ----
      setwd(paste0(workstation,pop,"/GLM"))  
      
      ##Run with 700k SNP array
      writeLines(capture.output(sessionInfo()), paste0(pop,"_sessionInfo.txt")) #save session info and check for compatibility later, if needed
      
      start.time <- Sys.time()
      myGAPIT <- GAPIT(
        Y=myY,
        #G=myG,
        file.G =paste0("700kSNP_",pop,"_chr_"),
        file.Ext.G="txt",
        file.from=1,
        file.to=12,
        #file.path="", 
        
        PCA.total=3,
        #CV=myCV 
        SNP.MAF = 0.05,
        
        model="GLM",
        Multiple_analysis= FALSE,
        Model.selection = FALSE,
        
        Random.model = FALSE
        
        
      )
      end.time <- Sys.time()
      time.taken <- end.time - start.time
      time.taken
      write(time.taken,paste0(pop,"_timetaken.txt"))
    }
  }
  {#*---- var ind - pc0 ----
    pop= "varind"
    myY= subpop_varind# change myY here
    
    {#*---- BLINK, FarmCPU, MLM - Random.model = T, mod-selection = F ----
      setwd(paste0(workstation,pop,"/Other"))
      
      ##Run with 700k SNP array
      writeLines(capture.output(sessionInfo()), paste0(pop,"_sessionInfo.txt")) #save session info and check for compatibility later, if needed
      
      start.time <- Sys.time()
      myGAPIT <- GAPIT(
        Y=myY,
        #G=myG,
        file.G =paste0("700kSNP_",pop,"_chr_"),
        file.Ext.G="txt",
        file.from=1,
        file.to=12,
        #file.path="", 
        
        PCA.total=0,
        #CV=myCV 
        SNP.MAF = 0.05,
        
        model= c("BLINK","FarmCPU","MLM"),
        Multiple_analysis= TRUE,
        Model.selection = FALSE,
        
        Random.model = TRUE
        
        
      )
      end.time <- Sys.time()
      time.taken <- end.time - start.time
      time.taken
      write(time.taken,paste0(pop,"_timetaken.txt"))
    }
    {#*---- GLM, Random.model = F, mod.selection = F ----
      setwd(paste0(workstation,pop,"/GLM"))  
      
      ##Run with 700k SNP array
      writeLines(capture.output(sessionInfo()), paste0(pop,"_sessionInfo.txt")) #save session info and check for compatibility later, if needed
      
      start.time <- Sys.time()
      myGAPIT <- GAPIT(
        Y=myY,
        #G=myG,
        file.G =paste0("700kSNP_",pop,"_chr_"),
        file.Ext.G="txt",
        file.from=1,
        file.to=12,
        #file.path="", 
        
        PCA.total=0,
        #CV=myCV 
        SNP.MAF = 0.05,
        
        model="GLM",
        Multiple_analysis= FALSE,
        Model.selection = FALSE,
        
        Random.model = FALSE
        
        
      )
      end.time <- Sys.time()
      time.taken <- end.time - start.time
      time.taken
      write(time.taken,paste0(pop,"_timetaken.txt"))
    }
  }
  {#*---- var jap - pc6 ----
    pop= "varjap"
    myY= subpop_varjap# change myY here
    
    {#*---- BLINK, FarmCPU, MLM - Random.model = T, mod-selection = F ----
      setwd(paste0(workstation,pop,"/Other"))
      
      ##Run with 700k SNP array
      writeLines(capture.output(sessionInfo()), paste0(pop,"_sessionInfo.txt")) #save session info and check for compatibility later, if needed
      
      start.time <- Sys.time()
      myGAPIT <- GAPIT(
        Y=myY,
        #G=myG,
        file.G =paste0("700kSNP_",pop,"_chr_"),
        file.Ext.G="txt",
        file.from=1,
        file.to=12,
        #file.path="", 
        
        PCA.total=6,
        #CV=myCV 
        SNP.MAF = 0.05,
        
        model= c("BLINK","FarmCPU","MLM"),
        Multiple_analysis= TRUE,
        Model.selection = FALSE,
        
        Random.model = TRUE
        
        
      )
      end.time <- Sys.time()
      time.taken <- end.time - start.time
      time.taken
      write(time.taken,paste0(pop,"_timetaken.txt"))
    }
    {#*---- GLM, Random.model = F, mod.selection = F ----
      setwd(paste0(workstation,pop,"/GLM"))  
      
      ##Run with 700k SNP array
      writeLines(capture.output(sessionInfo()), paste0(pop,"_sessionInfo.txt")) #save session info and check for compatibility later, if needed
      
      start.time <- Sys.time()
      myGAPIT <- GAPIT(
        Y=myY,
        #G=myG,
        file.G =paste0("700kSNP_",pop,"_chr_"),
        file.Ext.G="txt",
        file.from=1,
        file.to=12,
        #file.path="", 
        
        PCA.total=6,
        #CV=myCV 
        SNP.MAF = 0.05,
        
        model="GLM",
        Multiple_analysis= FALSE,
        Model.selection = FALSE,
        
        Random.model = FALSE
        
        
      )
      end.time <- Sys.time()
      time.taken <- end.time - start.time
      time.taken
      write(time.taken,paste0(pop,"_timetaken.txt"))
    }
  }
  {#*---- aus ----
    pop= "aus"
    myY= subpop_aus# change myY here
    
    {#*---- BLINK, FarmCPU, MLM - Random.model = T, mod-selection = F ----
      setwd(paste0(workstation,pop,"/Other"))
      
      ##Run with 700k SNP array
      writeLines(capture.output(sessionInfo()), paste0(pop,"_sessionInfo.txt")) #save session info and check for compatibility later, if needed
      
      start.time <- Sys.time()
      myGAPIT <- GAPIT(
        Y=myY,
        #G=myG,
        file.G =paste0("700kSNP_",pop,"_chr_"),
        file.Ext.G="txt",
        file.from=1,
        file.to=12,
        #file.path="", 
        
        PCA.total=0,
        #CV=myCV 
        SNP.MAF = 0.05,
        
        model= c("BLINK","FarmCPU","MLM"),
        Multiple_analysis= TRUE,
        Model.selection = FALSE,
        
        Random.model = TRUE
        
        
      )
      end.time <- Sys.time()
      time.taken <- end.time - start.time
      time.taken
      write(time.taken,paste0(pop,"_timetaken.txt"))
    }
    {#*---- GLM, Random.model = F, mod.selection = F ----
      setwd(paste0(workstation,pop,"/GLM"))  
      
      ##Run with 700k SNP array
      writeLines(capture.output(sessionInfo()), paste0(pop,"_sessionInfo.txt")) #save session info and check for compatibility later, if needed
      
      start.time <- Sys.time()
      myGAPIT <- GAPIT(
        Y=myY,
        #G=myG,
        file.G =paste0("700kSNP_",pop,"_chr_"),
        file.Ext.G="txt",
        file.from=1,
        file.to=12,
        #file.path="", 
        
        PCA.total=0,
        #CV=myCV 
        SNP.MAF = 0.05,
        
        model="GLM",
        Multiple_analysis= FALSE,
        Model.selection = FALSE,
        
        Random.model = FALSE
        
        
      )
      end.time <- Sys.time()
      time.taken <- end.time - start.time
      time.taken
      write(time.taken,paste0(pop,"_timetaken.txt"))
    }
  }
  {#*---- ind ----
    pop= "ind"
    myY= subpop_ind# change myY here
    
    {#*---- BLINK, FarmCPU, MLM - Random.model = T, mod-selection = F ----
      setwd(paste0(workstation,pop,"/Other"))
      
      ##Run with 700k SNP array
      writeLines(capture.output(sessionInfo()), paste0(pop,"_sessionInfo.txt")) #save session info and check for compatibility later, if needed
      
      start.time <- Sys.time()
      myGAPIT <- GAPIT(
        Y=myY,
        #G=myG,
        file.G =paste0("700kSNP_",pop,"_chr_"),
        file.Ext.G="txt",
        file.from=1,
        file.to=12,
        #file.path="", 
        
        PCA.total=0,
        #CV=myCV 
        SNP.MAF = 0.05,
        
        model= c("BLINK","FarmCPU","MLM"),
        Multiple_analysis= TRUE,
        Model.selection = FALSE,
        
        Random.model = TRUE
        
        
      )
      end.time <- Sys.time()
      time.taken <- end.time - start.time
      time.taken
      write(time.taken,paste0(pop,"_timetaken.txt"))
    }
    {#*---- GLM, Random.model = F, mod.selection = F ----
      setwd(paste0(workstation,pop,"/GLM"))  
      
      ##Run with 700k SNP array
      writeLines(capture.output(sessionInfo()), paste0(pop,"_sessionInfo.txt")) #save session info and check for compatibility later, if needed
      
      start.time <- Sys.time()
      myGAPIT <- GAPIT(
        Y=myY,
        #G=myG,
        file.G =paste0("700kSNP_",pop,"_chr_"),
        file.Ext.G="txt",
        file.from=1,
        file.to=12,
        #file.path="", 
        
        PCA.total=0,
        #CV=myCV 
        SNP.MAF = 0.05,
        
        model="GLM",
        Multiple_analysis= FALSE,
        Model.selection = FALSE,
        
        Random.model = FALSE
        
        
      )
      end.time <- Sys.time()
      time.taken <- end.time - start.time
      time.taken
      write(time.taken,paste0(pop,"_timetaken.txt"))
    }
  }
  {#*---- tej ----
    pop= "tej"
    myY= subpop_tej# change myY here
    
    {#*---- BLINK, FarmCPU, MLM - Random.model = T, mod-selection = F ----
      setwd(paste0(workstation,pop,"/Other"))
      
      ##Run with 700k SNP array
      writeLines(capture.output(sessionInfo()), paste0(pop,"_sessionInfo.txt")) #save session info and check for compatibility later, if needed
      
      start.time <- Sys.time()
      myGAPIT <- GAPIT(
        Y=myY,
        #G=myG,
        file.G =paste0("700kSNP_",pop,"_chr_"),
        file.Ext.G="txt",
        file.from=1,
        file.to=12,
        #file.path="", 
        
        PCA.total=0,
        #CV=myCV 
        SNP.MAF = 0.05,
        
        model= c("BLINK","FarmCPU","MLM"),
        Multiple_analysis= TRUE,
        Model.selection = FALSE,
        
        Random.model = TRUE
        
        
      )
      end.time <- Sys.time()
      time.taken <- end.time - start.time
      time.taken
      write(time.taken,paste0(pop,"_timetaken.txt"))
    }
    {#*---- GLM, Random.model = F, mod.selection = F ----
      setwd(paste0(workstation,pop,"/GLM"))  
      
      ##Run with 700k SNP array
      writeLines(capture.output(sessionInfo()), paste0(pop,"_sessionInfo.txt")) #save session info and check for compatibility later, if needed
      
      start.time <- Sys.time()
      myGAPIT <- GAPIT(
        Y=myY,
        #G=myG,
        file.G =paste0("700kSNP_",pop,"_chr_"),
        file.Ext.G="txt",
        file.from=1,
        file.to=12,
        #file.path="", 
        
        PCA.total=0,
        #CV=myCV 
        SNP.MAF = 0.05,
        
        model="GLM",
        Multiple_analysis= FALSE,
        Model.selection = FALSE,
        
        Random.model = FALSE
        
        
      )
      end.time <- Sys.time()
      time.taken <- end.time - start.time
      time.taken
      write(time.taken,paste0(pop,"_timetaken.txt"))
    }
  }
  {#*---- trj ----
    pop= "trj"
    myY= subpop_trj# change myY here
    
    {#*---- BLINK, FarmCPU, MLM - Random.model = T, mod-selection = F ----
      setwd(paste0(workstation,pop,"/Other"))
      
      ##Run with 700k SNP array
      writeLines(capture.output(sessionInfo()), paste0(pop,"_sessionInfo.txt")) #save session info and check for compatibility later, if needed
      
      start.time <- Sys.time()
      myGAPIT <- GAPIT(
        Y=myY,
        #G=myG,
        file.G =paste0("700kSNP_",pop,"_chr_"),
        file.Ext.G="txt",
        file.from=1,
        file.to=12,
        #file.path="", 
        
        PCA.total=0,
        #CV=myCV 
        SNP.MAF = 0.05,
        
        model= c("BLINK","FarmCPU","MLM"),
        Multiple_analysis= TRUE,
        Model.selection = FALSE,
        
        Random.model = TRUE
        
        
      )
      end.time <- Sys.time()
      time.taken <- end.time - start.time
      time.taken
      write(time.taken,paste0(pop,"_timetaken.txt"))
    }
    {#*---- GLM, Random.model = F, mod.selection = F ----
      setwd(paste0(workstation,pop,"/GLM"))  
      
      ##Run with 700k SNP array
      writeLines(capture.output(sessionInfo()), paste0(pop,"_sessionInfo.txt")) #save session info and check for compatibility later, if needed
      
      start.time <- Sys.time()
      myGAPIT <- GAPIT(
        Y=myY,
        #G=myG,
        file.G =paste0("700kSNP_",pop,"_chr_"),
        file.Ext.G="txt",
        file.from=1,
        file.to=12,
        #file.path="", 
        
        PCA.total=0,
        #CV=myCV 
        SNP.MAF = 0.05,
        
        model="GLM",
        Multiple_analysis= FALSE,
        Model.selection = FALSE,
        
        Random.model = FALSE
        
        
      )
      end.time <- Sys.time()
      time.taken <- end.time - start.time
      time.taken
      write(time.taken,paste0(pop,"_timetaken.txt"))
    }
  }
}
#### Archived code ####
# Need to run this once to prepare genotype files for each subpopulation
{#*---- Step 3.1.2 - Make genotype files for each subpopulation based ONLY on accessions that we have phenotype data for  ----
  library(bigreadr)
  #Genotype file (700k SNP array) was too large; need to split the file into smaller ones FOR EACH POPULATION
  {#*----  1. Make index for gen list ----
  genlist <- read.csv("RDP1_Nform_accession.csv")
  str(genlist)

  GWASgenlist<-genlist %>% dplyr::select(gen:varietal)
  str(GWASgenlist) #trimmed down to only 3 columns I need; this is an index genotype list to make myG for each subpopulation
  }
  {#*---- 2. Load genotype file  ----
  myG700k <- read.delim("HDRA-G6-4-RDP1-RDP2-NIAS.AGCT.hmp.txt", head = TRUE)
  summary(myG700k$chrom)

  # Select only NSFTV (RDP1)
  myG700kNSFTV<-myG700k %>% dplyr::select((rs.:QCcode)|starts_with("NSFTV"))

  head(myG700kNSFTV)
  dim(myG700kNSFTV)
  
  {#*---- Step 3.1: SUBPOP - Load phenotypes and genotypes ----
  # Prepare data files
  subpop <- read.csv("RDP1_Nform_accession.csv")
  rdp1 <- read.csv("RDP1_Nform_Table_S2_QNbyOrganTrtGen_GWAS.csv", head = TRUE) #QN

  # Merge df together. Keep only gen + phenotypes + subpop
  subpopdw  <- merge.data.table(rdp1,subpop, by = "gen")

  drops1 <- c("GSOR_ID","NSFTV_ID", "IRGC_ID" , "HDRA_ID", "varietal", "accession", "sampleset")
  subpopdata <- subpopdw [ , !(names(subpopdw) %in% drops1)]

  drops2 <- c("GSOR_ID","NSFTV_ID", "IRGC_ID" , "HDRA_ID", "subpop","accession", "sampleset")
  varietaldata <- subpopdw [ , !(names(subpopdw) %in% drops2)]

  # Filter and keep only each subpop together
  subpop_tej<- subset(subpopdata, subpop == "temperate-japonica" )
  subpop_trj<- subset(subpopdata, subpop == "tropical-japonica" )
  subpop_ind<- subset(subpopdata, subpop == "indica" )
  subpop_aus<- subset(subpopdata, subpop == "aus" )

  subpop_varind<- subset(varietaldata,  varietal == "indica" )
  subpop_varjap<- subset(varietaldata,  varietal == "japonica" )

  # Drop subpop column
  subpop_tej$subpop <- NULL
  subpop_trj$subpop <- NULL
  subpop_ind$subpop <- NULL
  subpop_aus$subpop <- NULL

  subpop_varind$varietal<- NULL
  subpop_varjap$varietal<- NULL
  #updated phenotype and genotype files, so that the genotype lists are matched.
}
  }
  {#*---- 3. Filter by accessions and write new files by chromosome for each population  ----
  myG700k_rdp1<-myG700kNSFTV %>% dplyr::select((rs.:QCcode),dplyr::any_of(GWASgenlist$gen)) 
    a = myG700k_rdp1 #change this file for each subpop
    for(i in unique(a$chrom)){write.table(a[a$chrom==i,],paste("700kSNP_rdp1_chr_",i,".txt",sep=""),sep="\t",quote=F,row.names=F)} #write 12 new files based on chromosome

  myG700k_aus<-myG700kNSFTV %>% dplyr::select((rs.:QCcode),dplyr::any_of(subpop_aus$gen))
    a = myG700k_aus #change this file for each subpop
    for(i in unique(a$chrom)){write.table(a[a$chrom==i,],paste("700kSNP_aus_chr_",i,".txt",sep=""),sep="\t",quote=F,row.names=F)} #write 12 new files based on chromosome
  myG700k_ind<-myG700kNSFTV %>% dplyr::select((rs.:QCcode),dplyr::any_of(subpop_ind$gen)) 
    a = myG700k_ind #change this file for each subpop
    for(i in unique(a$chrom)){write.table(a[a$chrom==i,],paste("700kSNP_ind_chr_",i,".txt",sep=""),sep="\t",quote=F,row.names=F)} #write 12 new files based on chromosome
  myG700k_tej<-myG700kNSFTV %>% dplyr::select((rs.:QCcode),dplyr::any_of(subpop_tej$gen)) 
    a = myG700k_tej #change this file for each subpop
    for(i in unique(a$chrom)){write.table(a[a$chrom==i,],paste("700kSNP_tej_chr_",i,".txt",sep=""),sep="\t",quote=F,row.names=F)} #write 12 new files based on chromosome
  myG700k_trj<-myG700kNSFTV %>% dplyr::select((rs.:QCcode),dplyr::any_of(subpop_trj$gen))
  a = myG700k_trj #change this file for each subpop
  for(i in unique(a$chrom)){write.table(a[a$chrom==i,],paste("700kSNP_trj_chr_",i,".txt",sep=""),sep="\t",quote=F,row.names=F)} #write 12 new files based on chromosome
  
  myG700k_varind<-myG700kNSFTV %>% dplyr::select((rs.:QCcode),dplyr::any_of(subpop_varind$gen)) 
    a = myG700k_varind #change this file for each subpop
    for(i in unique(a$chrom)){write.table(a[a$chrom==i,],paste("700kSNP_varind_chr_",i,".txt",sep=""),sep="\t",quote=F,row.names=F)} #write 12 new files based on chromosome
  myG700k_varjap<-myG700kNSFTV %>% dplyr::select((rs.:QCcode),dplyr::any_of(subpop_varjap$gen))
  a = myG700k_varjap #change this file for each subpop
  for(i in unique(a$chrom)){write.table(a[a$chrom==i,],paste("700kSNP_varjap_chr_",i,".txt",sep=""),sep="\t",quote=F,row.names=F)} #write 12 new files based on chromosome
  }
}