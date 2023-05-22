+++++++++++++++++++++++++++++++++
#### Mission Accomplished: v 1.3 (May 2023).
#### Goals:   
  # 1) Analyze data 
  # 2) Visualize data for publication
#### Last updated: Based on the source file "figure.r" last committed on 5/22/2023.
+++++++++++++++++++++++++++++++++
#### Prerequisite ####
{#OK*---- Install & Load Required Packages ----
  packagelist <- c(
    "tidyverse",
    "MASS",
    "lattice",
    "lme4",
    "car",
    "nlme",
    "hrbrthemes",
    "viridis", #https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html - robust color palette for color blindness
    "ggridges",
    "ggplot2",
    "patchwork",
    "gridExtra",
    "grid",
    "data.table",
    "cowplot",
    "Hmisc",
    "PerformanceAnalytics",
    "corrplot",
    "ggcorrplot",
    "psych"
  ) #Add packages into this list for ease of installation and loading
  #install.packages(packagelist) #Only use once
  lapply(packagelist, library, character.only = TRUE) #Load the required packages
}
{#OK*---- Prepare data ----
{#OK*---- Load raw data ----
  data <- read.csv("RDP1_Nform_data.csv", head = TRUE) #biomass data 

# examine data
  attach(data)
  head(data)
  dim(data) #must have 19756 samples (rows) and 15 variables (columns)
  summary(data)
  str(data)

# Make variables a factor
  data$accession<-as.factor(data$accession)
  data$GSOR_ID<-as.factor(data$GSOR_ID)
  data$IRGC_ID<-as.factor(data$IRGC_ID)
  data$NSFTV_ID<-as.factor(data$NSFTV_ID)
  data$HDRA_ID<-as.factor(data$HDRA_ID)
  data$GWAS_ID<-as.factor(data$GWAS_ID)
  data$subpop<-as.factor(data$subpop)
  data$varietal<-as.factor(data$varietal)
  data$trt<-as.factor(data$trt)
  data$exp<-as.factor(data$exp)
  data$subexp<-as.factor(data$subexp)
  data$block<-as.factor(data$block)
  data$organ<-as.factor(data$organ)
  str(data)
}
{#OK*---- Calculation - biomass traits for analysis ----
# Calculate whole plant total dw and fraction dw by organ for each plant replicate
  
  # Calculate whole plant total dw
  whole<-data %>%
    group_by(GWAS_ID,trt,exp,subexp,block,subpop,varietal)%>% 
    dplyr::summarise(dw=sum(dw))%>% 
    mutate (organ = "whole")%>%
    as.data.frame()
  
  wholedf<-rename(whole, totaldw = dw)
  
  # Calculate froot
  rootdf<- data %>% 
    group_by(GWAS_ID,trt,exp,subexp,block,subpop,varietal)%>%
    filter(organ=="root")%>%
    as.data.frame()%>%
    rename(rootdw = dw)

  mergeCols <- c("GWAS_ID", "trt", "exp", "subexp", "block", "subpop", "varietal")
  wholeroot<- merge.data.table(wholedf,rootdf, by = mergeCols)
  wholeroot$froot = wholeroot$rootdw/wholeroot$totaldw #use this df to test fraction of dw partitioned to root

  # Calculate shootrootratio
  shootdf<- data %>% 
    group_by(GWAS_ID,trt,exp,subexp,block,subpop,varietal)%>%
    filter(organ=="shoot")%>%
    as.data.frame()%>%
    rename(shootdw = dw)
  
  mergeCols <- c("GWAS_ID", "trt", "exp", "subexp", "block", "subpop", "varietal", "sampleset","accession","GSOR_ID","IRGC_ID","NSFTV_ID","HDRA_ID") 
  shootroot<- merge.data.table(shootdf,rootdf, by = mergeCols)
  shootroot$shootrootratio = shootroot$shootdw/shootroot$rootdw #use this df to test shoot to root ratio
}
}
getDTthreads() #check number of threads
writeLines(capture.output(sessionInfo()), "sessionInfo.txt") #save session info and check for compatibility later, if needed
+++++++++++++++++++++++++++++++++
#### Working Code: Preparing biomass data files for analyses and reporting ####
{#OK*---- DATA FILE "RDP1_Nform_gDW_mean" 1067 KB (precursor for GWAS inputs) - Calculate means across 4 biological replicates from 4 blocks  ----
  # Calculate mean shoot, root, total dw for each trt across 4 blocks; do it separately for each EXP
  expdwdf<-data %>% 
    group_by (GWAS_ID,trt,exp,subexp,subpop,varietal,organ) %>%
    dplyr::summarize (expdw=mean(dw, na.rm=TRUE))%>%
    as.data.frame()
    
  wholeexpdw<-whole %>% 
    group_by (GWAS_ID,trt,exp,subexp,subpop,varietal,organ)%>% 
    dplyr::summarize (expdw=mean(dw, na.rm=TRUE))
  wholeexpdw
  
  meandwbyorgan<-merge(expdwdf, wholeexpdw, all = TRUE)
  str(meandwbyorgan)
  summary(meandwbyorgan)
  
  # Calculate overall mean across 2 experiments
  meandwbyorganexpdf<-meandwbyorgan %>% 
    group_by (GWAS_ID,trt,organ, subpop, subexp,varietal) %>%
    dplyr::summarize (expdw=mean(expdw, na.rm=TRUE)) %>% 
    mutate(exp="exp")%>%
    as.data.frame()
  
  meandwbyorganexpdf$exp<-as.factor(meandwbyorganexpdf$exp)
  str(meandwbyorganexpdf)
  summary(meandwbyorganexpdf)
    
  totalmeandwbyorganexpdf<-merge(meandwbyorgan, meandwbyorganexpdf, all = TRUE)
  str(totalmeandwbyorganexpdf)
  summary(totalmeandwbyorganexpdf)
  
  RDP1_Nform_gDW_mean<-totalmeandwbyorganexpdf
  RDP1_Nform_gDW_mean %>% write.csv('RDP1_Nform_gDW_mean.csv', row.names = F) 
}  
{#OK*---- DATA FILE "RDP1_Nform_gDW_mean_ratio" 904 KB - Calculate biomass ratio to 3 mM NH4+  ----
# Calculate for each organ, ratio of meandw across 4 blocks of each N trt to that of 3 mM NH4+ for each accession
  RDP1_Nform_gDW_mean_ratio <-
  RDP1_Nform_gDW_mean%>% 
    pivot_wider(
      names_from = c(trt,organ,exp), 
      values_from = expdw
    ) %>%   
  #calculate ratio - there must be a way to do this more efficiently but I don't know that yet...

    mutate(`ratio_NH4+_0.3mM_root_1`= `NH4+_0.3mM_root_1`/`NH4+_3mM_root_1`)%>%
    mutate(`ratio_NH4+_0.3mM_shoot_1`= `NH4+_0.3mM_shoot_1`/`NH4+_3mM_shoot_1`)%>%
    mutate(`ratio_NH4+_0.3mM_whole_1`= `NH4+_0.3mM_whole_1`/`NH4+_3mM_whole_1`)%>%
    
    mutate(`ratio_NH4+_0.3mM_root_2`= `NH4+_0.3mM_root_2`/`NH4+_3mM_root_2`)%>%
    mutate(`ratio_NH4+_0.3mM_shoot_2`= `NH4+_0.3mM_shoot_2`/`NH4+_3mM_shoot_2`)%>%
    mutate(`ratio_NH4+_0.3mM_whole_2`= `NH4+_0.3mM_whole_2`/`NH4+_3mM_whole_2`)%>%
    
    mutate(`ratio_NH4+_0.3mM_root_exp`= `NH4+_0.3mM_root_exp`/`NH4+_3mM_root_exp`)%>%
    mutate(`ratio_NH4+_0.3mM_shoot_exp`= `NH4+_0.3mM_shoot_exp`/`NH4+_3mM_shoot_exp`)%>%
    mutate(`ratio_NH4+_0.3mM_whole_exp`= `NH4+_0.3mM_whole_exp`/`NH4+_3mM_whole_exp`)%>%
    
    mutate(`ratio_NH4+_10mM_root_1`= `NH4+_10mM_root_1`/`NH4+_3mM_root_1`)%>%
    mutate(`ratio_NH4+_10mM_shoot_1`= `NH4+_10mM_shoot_1`/`NH4+_3mM_shoot_1`)%>%
    mutate(`ratio_NH4+_10mM_whole_1`= `NH4+_10mM_whole_1`/`NH4+_3mM_whole_1`)%>%
    
    mutate(`ratio_NH4+_10mM_root_2`= `NH4+_10mM_root_2`/`NH4+_3mM_root_2`)%>%
    mutate(`ratio_NH4+_10mM_shoot_2`= `NH4+_10mM_shoot_2`/`NH4+_3mM_shoot_2`)%>%
    mutate(`ratio_NH4+_10mM_whole_2`= `NH4+_10mM_whole_2`/`NH4+_3mM_whole_2`)%>%
    
    mutate(`ratio_NH4+_10mM_root_exp`= `NH4+_10mM_root_exp`/`NH4+_3mM_root_exp`)%>%
    mutate(`ratio_NH4+_10mM_shoot_exp`= `NH4+_10mM_shoot_exp`/`NH4+_3mM_shoot_exp`)%>%
    mutate(`ratio_NH4+_10mM_whole_exp`= `NH4+_10mM_whole_exp`/`NH4+_3mM_whole_exp`)%>%
    
    mutate(`ratio_NO3-_3mM_root_1`= `NO3-_3mM_root_1`/`NH4+_3mM_root_1`)%>%
    mutate(`ratio_NO3-_3mM_shoot_1`= `NO3-_3mM_shoot_1`/`NH4+_3mM_shoot_1`)%>%
    mutate(`ratio_NO3-_3mM_whole_1`= `NO3-_3mM_whole_1`/`NH4+_3mM_whole_1`)%>%
    
    mutate(`ratio_NO3-_3mM_root_2`= `NO3-_3mM_root_2`/`NH4+_3mM_root_2`)%>%
    mutate(`ratio_NO3-_3mM_shoot_2`= `NO3-_3mM_shoot_2`/`NH4+_3mM_shoot_2`)%>%
    mutate(`ratio_NO3-_3mM_whole_2`= `NO3-_3mM_whole_2`/`NH4+_3mM_whole_2`)%>%
    
    mutate(`ratio_NO3-_3mM_root_exp`= `NO3-_3mM_root_exp`/`NH4+_3mM_root_exp`)%>%
    mutate(`ratio_NO3-_3mM_shoot_exp`= `NO3-_3mM_shoot_exp`/`NH4+_3mM_shoot_exp`)%>%
    mutate(`ratio_NO3-_3mM_whole_exp`= `NO3-_3mM_whole_exp`/`NH4+_3mM_whole_exp`)%>%
    
    dplyr::select(c(GWAS_ID,subpop, varietal,subexp,starts_with("ratio")))%>%
    
    pivot_longer(
      cols = starts_with("ratio"),
      names_to = c("prefix","n", "conc","organ","exp"),
      names_sep = "_",
      values_to = "ratio",
      values_drop_na = F
    )%>%
    dplyr::select(!(prefix))%>%
    mutate(`to_NH4+_3mM` = paste(n, conc, sep = '_'))%>%
  as.data.frame()%>%
    mutate(n=as.factor(n))%>%
    mutate(conc=as.factor(conc))%>%
    mutate(organ=as.factor(organ))%>%
    mutate(exp=as.factor(exp))%>%
    mutate(`to_NH4+_3mM`=as.factor(`to_NH4+_3mM`))

RDP1_Nform_gDW_mean_ratio%>%
    write.csv(., file=paste("RDP1_Nform_gDW_mean_ratio.csv",sep="."),quote=F,row.names=F)
}
+++++++++++++++++++++++++++++++++
#### Working Code: Data Reporting ####
{#OK*---- Load biomass data for reporting ----
# Prior to running this step, all "DATA FILE" sections in Analysis must be run first! The starting point of ALL three datasets used to make figures is "RDP1_Nform_gDW_mean" which is the means across 4 blocks in each experiment replicate.
  
  #1 Mean
  data_mean<- RDP1_Nform_gDW_mean
  data_mean$dw <-data_mean$expdw #add a new column so that I can use the old code that analyze data based on var name "dw"
  
  #2 Ratio to the reference N trt
  data_ratio<- RDP1_Nform_gDW_mean_ratio<- RDP1_Nform_gDW_mean%>% 
    pivot_wider(
      names_from = c(trt,organ,exp), 
      values_from = expdw
    )%>%
    mutate(`ratio_NH4+_0.3mM_root_1`= `NH4+_0.3mM_root_1`/`NH4+_3mM_root_1`)%>%
    mutate(`ratio_NH4+_0.3mM_shoot_1`= `NH4+_0.3mM_shoot_1`/`NH4+_3mM_shoot_1`)%>%
    mutate(`ratio_NH4+_0.3mM_whole_1`= `NH4+_0.3mM_whole_1`/`NH4+_3mM_whole_1`)%>%
    
    mutate(`ratio_NH4+_0.3mM_root_2`= `NH4+_0.3mM_root_2`/`NH4+_3mM_root_2`)%>%
    mutate(`ratio_NH4+_0.3mM_shoot_2`= `NH4+_0.3mM_shoot_2`/`NH4+_3mM_shoot_2`)%>%
    mutate(`ratio_NH4+_0.3mM_whole_2`= `NH4+_0.3mM_whole_2`/`NH4+_3mM_whole_2`)%>%
    
    mutate(`ratio_NH4+_0.3mM_root_exp`= `NH4+_0.3mM_root_exp`/`NH4+_3mM_root_exp`)%>%
    mutate(`ratio_NH4+_0.3mM_shoot_exp`= `NH4+_0.3mM_shoot_exp`/`NH4+_3mM_shoot_exp`)%>%
    mutate(`ratio_NH4+_0.3mM_whole_exp`= `NH4+_0.3mM_whole_exp`/`NH4+_3mM_whole_exp`)%>%
    
    mutate(`ratio_NH4+_10mM_root_1`= `NH4+_10mM_root_1`/`NH4+_3mM_root_1`)%>%
    mutate(`ratio_NH4+_10mM_shoot_1`= `NH4+_10mM_shoot_1`/`NH4+_3mM_shoot_1`)%>%
    mutate(`ratio_NH4+_10mM_whole_1`= `NH4+_10mM_whole_1`/`NH4+_3mM_whole_1`)%>%
    
    mutate(`ratio_NH4+_10mM_root_2`= `NH4+_10mM_root_2`/`NH4+_3mM_root_2`)%>%
    mutate(`ratio_NH4+_10mM_shoot_2`= `NH4+_10mM_shoot_2`/`NH4+_3mM_shoot_2`)%>%
    mutate(`ratio_NH4+_10mM_whole_2`= `NH4+_10mM_whole_2`/`NH4+_3mM_whole_2`)%>%
    
    mutate(`ratio_NH4+_10mM_root_exp`= `NH4+_10mM_root_exp`/`NH4+_3mM_root_exp`)%>%
    mutate(`ratio_NH4+_10mM_shoot_exp`= `NH4+_10mM_shoot_exp`/`NH4+_3mM_shoot_exp`)%>%
    mutate(`ratio_NH4+_10mM_whole_exp`= `NH4+_10mM_whole_exp`/`NH4+_3mM_whole_exp`)%>%
    
    mutate(`ratio_NO3-_3mM_root_1`= `NO3-_3mM_root_1`/`NH4+_3mM_root_1`)%>%
    mutate(`ratio_NO3-_3mM_shoot_1`= `NO3-_3mM_shoot_1`/`NH4+_3mM_shoot_1`)%>%
    mutate(`ratio_NO3-_3mM_whole_1`= `NO3-_3mM_whole_1`/`NH4+_3mM_whole_1`)%>%
    
    mutate(`ratio_NO3-_3mM_root_2`= `NO3-_3mM_root_2`/`NH4+_3mM_root_2`)%>%
    mutate(`ratio_NO3-_3mM_shoot_2`= `NO3-_3mM_shoot_2`/`NH4+_3mM_shoot_2`)%>%
    mutate(`ratio_NO3-_3mM_whole_2`= `NO3-_3mM_whole_2`/`NH4+_3mM_whole_2`)%>%
    
    mutate(`ratio_NO3-_3mM_root_exp`= `NO3-_3mM_root_exp`/`NH4+_3mM_root_exp`)%>%
    mutate(`ratio_NO3-_3mM_shoot_exp`= `NO3-_3mM_shoot_exp`/`NH4+_3mM_shoot_exp`)%>%
    mutate(`ratio_NO3-_3mM_whole_exp`= `NO3-_3mM_whole_exp`/`NH4+_3mM_whole_exp`)%>%
    
    dplyr::select(c(GWAS_ID,subpop,varietal,subexp,starts_with("ratio")))%>%
    
    pivot_longer(
      cols = starts_with("ratio"),
      names_to = c("prefix","n", "conc","organ","exp"),
      names_sep = "_",
      values_to = "ratio",
      values_drop_na = F
    )%>%
    dplyr::select(!(prefix))%>%
    mutate(`to_NH4+_3mM` = paste(n, conc, sep = '_'))%>%
    as.data.frame()%>%
    mutate(n=as.factor(n))%>%
    mutate(conc=as.factor(conc))%>%
    mutate(organ=as.factor(organ))%>%
    mutate(exp=as.factor(exp))%>%
    mutate(`to_NH4+_3mM`=as.factor(`to_NH4+_3mM`))
    
  #3 Mean, but in a "wide" format with each variable in separate columns
  databycolumn<- RDP1_Nform_gDW_mean%>% 
    group_by(GWAS_ID,exp,subpop, subexp,varietal,organ)%>% 
    dplyr::summarise(expdw=mean(expdw, na.rm=TRUE))%>%
    mutate (trt = as.factor("N_all"))%>%
    as.data.frame()%>%
    full_join(RDP1_Nform_gDW_mean)%>%
    pivot_wider(
      names_from = c(trt,organ,exp), 
      values_from = expdw
    )%>%
    as.data.frame() #Essentially the same numbers as "RDP1_Nform_gDW_mean", just 'wider' and contain averages across N trts; can visualize data better when root and shoot are in separate columns!

  attach(data_mean)
  attach(data_ratio)
  attach(databycolumn)
  
  head(data_mean)
  head(data_ratio)
  head(databycolumn)
  
  dim(data_mean) #12729 rows, 9 columns
  dim(data_ratio) #10530 rows, 10 columns
  dim(databycolumn) #390 rows, 49 columns
}
{#OK*---- FIGURES & TABLES following the narrative of the manuscript ----
{#OK*---- TABLE S1 - "RDP1_Nform_Table_S1_accession" 38 KB - Create a list of accessions  ----
  data %>%
    group_by(sampleset,accession,GSOR_ID,IRGC_ID,NSFTV_ID,HDRA_ID,GWAS_ID,subpop,varietal)%>% 
    distinct(GWAS_ID)%>% 
    write.csv(., file=paste("RDP1_Nform_Table_S1_accession_list.csv",sep="."),quote=F,row.names=F)
}
{#OK*---- TABLE S4 and TABLE S5 - Analyses of variance (ANOVA) and Tukey's HSD tests with R/emmeans ----
{#OK*---- ANALYSIS - Analyses of variance (ANOVA) for biomass trait ----
# Testing ANOVA for dw (shoot,root,total, froot), shootrootratio
  # Normal distribution - Shapiro wilk test is limited for sample size < 5000; https://stats.stackexchange.com/questions/446262/can-a-sample-larger-than-5-000-data-points-be-tested-for-normality-using-shapiro. However, when the sample sizes are large enough (>30-40), violations of to the normality assumption (for lm residuals) should not have significant impact on the conclusions (Ghasemi & Zahediasl, 2012; 10.5812/ijem.3505) and (Lumley et al, 2002; 10.1146/annurev.publhealth.23.100901.140546). Inspected Q-Q plots.
  {#OK*---- 1.1 ANOVA dw by genotype ----
    test_shoot<-lm(dw ~ trt*GWAS_ID+exp+exp/subexp+block, data = subset(data_allorgan,organ =="shoot"))
    test_root<-lm(dw ~ trt*GWAS_ID+exp+exp/subexp+block, data = subset(data_allorgan,organ =="root"))
    test_total<-lm(dw ~ trt*GWAS_ID+exp+exp/subexp+block, data = subset(data_allorgan,organ =="whole"))
    test_froot<-lm(froot ~ trt*GWAS_ID+exp+exp/subexp+block, data = wholeroot)
  }
  {#OK*---- 1.2 saving qqplots and anova tables ----
    pdf(file=paste("aov_dw_gen","root","qqplot.pdf",sep=".") )
    qqnorm(test_root$residuals,plot.it = T);  qqline(test_root$residuals, col = "red")
    dev.off()

    pdf(file=paste("aov_dw_gen","shoot","qqplot.pdf",sep=".") )
    qqnorm(test_shoot$residuals,plot.it = T);  qqline(test_shoot$residuals, col = "red")
    dev.off()

    pdf(file=paste("aov_dw_gen","whole","qqplot.pdf",sep=".") )
    qqnorm(test_total$residuals,plot.it = T);  qqline(test_total$residuals, col = "red")
    dev.off()

    pdf(file=paste("aov_dw_gen","froot","qqplot.pdf",sep=".") )
    qqnorm(test_froot$residuals,plot.it = T);  qqline(test_froot$residuals, col = "red")
    dev.off()
    }
  {#OK only need SS*----1.3 Tukey tests for mean comparisons ----
    require(multcomp)
    require(emmeans) #https://cran.r-project.org/web/packages/emmeans/vignettes/messy-data.html
    require(multcompView)

    shootresult <-  car::Anova(test_shoot)
    rootresult<-    car::Anova(test_root)
    totresult <-    car::Anova(test_total)
    frootresult<-    car::Anova(test_froot)

    {#OK*---- Make TABLE S4a - Extract P-value from ANOVA (Type II SS) ----
      Pr<-data.frame(
        Shoot = shootresult$"Pr(>F)",
        Root = rootresult$"Pr(>F)",
        Whole = totresult$"Pr(>F)",
        Fraction_root = frootresult$"Pr(>F)")
      row.names(Pr) <- c("Nitrogen", "Genotype", "Experiment replicate", "Block", "Genetic x Nitrogen", "Group within experiment replicate", "Residuals")

      t_Pr<-as.data.frame(t(Pr))
      tt_Pr_gen<-as.data.frame(t(t_Pr))
      tt_Pr_gen<-rownames_to_column(tt_Pr_gen, var = "Source of Variation")

      #Save as a csv
      tt_Pr_gen%>%
        write.csv(., file=paste("RDP1_Nform_Table_S4a_ANOVA_Pvalue_genotype.csv",sep="."),quote=F,row.names=F)
    }
    {#*---- Tukey HSD tests ----
    # {#*---- Shoot - No interaction ----
    #   {#*---- trt  ----
    #     tukeyemmeans = emmeans::emmeans(test_shoot,pairwise~trt)
    #     tukeyletter_trt<- multcomp::cld(tukeyemmeans[[1]],adjust="fdr",alpha=0.05,Letters=LETTERS, decreasing = T) 
    #   }
    #   {#*---- GWAS_ID  ----
    #     tukeyemmeans = emmeans::emmeans(test_shoot,pairwise~GWAS_ID)
    #     tukeyletter_GWAS_ID<- multcomp::cld(tukeyemmeans[[1]],adjust="fdr",alpha=0.05,Letters=LETTERS, decreasing = T) 
    #   }
    #   tukeyletter_shoot<-plyr::rbind.fill(tukeyletter_trt,tukeyletter_GWAS_ID)
    # 
    #   tukeyletter_shoot$organ<-"shoot"
    # }
    # {#*---- Root - interaction trt x gen ----
    #   # Find interaction terms for Tukey HSD following https://stats.stackexchange.com/questions/5250/multiple-comparisons-on-a-mixed-effects-model
    #   data_allorgan$trtGWAS_ID <- interaction(data_allorgan$trt, data_allorgan$GWAS_ID)
    #   test_root_interaction<-lm(dw ~ trtGWAS_ID+exp+exp/subexp+block, data = subset(data_allorgan,organ =="root")) #specify interaction term
    # 
    #   tukeyemmeans = emmeans::emmeans(test_root_interaction,pairwise~trtGWAS_ID) # computing all pairwise comparisons; seems like tukey is a default. #https://cran.r-project.org/web/packages/emmeans/vignettes/confidence-intervals.html. Confirmed when I look at the results. This used Tukey! 
    # 
    #   tukeyletter_root<- multcomp::cld(tukeyemmeans[[1]],adjust="fdr",alpha=0.05,Letters=LETTERS, decreasing = T) 
    #   tukeyletter_root$organ<-"root"
    # 
    #   tukeyletter_root<-plyr::rbind.fill(tukeyletter_trt,tukeyletter_GWAS_ID)
    # 
    #   tukeyletter_root$organ<-"root"
    # }
    # {#*---- Whole plant - No interaction ----
    #   {#OK*---- trt  ----
    #     tukeyemmeans = emmeans::emmeans(test_total,pairwise~trt)
    #     tukeyletter_trt<- multcomp::cld(tukeyemmeans[[1]],adjust="fdr",alpha=0.05,Letters=LETTERS, decreasing = T) 
    #   }
    #   {#OK*---- GWAS_ID  ----
    #     tukeyemmeans = emmeans::emmeans(test_total,pairwise~GWAS_ID)
    #     tukeyletter_GWAS_ID<- multcomp::cld(tukeyemmeans[[1]],adjust="fdr",alpha=0.05,Letters=LETTERS, decreasing = T) 
    #   }
    #   tukeyletter_whole<-plyr::rbind.fill(tukeyletter_trt,tukeyletter_GWAS_ID)
    # 
    #   tukeyletter_whole$organ<-"whole"
    # }
    # {#*---- Froot - No interaction ----
    #   {#OK*---- trt  ----
    #     tukeyemmeans = emmeans::emmeans(test_froot,pairwise~trt)
    #     tukeyletter_trt<- multcomp::cld(tukeyemmeans[[1]],adjust="fdr",alpha=0.05,Letters=LETTERS, decreasing = T) 
    #   }
    #   {#OK*---- GWAS_ID  ----
    #     tukeyemmeans = emmeans::emmeans(test_froot,pairwise~GWAS_ID)
    #     tukeyletter_GWAS_ID<- multcomp::cld(tukeyemmeans[[1]],adjust="fdr",alpha=0.05,Letters=LETTERS, decreasing = T) 
    #   }
    #   tukeyletter_froot<-plyr::rbind.fill(tukeyletter_trt,tukeyletter_GWAS_ID)
    # 
    #   tukeyletter_froot$organ<-"froot"
    # }
    # tukeyletter_GWAS_ID<-plyr::rbind.fill(tukeyletter_whole,tukeyletter_shoot,tukeyletter_root,tukeyletter_froot)
    # 
    # tukeyletter_GWAS_ID%>%
    #   write.csv(., file=paste("RDP1_Nform_Table_S5_emmeansTukey_GWAS_ID.csv",sep="."),quote=F,row.names=F)
    }
  }
  {#OK*---- 2.1 ANOVA dw by subpopulation  ----
    test_shoot<-lm(dw ~ trt*subpop+exp+exp/subexp+block, data = subset(data_allorgan,organ =="shoot"))
    test_root<-lm(dw ~ trt*subpop+exp+exp/subexp+block, data = subset(data_allorgan,organ =="root"))
    test_total<-lm(dw ~ trt*subpop+exp+exp/subexp+block, data = subset(data_allorgan,organ =="whole"))
    test_froot<-lm(froot ~ trt*subpop+exp+exp/subexp+block, data = wholeroot)
  }
  {#OK*---- 2.2 saving qqplots and anova tables ----
    pdf(file=paste("aov_dw_subpop","root","qqplot.pdf",sep=".") )
    qqnorm(test_root$residuals,plot.it = T);  qqline(test_root$residuals, col = "red")
    dev.off()

    pdf(file=paste("aov_dw_subpop","shoot","qqplot.pdf",sep=".") )
    qqnorm(test_shoot$residuals,plot.it = T);  qqline(test_shoot$residuals, col = "red")
    dev.off()

    pdf(file=paste("aov_dw_subpop","whole","qqplot.pdf",sep=".") )
    qqnorm(test_total$residuals,plot.it = T);  qqline(test_total$residuals, col = "red")
    dev.off()

    pdf(file=paste("aov_dw_subpop","froot","qqplot.pdf",sep=".") )
    qqnorm(test_froot$residuals,plot.it = T);  qqline(test_froot$residuals, col = "red")
    dev.off()

    # capture.output(summary(test_root), file = paste("aov_dw_subpop","root","lm","txt",sep="."))
    # capture.output(car::Anova(test_root, type ="II"), file = paste("aovII_dw_subpop","root","table","txt",sep="."))
    # capture.output(car::Anova(test_root, type ="III"), file = paste("aovIII_dw_subpop","root","table","txt",sep="."))
    #
    # capture.output(summary(test_shoot), file = paste("aov_dw_subpop","shoot","lm","txt",sep="."))
    # capture.output(car::Anova(test_shoot, type ="II"), file = paste("aovII_dw_subpop","shoot","table","txt",sep="."))
    # capture.output(car::Anova(test_shoot, type ="III"), file = paste("aovIII_dw_subpop","shoot","table","txt",sep="."))
    #
    # capture.output(summary(test_total), file = paste("aov_dw_subpop","whole","lm","txt",sep="."))
    # capture.output(car::Anova(test_total, type ="II"), file = paste("aovII_dw_subpop","whole","table","txt",sep="."))
    # capture.output(car::Anova(test_total, type ="III"), file = paste("aovIII_dw_subpop","whole","table","txt",sep="."))
    #
    # capture.output(summary(test_froot), file = paste("aov_dw_subpop","froot","lm","txt",sep="."))
    # capture.output(car::Anova(test_froot, type ="II"), file = paste("aovII_dw_subpop","froot","table","txt",sep="."))
    # capture.output(car::Anova(test_froot, type ="III"), file = paste("aovIII_dw_subpop","froot","table","txt",sep="."))
  }
  {#OK*---- 2.3 Tukey tests for mean comparisons ----
    require(multcomp)
    require(emmeans)
    require(multcompView)

    shootresult <-  car::Anova(test_shoot)
    rootresult<-    car::Anova(test_root)
    totresult <-    car::Anova(test_total)
    frootresult<-    car::Anova(test_froot)

    {#OK*---- Make TABLE S4b - Extract P-value from ANOVA (Type II SS) ----
      Pr<-data.frame(
        Shoot = shootresult$"Pr(>F)",
        Root = rootresult$"Pr(>F)",
        Whole = totresult$"Pr(>F)",
        Fraction_root = frootresult$"Pr(>F)")
      row.names(Pr) <- c("Nitrogen", "Subpopulation", "Experiment replicate", "Block", "Genetic x Nitrogen", "Group within experiment replicate", "Residuals")

      t_Pr<-as.data.frame(t(Pr))
      tt_Pr_subpop<-as.data.frame(t(t_Pr))
      tt_Pr_subpop<-rownames_to_column(tt_Pr_subpop, var = "Source of Variation")

      #Save as a csv
      tt_Pr_subpop%>%
        write.csv(., file=paste("RDP1_Nform_Table_S4b_ANOVA_Pvalue_subpop.csv",sep="."),quote=F,row.names=F)
    }
    {#OK*---- Tukey HSD tests ----
    {#OK*---- Shoot - No interaction ----
      {#OK*---- trt  ----
        tukeyemmeans = emmeans::emmeans(test_shoot,pairwise~trt)
        tukeyletter_trt<- multcomp::cld(tukeyemmeans[[1]],adjust="fdr",alpha=0.05,Letters=LETTERS, decreasing = T) 
      }
      {#OK*---- subpop  ----
        tukeyemmeans = emmeans::emmeans(test_shoot,pairwise~subpop)
        tukeyletter_subpop<- multcomp::cld(tukeyemmeans[[1]],adjust="fdr",alpha=0.05,Letters=LETTERS, decreasing = T) 
      }
      tukeyletter_shoot<-plyr::rbind.fill(tukeyletter_trt,tukeyletter_subpop)

      tukeyletter_shoot$organ<-"shoot"
    }
    {#OK*---- Root - No interaction ----
      {#OK*---- trt  ----
        tukeyemmeans = emmeans::emmeans(test_root,pairwise~trt)
        tukeyletter_trt<- multcomp::cld(tukeyemmeans[[1]],adjust="fdr",alpha=0.05,Letters=LETTERS, decreasing = T) 
      }
      {#OK*---- subpop  ----
        tukeyemmeans = emmeans::emmeans(test_root,pairwise~subpop)
        tukeyletter_subpop<- multcomp::cld(tukeyemmeans[[1]],adjust="fdr",alpha=0.05,Letters=LETTERS, decreasing = T) 
      }
      tukeyletter_root<-plyr::rbind.fill(tukeyletter_trt,tukeyletter_subpop)

      tukeyletter_root$organ<-"root"

    }
    {#OK*---- Whole plant - No interaction ----
      {#OK*---- trt  ----
        tukeyemmeans = emmeans::emmeans(test_total,pairwise~trt)
        tukeyletter_trt<- multcomp::cld(tukeyemmeans[[1]],adjust="fdr",alpha=0.05,Letters=LETTERS, decreasing = T) 
      }
      {#OK*---- subpop  ----
        tukeyemmeans = emmeans::emmeans(test_total,pairwise~subpop)
        tukeyletter_subpop<- multcomp::cld(tukeyemmeans[[1]],adjust="fdr",alpha=0.05,Letters=LETTERS, decreasing = T) 
      }
      tukeyletter_whole<-plyr::rbind.fill(tukeyletter_trt,tukeyletter_subpop)

      tukeyletter_whole$organ<-"whole"
    }
    {#OK*---- Froot - interaction trt x subpop ----
      # Find interaction terms for Tukey HSD following https://stats.stackexchange.com/questions/5250/multiple-comparisons-on-a-mixed-effects-model
      wholeroot$trtsubpop <- interaction(wholeroot$trt, wholeroot$subpop)
      test_froot<-lm(froot ~ trtsubpop+exp+exp/subexp+block, data = wholeroot) #specify interaction term

      tukeyemmeans = emmeans::emmeans(test_froot,pairwise~trtsubpop) # computing all pairwise comparisons; seems like tukey is a default. #https://cran.r-project.org/web/packages/emmeans/vignettes/confidence-intervals.html. Confirmed when I look at the results. This used Tukey! 

      tukeyletter_froot<- multcomp::cld(tukeyemmeans[[1]],adjust="fdr",alpha=0.05,Letters=LETTERS, decreasing = T) 
      tukeyletter_froot$organ<-"froot"
    }
    tukeyletter_subpop<-plyr::rbind.fill(tukeyletter_whole,tukeyletter_shoot,tukeyletter_root,tukeyletter_froot)

    tukeyletter_subpop%>%
      write.csv(., file=paste("RDP1_Nform_Table_S5_emmeansTukey_subpop.csv",sep="."),quote=F,row.names=F)
    }
  }
  {#OK*---- 3.1 ANOVA dw by varietal group  ----
    test_shoot<-lm(dw ~ trt*varietal+exp+exp/subexp+block, data = subset(data_allorgan,organ =="shoot"))
    test_root<-lm(dw ~ trt*varietal+exp+exp/subexp+block, data = subset(data_allorgan,organ =="root"))
    test_total<-lm(dw ~ trt*varietal+exp+exp/subexp+block, data = subset(data_allorgan,organ =="whole"))
    test_froot<-lm(froot ~ trt*varietal+exp+exp/subexp+block, data = wholeroot)
  }
  {#OK*---- 3.2 saving qqplots and anova tables ----
    pdf(file=paste("aov_dw_varietal","root","qqplot.pdf",sep=".") )
    qqnorm(test_root$residuals,plot.it = T);  qqline(test_root$residuals, col = "red")
    dev.off()

    pdf(file=paste("aov_dw_varietal","shoot","qqplot.pdf",sep=".") )
    qqnorm(test_shoot$residuals,plot.it = T);  qqline(test_shoot$residuals, col = "red")
    dev.off()

    pdf(file=paste("aov_dw_varietal","whole","qqplot.pdf",sep=".") )
    qqnorm(test_total$residuals,plot.it = T);  qqline(test_total$residuals, col = "red")
    dev.off()

    pdf(file=paste("aov_dw_varietal","froot","qqplot.pdf",sep=".") )
    qqnorm(test_froot$residuals,plot.it = T);  qqline(test_froot$residuals, col = "red")
    dev.off()

    # capture.output(summary(test_root), file = paste("aov_dw_varietal","root","lm","txt",sep="."))
    # capture.output(car::Anova(test_root, type="II"), file = paste("aovII_dw_varietal","root","table","txt",sep="."))
    # capture.output(car::Anova(test_root, type="III"), file = paste("aovIII_dw_varietal","root","table","txt",sep="."))
    #
    #
    # capture.output(summary(test_shoot), file = paste("aov_dw_varietal","shoot","lm","txt",sep="."))
    # capture.output(car::Anova(test_shoot, type="II"), file = paste("aovII_dw_varietal","shoot","table","txt",sep="."))
    # capture.output(car::Anova(test_shoot, type="III"), file = paste("aovIII_dw_varietal","shoot","table","txt",sep="."))
    #
    # capture.output(summary(test_total), file = paste("aov_dw_varietal","whole","lm","txt",sep="."))
    # capture.output(car::Anova(test_total, type="II"), file = paste("aovII_dw_varietal","whole","table","txt",sep="."))
    # capture.output(car::Anova(test_total, type="III"), file = paste("aovIII_dw_varietal","whole","table","txt",sep="."))
    #
    # capture.output(summary(test_froot), file = paste("aov_dw_varietal","froot","lm","txt",sep="."))
    # capture.output(car::Anova(test_froot, type="II"), file = paste("aovII_dw_varietal","froot","table","txt",sep="."))
    # capture.output(car::Anova(test_froot, type="III"), file = paste("aovIII_dw_varietal","froot","table","txt",sep="."))
  }
  {#OK*---- 3.3 Tukey tests for mean comparisons ----
    require(multcomp)
    require(emmeans)
    require(multcompView)

    shootresult <-  car::Anova(test_shoot)
    rootresult<-    car::Anova(test_root)
    totresult <-    car::Anova(test_total)
    frootresult<-    car::Anova(test_froot)

    {#OK*---- Make TABLE S4c - Extract P-value from ANOVA (Type II SS) ----
      Pr<-data.frame(
        Shoot = shootresult$"Pr(>F)",
        Root = rootresult$"Pr(>F)",
        Whole = totresult$"Pr(>F)",
        Fraction_root = frootresult$"Pr(>F)")
      row.names(Pr) <- c("Nitrogen", "Varietal", "Experiment replicate", "Block", "Genetic x Nitrogen", "Group within experiment replicate", "Residuals")

      t_Pr<-as.data.frame(t(Pr))
      tt_Pr_varietal<-as.data.frame(t(t_Pr))
      tt_Pr_varietal<-rownames_to_column(tt_Pr_varietal, var = "Source of Variation")

      #Save as a csv
      tt_Pr_varietal%>%
        write.csv(., file=paste("RDP1_Nform_Table_S4c_ANOVA_Pvalue_varietal.csv",sep="."),quote=F,row.names=F)
    }
    {#OK*---- Tukey HSD tests ----
    {#OK*---- Shoot - interaction trt x varietal----
      # Find interaction terms for Tukey HSD following https://stats.stackexchange.com/questions/5250/multiple-comparisons-on-a-mixed-effects-model
      data_allorgan$trtvarietal <- interaction(data_allorgan$trt, data_allorgan$varietal)
      test_shoot_interaction<-lm(dw ~ trtvarietal+exp+exp/subexp+block, data = subset(data_allorgan,organ =="shoot")) #specify interaction term

      tukeyemmeans = emmeans::emmeans(test_shoot_interaction,pairwise~trtvarietal) # computing all pairwise comparisons; seems like tukey is a default. #https://cran.r-project.org/web/packages/emmeans/vignettes/confidence-intervals.html. Confirmed when I look at the results. This used Tukey! 

      tukeyletter_shoot<- multcomp::cld(tukeyemmeans[[1]],adjust="fdr",alpha=0.05,Letters=LETTERS, decreasing = T) 

      tukeyletter_shoot$organ<-"shoot"
    }
    {#OK*---- Root - No interaction ----
      {#OK*---- trt  ----
        tukeyemmeans = emmeans::emmeans(test_root,pairwise~trt)
        tukeyletter_trt<- multcomp::cld(tukeyemmeans[[1]],adjust="fdr",alpha=0.05,Letters=LETTERS, decreasing = T) 
      }
      {#OK*---- varietal  ----
        tukeyemmeans = emmeans::emmeans(test_root,pairwise~varietal)
        tukeyletter_varietal<- multcomp::cld(tukeyemmeans[[1]],adjust="fdr",alpha=0.05,Letters=LETTERS, decreasing = T) 
      }
      tukeyletter_root<-plyr::rbind.fill(tukeyletter_trt,tukeyletter_varietal)

      tukeyletter_root$organ<-"root"

    }
    {#OK*---- Whole plant - No interaction ----
      {#OK*---- trt  ----
        tukeyemmeans = emmeans::emmeans(test_total,pairwise~trt)
        tukeyletter_trt<- multcomp::cld(tukeyemmeans[[1]],adjust="fdr",alpha=0.05,Letters=LETTERS, decreasing = T) 
      }
      {#OK*---- varietal  ----
        tukeyemmeans = emmeans::emmeans(test_total,pairwise~varietal)
        tukeyletter_varietal<- multcomp::cld(tukeyemmeans[[1]],adjust="fdr",alpha=0.05,Letters=LETTERS, decreasing = T) 
      }
      tukeyletter_whole<-plyr::rbind.fill(tukeyletter_trt,tukeyletter_varietal)

      tukeyletter_whole$organ<-"whole"
    }
    {#OK*---- Froot - interaction trt x varietal ----
      # Find interaction terms for Tukey HSD following https://stats.stackexchange.com/questions/5250/multiple-comparisons-on-a-mixed-effects-model
      wholeroot$trtvarietal <- interaction(wholeroot$trt, wholeroot$varietal)
      test_froot_interaction<-lm(froot ~ trtvarietal+exp+exp/subexp+block, data = wholeroot)  #specify interaction term

      tukeyemmeans = emmeans::emmeans(test_froot_interaction,pairwise~trtvarietal) # computing all pairwise comparisons; seems like tukey is a default. #https://cran.r-project.org/web/packages/emmeans/vignettes/confidence-intervals.html. Confirmed when I look at the results. This used Tukey!

      tukeyletter_froot<- multcomp::cld(tukeyemmeans[[1]],adjust="fdr",alpha=0.05,Letters=LETTERS, decreasing = T) 
      tukeyletter_froot$organ<-"froot"
    }
    tukeyletter_varietal<-plyr::rbind.fill(tukeyletter_whole,tukeyletter_shoot,tukeyletter_root,tukeyletter_froot)

    tukeyletter_varietal%>%
      write.csv(., file=paste("RDP1_Nform_Table_S5_emmeansTukey_varietal.csv",sep="."),quote=F,row.names=F)
    }
  }
  {#OK*---- 4.1 ANOVA ratio shoot to root  ----
    test_shootrootratio<-lm(shootrootratio ~ trt*GWAS_ID+exp+exp/subexp+block, data = subset(shootroot))
  }
  {#OK*---- 4.2 saving qqplots and car::Anova tables ----
    pdf(file=paste("aov_dw_gen","shootroot","qqplot.pdf",sep=".") )
    qqnorm(test_shootrootratio$residuals,plot.it = T);  qqline(test_shootrootratio$residuals, col = "red")
    dev.off()
    #shapiro.test(residuals(test_root))

    # capture.output(summary(test_shootrootratio), file = paste("aov_dw_gen","shootroot","lm","txt",sep="."))
    # capture.output(anova(test_shootrootratio, test="II"), file = paste("aovII_dw_gen","shootroot","table","txt",sep="."))
    # capture.output(anova(test_shootrootratio, test="III"), file = paste("aovIII_dw_gen","shootroot","table","txt",sep="."))
  }
  {#OK only SS*----4.3 Tukey tests for mean comparisons ----
    require(multcomp)
    require(emmeans)
    require(multcompView)

    shootrootresult <-  car::Anova(test_shootrootratio)

    {#OK*---- Make TABLE S4d - Extract P-value from ANOVA (Type II SS) ----
      Pr<-data.frame(Shoot_Root_ratio = shootrootresult$"Pr(>F)")
      row.names(Pr) <- c("Nitrogen", "Genotype", "Experiment replicate", "Block", "Genetic x Nitrogen", "Group within experiment replicate", "Residuals")

      t_Pr<-as.data.frame(t(Pr))
      tt_Pr_shootroot<-as.data.frame(t(t_Pr))
      tt_Pr_shootroot<-rownames_to_column(tt_Pr_shootroot, var = "Source of Variation")

      #Save as a csv
      tt_Pr_shootroot%>%
        write.csv(., file=paste("RDP1_Nform_Table_S4d_ANOVA_Pvalue_shootroot.csv",sep="."),quote=F,row.names=F)
    }
    {#*---- Tukey HSD tests ----
    # {#*---- Shoot to Root ratio - No interaction ----
    #   {#*---- trt  ----
    #     tukeyemmeans = emmeans::emmeans(test_shootrootratio,pairwise~trt)
    #     tukeyletter_trt<- multcomp::cld(tukeyemmeans[[1]],adjust="fdr",alpha=0.05,Letters=LETTERS, decreasing = T) 
    #   }
    #   {#*---- GWAS_ID  ----
    #     tukeyemmeans = emmeans::emmeans(test_shootrootratio,pairwise~GWAS_ID)
    #     tukeyletter_GWAS_ID<- multcomp::cld(tukeyemmeans[[1]],adjust="fdr",alpha=0.05,Letters=LETTERS, decreasing = T) 
    #   }
    #   tukeyletter_shootroot<-plyr::rbind.fill(tukeyletter_trt,tukeyletter_GWAS_ID)
    # 
    #   tukeyletter_shootroot$organ<-"shootroot"
    # }
    # tukeyletter_shootroot%>%
    #   write.csv(., file=paste("RDP1_Nform_Table_S5_emmeansTukey_shootroot_genotype.csv",sep="."),quote=F,row.names=F)
    }
  }
}
{#OK *---- ANALYSIS -  Analyses of variance (ANOVA) for biomass ratio to the reference N trt (3 mM NH4+) ----
  # **Can only test ANOVA for subpop and varietal because these are means across reps; no replicate for each genotype on its own**
  df_ratio<-
    RDP1_Nform_gDW_mean_ratio %>%
    dplyr::filter ( exp == ("1")|exp == ("2"))%>%
    mutate(trt_ratio =`to_NH4+_3mM`)


  # df_ratio<-df_ratio[complete.cases(df_ratio),] #remove NA's

  {#OK*----2.1 ANOVA ratio by subpopulation  ----
    test_shoot<-lm(ratio ~ trt_ratio *subpop+exp+exp/subexp, data = subset(df_ratio,organ =="shoot"))
    test_root<-lm(ratio ~ trt_ratio *subpop+exp+exp/subexp, data = subset(df_ratio,organ =="root"))
    test_total<-lm(ratio ~ trt_ratio *subpop+exp+exp/subexp, data = subset(df_ratio,organ =="whole"))
    }
  {#OK only QQ plots*---- 2.2 saving qqplots and anova tables ----
    pdf(file=paste("aov_ratio_subpop","root","qqplot.pdf",sep=".") )
    qqnorm(test_root$residuals,plot.it = T);  qqline(test_root$residuals, col = "red")
    dev.off()
    #shapiro.test(residuals(test_root))

    pdf(file=paste("aov_ratio_subpop","shoot","qqplot.pdf",sep=".") )
    qqnorm(test_shoot$residuals,plot.it = T);  qqline(test_shoot$residuals, col = "red")
    dev.off()
    #shapiro.test(residuals(test_shoot))

    pdf(file=paste("aov_ratio_subpop","whole","qqplot.pdf",sep=".") )
    qqnorm(test_total$residuals,plot.it = T);  qqline(test_total$residuals, col = "red")
    dev.off()
    #shapiro.test(residuals(test_total))

    #"Error in match.arg(test) : 'arg' should be one of “F”, “Chisq”", but can use results type II from the 2.3 instead.
    # capture.output(summary(test_root), file = paste("aov_ratio_subpop","root","lm","txt",sep="."))
    # capture.output(car::Anova(test_root, test="II"), file = paste("aovII_ratio_subpop","root","table","txt",sep="."))
    # capture.output(car::Anova(test_root, test="III"), file = paste("aovIII_ratio_subpop","root","table","txt",sep="."))
    # 
    # capture.output(summary(test_shoot), file = paste("aov_ratio_subpop","shoot","lm","txt",sep="."))
    # capture.output(car::Anova(test_shoot, test="II"), file = paste("aovII_ratio_subpop","shoot","table","txt",sep="."))
    # capture.output(car::Anova(test_shoot, test="III"), file = paste("aovIII_ratio_subpop","shoot","table","txt",sep="."))
    # 
    # capture.output(summary(test_total), file = paste("aov_ratio_subpop","whole","lm","txt",sep="."))
    # capture.output(car::Anova(test_total,test="III"), file = paste("aovII_ratio_subpop","whole","table","txt",sep="."))
    # capture.output(car::Anova(test_total,test="II"), file = paste("aovIII_ratio_subpop","whole","table","txt",sep="."))
  }
  {#OK*----2.3 Tukey tests for mean comparisons ----
    require(multcomp)
    require(emmeans)
    require(multcompView)

    shootresult <-  car::Anova(test_shoot)
    rootresult<-    car::Anova(test_root)
    totresult <-    car::Anova(test_total)

    {#OK*---- Make TABLE S4e - Extract P-value from ANOVA (Type II SS) ----
      Pr<-data.frame(
        Shoot = shootresult$"Pr(>F)",
        Root = rootresult$"Pr(>F)",
        Whole = totresult$"Pr(>F)")
      row.names(Pr) <- c("Nitrogen ratio", "Subpopulation", "Experiment replicate", "Genetic x Nitrogen", "Group within experiment replicate", "Residuals")

      t_Pr<-as.data.frame(t(Pr))
      tt_Pr_subpop_ratio<-as.data.frame(t(t_Pr))
      tt_Pr_subpop_ratio<-rownames_to_column(tt_Pr_subpop_ratio, var = "Source of Variation")

      #Save as a csv
      tt_Pr_subpop_ratio%>%
        write.csv(., file=paste("RDP1_Nform_Table_S4e_ANOVA_Pvalue_Subpopulation_ratio.csv",sep="."),quote=F,row.names=F)
    }
    {#OK*---- Tukey HSD tests ----
    {#OK*---- Shoot - N only ----
      {#OK*---- trt_ratio  ----
        tukeyemmeans = emmeans::emmeans(test_shoot,pairwise~trt_ratio)
        tukeyletter_trt_ratio<- multcomp::cld(tukeyemmeans[[1]],adjust="fdr",alpha=0.05,Letters=LETTERS, decreasing = T) 
      }
      tukeyletter_shoot<-tukeyletter_trt_ratio

      tukeyletter_shoot$organ<-"shoot"
    }
    {#OK*---- Root - N only  ----
      {#OK*---- trt_ratio  ----
        tukeyemmeans = emmeans::emmeans(test_root,pairwise~trt_ratio)
        tukeyletter_trt_ratio<- multcomp::cld(tukeyemmeans[[1]],adjust="fdr",alpha=0.05,Letters=LETTERS, decreasing = T) 
      }
      tukeyletter_root<-tukeyletter_trt_ratio

      tukeyletter_root$organ<-"root"

    }
    {#OK*---- Whole plant - N only  ----
      {#OK*---- trt_ratio  ----
        tukeyemmeans = emmeans::emmeans(test_total,pairwise~trt_ratio)
        tukeyletter_trt_ratio<- multcomp::cld(tukeyemmeans[[1]],adjust="fdr",alpha=0.05,Letters=LETTERS, decreasing = T) 
      }
      tukeyletter_whole<-tukeyletter_trt_ratio

      tukeyletter_whole$organ<-"whole"
    }

    tukeyletter_subpop_ratio<-plyr::rbind.fill(tukeyletter_whole,tukeyletter_shoot,tukeyletter_root)

    tukeyletter_subpop_ratio%>%
      write.csv(., file=paste("RDP1_Nform_Table_S5_emmeansTukey_subpop_ratio.csv",sep="."),quote=F,row.names=F)
    }
  }
  {#OK*----3.1 ANOVA ratio by varietal group  ----
    test_shoot<-lm(ratio ~ trt_ratio*varietal+exp+exp/subexp, data = subset(df_ratio,organ =="shoot"))
    test_root<-lm(ratio ~ trt_ratio*varietal+exp+exp/subexp, data = subset(df_ratio,organ =="root"))
    test_total<-lm(ratio ~ trt_ratio*varietal+exp+exp/subexp, data = subset(df_ratio,organ =="whole"))
  }
  {#OK only QQ plots*---- 3.2 saving qqplots and anova tables ----
    pdf(file=paste("aov_ratio_varietal","root","qqplot.pdf",sep=".") )
    qqnorm(test_root$residuals,plot.it = T);  qqline(test_root$residuals, col = "red")
    dev.off()
    #shapiro.test(residuals(test_root))

    pdf(file=paste("aov_ratio_varietal","shoot","qqplot.pdf",sep=".") )
    qqnorm(test_shoot$residuals,plot.it = T);  qqline(test_shoot$residuals, col = "red")
    dev.off()
    #shapiro.test(residuals(test_shoot))

    pdf(file=paste("aov_ratio_varietal","whole","qqplot.pdf",sep=".") )
    qqnorm(test_total$residuals,plot.it = T);  qqline(test_total$residuals, col = "red")
    dev.off()
    #shapiro.test(residuals(test_total))

    #"Error in match.arg(test) : 'arg' should be one of “F”, “Chisq”", but can use results type II from the 3.3 instead.
    # capture.output(summary(test_root), file = paste("aov_ratio_varietal","root","lm","txt",sep="."))
    # capture.output(car::Anova(test_root,test="II"), file = paste("aovII_ratio_varietal","root","table","txt",sep="."))
    # capture.output(car::Anova(test_root,test="III"), file = paste("aovIII_ratio_varietal","root","table","txt",sep="."))
    # 
    # capture.output(summary(test_shoot), file = paste("aov_ratio_varietal","shoot","lm","txt",sep="."))
    # capture.output(car::Anova(test_shoot, test="II"), file = paste("aovII_ratio_varietal","shoot","table","txt",sep="."))
    # capture.output(car::Anova(test_shoot, test="III"), file = paste("aovIII_ratio_varietal","shoot","table","txt",sep="."))
    # 
    # capture.output(summary(test_total), file = paste("aov_ratio_varietal","whole","lm","txt",sep="."))
    # capture.output(car::Anova(test_total, test="II"), file = paste("aovII_ratio_varietal","whole","table","txt",sep="."))
    # capture.output(car::Anova(test_total, test="III"), file = paste("aovIII_ratio_varietal","whole","table","txt",sep="."))
  }
  {#OK*----3.3 Tukey tests for mean comparisons ----
    require(multcomp)
    require(emmeans)
    require(multcompView)

    shootresult <-  car::Anova(test_shoot)
    rootresult<-    car::Anova(test_root)
    totresult <-    car::Anova(test_total)

    {#OK*---- Make TABLE S4e - Extract P-value from ANOVA (Type II SS) ----
      Pr<-data.frame(
        Shoot = shootresult$"Pr(>F)",
        Root = rootresult$"Pr(>F)",
        Whole = totresult$"Pr(>F)")
      row.names(Pr) <- c("Nitrogen ratio", "Varietal", "Experiment replicate", "Genetic x Nitrogen", "Group within experiment replicate", "Residuals")

      t_Pr<-as.data.frame(t(Pr))
      tt_Pr_varietal_ratio<-as.data.frame(t(t_Pr))
      tt_Pr_varietal_ratio<-rownames_to_column(tt_Pr_varietal_ratio, var = "Source of Variation")

      #Save as a csv
      tt_Pr_varietal_ratio%>%
        write.csv(., file=paste("RDP1_Nform_Table_S4e_ANOVA_Pvalue_varietal_ratio.csv",sep="."),quote=F,row.names=F)
    }
    {#OK*---- Tukey HSD tests ----
    {#OK*---- Shoot - interaction ----
      # Find interaction terms for Tukey HSD following https://stats.stackexchange.com/questions/5250/multiple-comparisons-on-a-mixed-effects-model
      df_ratio$trtvarietal <- interaction(df_ratio$trt, df_ratio$varietal)
      test_shoot_interaction<-lm(ratio ~ trtvarietal+exp+exp/subexp, data = subset(df_ratio,organ =="shoot")) #specify interaction term

      tukeyemmeans = emmeans::emmeans(test_shoot_interaction,pairwise~trtvarietal) # computing all pairwise comparisons; seems like tukey is a default. #https://cran.r-project.org/web/packages/emmeans/vignettes/confidence-intervals.html. Confirmed when I look at the results. This used Tukey!

      tukeyletter_trt_ratio<- multcomp::cld(tukeyemmeans[[1]],adjust="fdr",alpha=0.05,Letters=LETTERS, decreasing = T) 

      tukeyletter_shoot<-tukeyletter_trt_ratio

      tukeyletter_shoot$organ<-"shoot"
    }
    {#OK*---- Root - N only  ----
      {#OK*---- trt_ratio  ----
        tukeyemmeans = emmeans::emmeans(test_root,pairwise~trt_ratio)
        tukeyletter_trt_ratio<- multcomp::cld(tukeyemmeans[[1]],adjust="fdr",alpha=0.05,Letters=LETTERS, decreasing = T) 
      }
      tukeyletter_root<-tukeyletter_trt_ratio

      tukeyletter_root$organ<-"root"

    }
    {#OK*---- Whole plant - N only  ----
      {#OK*---- trt_ratio  ----
        tukeyemmeans = emmeans::emmeans(test_total,pairwise~trt_ratio)
        tukeyletter_trt_ratio<- multcomp::cld(tukeyemmeans[[1]],adjust="fdr",alpha=0.05,Letters=LETTERS, decreasing = T) 
      }
      tukeyletter_whole<-tukeyletter_trt_ratio

      tukeyletter_whole$organ<-"whole"
    }

    tukeyletter_varietal_ratio<-plyr::rbind.fill(tukeyletter_whole,tukeyletter_shoot,tukeyletter_root)

    tukeyletter_varietal_ratio%>%
      write.csv(., file=paste("RDP1_Nform_Table_S5_emmeansTukey_varietal_ratio.csv",sep="."),quote=F,row.names=F)
    }
  }
}
{#OK*----TABLE S4 Full ----
  tt_Pr_gen$Trait<-"DW_GWAS_ID"
  tt_Pr_subpop$Trait<-"DW_subpop"
  tt_Pr_varietal$Trait<-"DW_varietal"
  tt_Pr_shootroot$Trait<-"DW_shootrootratio"
  tt_Pr_subpop_ratio$Trait<-"DWratio_subpop"
  tt_Pr_varietal_ratio$Trait<-"DWratio_varietal"
  
  FullTableS4<-plyr::rbind.fill(
    tt_Pr_gen,
    tt_Pr_subpop,
    tt_Pr_varietal,
    tt_Pr_shootroot,
    tt_Pr_subpop_ratio,
    tt_Pr_varietal_ratio)
  
  FullTableS4%>%
    write.csv(., file=paste("RDP1_Nform_Table_S4_ANOVA_Pvalue.csv",sep="."),quote=F,row.names=F)
}
{#OK*----TABLE S5 Full ----
  head(tukeyletter_subpop)
  head(tukeyletter_varietal)
  head(tukeyletter_subpop_ratio)
  head(tukeyletter_varietal_ratio)
  
  tukeyletter_subpop$trait<-"DW_subpop"
  tukeyletter_varietal$trait<-"DW_varietal"
  tukeyletter_subpop_ratio$trait<-"DWratio_subpop"
  tukeyletter_varietal_ratio$trait<-"DWratio_varietal"
  
  FullTableS5<-plyr::rbind.fill(
    tukeyletter_subpop,
    tukeyletter_varietal,
    tukeyletter_subpop_ratio,
    tukeyletter_varietal_ratio)
  
  FullTableS5%>%
    write.csv(., file=paste("RDP1_Nform_Table_S5_emmeansTukey.csv",sep="."),quote=F,row.names=F)
}
}
{#OK*---- TABLE S2 - "RDP1_Nform_gDW_QNbyOrganTrtGen_GWAS" 92 KB - Compute quantile normalized values of residuals for GWAS ----
{#OK*---- ANALYSIS "Quantile normalization" to retain only variations from genetics and nitrogen treatments, and address non-normal distribution of biomass models residuals----
  {#OK*---- Create a function to compute QN ----
# Function to find QN in this section is adapted from https://www.mv.helsinki.fi/home/mjxpirin/GWAS_course/material/GWAS1.html, accessed in March 2023.
  inv.normalise <- function(x) {
    return(qnorm( (rank(x, na.last = "keep") - 0.5) / sum(!is.na(x))))}
  } 
  {#OK*---- Fit new linear models  ----
  data_allorgan <-
    data[,c(7:15)]%>% #keep only GWAS_ID as accession names
        rbind(whole,use.names=T,fill=F)%>% #add whole plant dw
        filter(!is.na(organ))
      str(data_allorgan)
      summary(data_allorgan)
{#OK*---- shoot ----
  organtested="shoot"
      test_dw<-lm(dw ~ exp+exp/subexp+block, data = subset(data_allorgan,organ==organtested)) # I should not include organ!
      r = residuals(test_dw)
      q = inv.normalise(r)
      
      # Compare distribution
      plot(density(scale(data_allorgan$dw)), col = "black", lwd = 2, xlim = c(-4,4), ylim = c(0,1), 
           xlab = "trait", main = "" ) #see distribution of dif organs
      
      plot(density(scale(subset(data_allorgan,organ==organtested)$dw)), col = "black", lwd = 2, xlim = c(-4,4), ylim = c(0,1), 
           xlab = "trait", main = "" )
      plot(density(scale(test_dw$model$dw)), col = "black", lwd = 2, xlim = c(-4,4), ylim = c(0,1), 
           xlab = "trait", main = "" ) #should be similar to the previous one
      
      lines(density(scale(test_dw$residuals)), col = "blue", lwd = 2) #residual
      lines(density(scale(q)), col = "red", lwd = 2) #inversed residual
      
      plot(subset(data_allorgan,organ==organtested)$dw, q, col = "black")
      
      plot(subset(data_allorgan,organ==organtested)$dw, q, col = c("black","red")[test_dw$model$exp])
      plot(subset(data_allorgan,organ==organtested)$dw, q, col = c("black","red")[test_dw$model$subexp])
      plot(subset(data_allorgan,organ==organtested)$dw, q, col = c("black","red","blue","orange")[test_dw$model$block])
      
      shootQN<-subset(data_allorgan,organ==organtested)%>%
        cbind(q)
      }
{#OK*---- root ----
        organtested="root"
        test_dw<-lm(dw ~ exp+exp/subexp+block, data = subset(data_allorgan,organ==organtested)) # I should not include organ!
        r = residuals(test_dw)
        q = inv.normalise(r)
        
        # Compare distribution
        plot(density(scale(data_allorgan$dw)), col = "black", lwd = 2, xlim = c(-4,4), ylim = c(0,1), 
             xlab = "trait", main = "" ) #see distribution of dif organs
        
        plot(density(scale(subset(data_allorgan,organ==organtested)$dw)), col = "black", lwd = 2, xlim = c(-4,4), ylim = c(0,1), 
             xlab = "trait", main = "" )
        plot(density(scale(test_dw$model$dw)), col = "black", lwd = 2, xlim = c(-4,4), ylim = c(0,1), 
             xlab = "trait", main = "" ) #should be similar to the previous one
        
        lines(density(scale(test_dw$residuals)), col = "blue", lwd = 2) #residual
        lines(density(scale(q)), col = "red", lwd = 2) #inversed residual
        
        plot(subset(data_allorgan,organ==organtested)$dw, q, col = "black")
        
        plot(subset(data_allorgan,organ==organtested)$dw, q, col = c("black","red")[test_dw$model$exp])
        plot(subset(data_allorgan,organ==organtested)$dw, q, col = c("black","red")[test_dw$model$subexp])
        plot(subset(data_allorgan,organ==organtested)$dw, q, col = c("black","red","blue","orange")[test_dw$model$block])
        
        rootQN<-subset(data_allorgan,organ==organtested)%>%
          cbind(q)
}
{#OK*---- whole ----
        organtested="whole"
        test_dw<-lm(dw ~ exp+exp/subexp+block, data = subset(data_allorgan,organ==organtested)) # I should not include organ!
        r = residuals(test_dw)
        q = inv.normalise(r)
        
        # Compare distribution
        plot(density(scale(data_allorgan$dw)), col = "black", lwd = 2, xlim = c(-4,4), ylim = c(0,1), 
             xlab = "trait", main = "" ) #see distribution of dif organs
        
        plot(density(scale(subset(data_allorgan,organ==organtested)$dw)), col = "black", lwd = 2, xlim = c(-4,4), ylim = c(0,1), 
             xlab = "trait", main = "" )
        plot(density(scale(test_dw$model$dw)), col = "black", lwd = 2, xlim = c(-4,4), ylim = c(0,1), 
             xlab = "trait", main = "" ) #should be similar to the previous one
        
        lines(density(scale(test_dw$residuals)), col = "blue", lwd = 2) #residual
        lines(density(scale(q)), col = "red", lwd = 2) #inversed residual
        
        plot(subset(data_allorgan,organ==organtested)$dw, q, col = "black")
        
        plot(subset(data_allorgan,organ==organtested)$dw, q, col = c("black","red")[test_dw$model$exp])
        plot(subset(data_allorgan,organ==organtested)$dw, q, col = c("black","red")[test_dw$model$subexp])
        plot(subset(data_allorgan,organ==organtested)$dw, q, col = c("black","red","blue","orange")[test_dw$model$block])
        
        wholeQN<-subset(data_allorgan,organ==organtested)%>%
          cbind(q)
      }
        
  }
}
{#OK*---- DATA FILE "RDP1_Nform_gDW_QNbyOrganTrtGen_GWAS" - Preparing data file for GWAS ----
   wholeQN%>%
      rbind(shootQN)%>%
      rbind(rootQN)%>%
      group_by(GWAS_ID,trt,organ)%>%
      add_tally()%>%
      group_by(GWAS_ID,trt,organ)%>%
      dplyr::summarize (dw_q=mean(q, na.rm=TRUE))%>%
      mutate(trt = str_replace(trt, paste0("NH4","\\+","_0.3mM"), "L"))%>%
      mutate(trt = str_replace(trt, paste0("NH4","\\+","_10mM"), "H"))%>%
      mutate(trt = str_replace(trt, paste0("NH4","\\+","_3mM"), "M"))%>%
      mutate(trt = str_replace(trt, paste0("NO3","\\-","_3mM"), "N"))%>%
      pivot_wider(
        names_from = c(trt,organ),
        values_from = dw_q,
        values_fill = NaN
      )%>%
      #mutate(across(where(is.numeric), round, 4))%>%
      rename(gen =GWAS_ID)%>%
      write.csv(., file=paste("RDP1_Nform_Table_S2_QNbyOrganTrtGen_GWAS.csv",sep="."),quote=F,row.names=F)
}
}
{#OK*---- TABLE S3 - "RDP1_Nform_Table_S3_model" 53 KB - PDF file made in MS word after examining all Q-Q plots from GWAS results ----
}
{#OK*---- FIGURE 1 Heritability ----
{#OK*---- Re-run ANOVA ----
    test_shoot<-lm(dw ~ trt*GWAS_ID+exp+exp/subexp+block, data = subset(data_allorgan,organ =="shoot"))
    test_root<-lm(dw ~ trt*GWAS_ID+exp+exp/subexp+block, data = subset(data_allorgan,organ =="root"))
    test_total<-lm(dw ~ trt*GWAS_ID+exp+exp/subexp+block, data = subset(data_allorgan,organ =="whole"))
    test_froot<-lm(froot ~ trt*GWAS_ID+exp+exp/subexp+block, data = wholeroot)
 
    shootresult <-  car::Anova(test_shoot)
    rootresult<-    car::Anova(test_root)
    totresult <-    car::Anova(test_total)
    frootresult<-    car::Anova(test_froot)
}
{#OK*---- Extract SS from ANOVA to calculate heritability ----
    SS<-data.frame(
    Root = rootresult$"Sum Sq",
    Shoot = shootresult$"Sum Sq",
    Whole = totresult$"Sum Sq",
    Fraction_root = frootresult$"Sum Sq")
  row.names(SS) <- c("Nitrogen", "Genetic", "Experiment replicate", "Block", "Genetic x Nitrogen", "Group within experiment replicate", "Residuals")

  t_SS<-as.data.frame(t(SS))
  t_SS<-rownames_to_column( t_SS, var = "trait")  
  
  #Calculate %manually and save as a csv
  t_SS%>%
    mutate(TotalSS = rowSums(across(where(is.numeric))))%>%
    mutate(`%Nitrogen` = `Nitrogen`/TotalSS)%>%
    mutate(`%Genetic` = `Genetic`/TotalSS)%>%
    mutate(`%Experiment replicate` = `Experiment replicate`/TotalSS)%>%
    mutate(`%Block` = `Block`/TotalSS)%>%
    mutate(`%GeneticxNitrogen` = `Genetic x Nitrogen`/TotalSS)%>%
    mutate(`%Group within experiment replicate` = `Group within experiment replicate`/TotalSS)%>%
    mutate(`%Residual` = `Residuals`/TotalSS)%>%
    dplyr::select(trait|starts_with("%"))%>%
    pivot_longer(cols= !trait, names_to = "component", values_to = "Proportion_Variance")%>%
    write.csv(., file=paste("RDP1_Nform_Figure_1_SourceData.csv",sep="."),quote=F,row.names=F)
}
{#OK*---- Make FIGURE 1 ----
  #Make a bar plot
  t_SS%>%
    pivot_longer(cols= !trait, names_to = "Component", values_to = "SS")%>%
    mutate(trait = str_replace_all(trait, paste0("_"), " "))%>%
  #pdf(file="RDP1_Nform_Figure_1.pdf", width = 6, height = 3)
  mutate(Component = factor(Component, levels=c(
    "Residuals","Group within experiment replicate", "Genetic x Nitrogen", "Block","Experiment replicate","Nitrogen","Genetic"))) %>%
  mutate(trait = factor(trait, levels=c(
    "Root", "Shoot", "Whole","Fraction root"))) %>%
  ggplot(aes(fill=Component, y=SS, x=trait)) + 
    geom_bar(position="fill", stat="identity",width = 0.8)+
    xlab("Biomass")+
    ylab("Proportion of variance")+
    scale_y_continuous(breaks = seq(0, 1, by=0.1))+
    theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = .5, linetype = 1),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          legend.background = element_rect(fill = "white"),
          legend.position="right",
          axis.title=element_text(size=9),
          axis.ticks = element_line(colour = "black", linewidth = .5),
          axis.text = element_text(color = "black", size = 9)
          )+
    scale_fill_viridis(discrete = TRUE)
  #dev.off() 
  
  ggsave("RDP1_Nform_Figure_1.pdf", dpi=300, height=3, width=6, units="in")  
  ggsave("RDP1_Nform_Figure_1.png", dpi=300, height=3, width=6, units="in")
}
}
{#OK*---- FIGURE S1 - Correlation matrix ----
  # https://github.com/caijun/ggcorrplot2
  datacorall<-databycolumn %>%
    dplyr::select(ends_with("exp"))%>%
    dplyr::select (!subexp)%>%
    rename_with(~ gsub("_exp", "", .x, fixed = TRUE))%>%
    rename_with(~ gsub("_", " ", .x, fixed = TRUE))
  
  datacorall<-as.matrix(datacorall)
  ct <- corr.test(datacorall, adjust = "none")
  corr <- ct$r
  p.mat <- ct$p
  
  # Run the following code to save the final corr fig as pdf
  pdf(file="RDP1_Nform_Figure_S1.pdf", width = 5, height = 5)
  corrplot(corr , method = "color", 
           type = "upper", order = "original", number.cex = .6,
           addCoef.col = "white", # Add coefficient of correlation
           tl.col = "black", tl.srt = 90, # Text label color and rotation
           # Combine with significance
           p.mat = p.mat, sig.level = 0.01, insig = "blank",
           # hide correlation coefficient on the principal diagonal
           diag = T,
           tl.cex = 0.8,
           number.digits = 2,outline = T, cl.cex = 0.8,cl.pos = 'b',cl.ratio = 0.1)
  dev.off()
  #This one can't use ggsave
}
{#OK*---- FIGURE 2 - Distribution of biomass by subpopulation ----
#RATIONALE: To demonstrate the effects of different components that influenced the biomass, I followed the ANOVA results. Blocks and Experiment replicates were significant. However, we expected 1) Blocks to be significant, because that's why we had blocks to begin with; to alleviate unaccounted variations, 2) Experiment replicates - ideally should NOT be different, but given the two replicates were conducted at different times;thus it also makes sense to yield different biomass. Should we visualize these two components? Decision: Since we focus on N and GENETIC, including blocks and experiment replicates in ANOVA was necessary and was proven effective in accounting for variations. But for the visualization purpose, taking averages across blocks and experiment replicates reduce 'visual' burden/complexity. (In theory, it would be best to visualize with individual plants as each data point; displaying all blocks and experiment replicates). So, we took averages across 4 blocks and 2 experiment replicates -> means across blocks and experiment replicates to make figures. *Therefore, the number of data points for jittering are 1/(4*2) less than what's supposed to be here! Tukey's test results (mean separations) are based on individual plants, not means.
{#OK*---- FIGURE 2ABC- Distribution of biomass by N ----
     fig2A <-
      data_mean %>%
      filter(organ == "whole") %>% 
      filter(exp == "exp") %>% 
      mutate(trt2 = str_replace_all(trt, paste0("_"), " "))%>% 
      mutate(trt2 = factor(trt2, levels=c("NO3- 3mM","NH4+ 3mM","NH4+ 0.3mM","NH4+ 10mM"))) %>%
      ggplot(aes(x=factor(trt2, levels=c("NO3- 3mM","NH4+ 3mM","NH4+ 0.3mM","NH4+ 10mM"))
                 ,y=dw))+ 
      annotate("text", x="NO3- 3mM", y=0.23, label= "a")+
      annotate("text", x="NH4+ 3mM", y=0.23, label= "a")+
      annotate("text", x="NH4+ 0.3mM", y=0.23, label= "a")+
      annotate("text", x="NH4+ 10mM", y=0.23, label= "b")+
      
      xlab("")+
       ylab("Whole plant biomass (g)")+
      scale_y_continuous(breaks = seq(0, 0.25, by=0.05), limits=c(0,0.25))+
    geom_violin(aes(fill=trt2),
                alpha=0.8,
                size=0.2)+ 
    geom_boxplot(
      width=0.3, color="black", 
      #outlier.shape = NA,
      outlier.shape = 4,
      alpha=0.5) +
      # geom_jitter(height = 0, width = 0.1, alpha=0.2, color="black", size=0.2)+
      theme(panel.background = element_rect(fill = "white", color = "black",linewidth = .5, linetype = 1),
            panel.border = element_rect(fill = 'transparent', color = "black",linewidth = .5, linetype = 1),
            panel.grid.major.x =element_line(color = "gray",
                                    linewidth = 0.25,
                                    linetype = 2),
            legend.position="none",
            axis.title.y=element_blank(),
            axis.title=element_blank(),
            axis.text = element_text(color = "black", size = 7)
      )+
      # scale_fill_viridis(discrete = TRUE)+
    scale_fill_manual(values=c("#FDE725FF","#20A387FF", "#39568CFF", "#440154FF"))+
      coord_flip()
    fig2A 
   
    fig2B <- data_mean %>%
      filter(organ == "shoot") %>% 
      filter(exp == "exp") %>% 
      mutate(trt2 = str_replace_all(trt, paste0("_"), " "))%>% 
      mutate(trt2 = factor(trt2, levels=c("NO3- 3mM","NH4+ 3mM","NH4+ 0.3mM","NH4+ 10mM"))) %>%
      ggplot(aes(x=factor(trt2, levels=c("NO3- 3mM","NH4+ 3mM","NH4+ 0.3mM","NH4+ 10mM"))
                 ,y=dw))+ 
      annotate("text", x="NO3- 3mM", y=0.23, label= "bc")+
      annotate("text", x="NH4+ 3mM", y=0.23, label= "a")+
      annotate("text", x="NH4+ 0.3mM", y=0.23, label= "b")+
      annotate("text", x="NH4+ 10mM", y=0.23, label= "c")+
      xlab("")+
       ylab("Shoot biomass (g)")+
      scale_y_continuous(breaks = seq(0, 0.250, by=0.05), limits=c(0,0.25))+
      geom_violin(aes(fill=trt2),
                  alpha=0.8,
                  size=0.2)+ 
      geom_boxplot(
        width=0.3, color="black", 
        #outlier.shape = NA,
        outlier.shape = 4,
        alpha=0.5) +
      # geom_jitter(height = 0, width = 0.1, alpha=0.2, color="black", size=0.2)+
      theme(panel.background = element_rect(fill = "white", color = "black",linewidth = .5, linetype = 1),
            panel.border = element_rect(fill = 'transparent', color = "black",linewidth = .5, linetype = 1),
            panel.grid.major.x =element_line(color = "gray",
                                             linewidth = 0.25,
                                             linetype = 2),
            legend.position="none",
            axis.title.y=element_blank(),
            axis.title=element_blank(),
            axis.text.y=element_blank(), 
            axis.text = element_text(color = "black", size = 7)
      )+
      # scale_fill_viridis(discrete = TRUE)+
      scale_fill_manual(values=c("#FDE725FF","#20A387FF", "#39568CFF", "#440154FF"))+      coord_flip()+
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())
    fig2B
    
    fig2C <- data_mean %>%
      filter(organ == "root") %>% 
      filter(exp == "exp") %>% 
      mutate(trt2 = str_replace_all(trt, paste0("_"), " "))%>% 
      mutate(trt2 = factor(trt2, levels=c("NO3- 3mM","NH4+ 3mM","NH4+ 0.3mM","NH4+ 10mM"))) %>%
      ggplot(aes(x=factor(trt2, levels=c("NO3- 3mM","NH4+ 3mM","NH4+ 0.3mM","NH4+ 10mM"))
                 ,y=dw))+ 
      annotate("text", x="NO3- 3mM", y=0.09, label= "a")+
      annotate("text", x="NH4+ 3mM", y=0.09, label= "b")+
      annotate("text", x="NH4+ 0.3mM", y=0.09, label= "a")+
      annotate("text", x="NH4+ 10mM", y=0.09, label= "c")+
      
      xlab("")+
       ylab("Root biomass (g)")+
      scale_y_continuous(breaks = seq(0, 0.1, by=0.02), limits=c(0,0.1))+
      geom_violin(aes(fill=trt2),
                  alpha=0.8,
                  size=0.2)+ 
      geom_boxplot(
        width=0.3, color="black", 
        #outlier.shape = NA,
        outlier.shape = 4,
        alpha=0.5) +
      # geom_jitter(height = 0, width = 0.1, alpha=0.2, color="black", size=0.2)+
      theme(panel.background = element_rect(fill = "white", color = "black",linewidth = .5, linetype = 1),
            panel.border = element_rect(fill = 'transparent', color = "black",linewidth = .5, linetype = 1),
            panel.grid.major.x =element_line(color = "gray",
                                             linewidth = 0.25,
                                             linetype = 2),
            legend.position="none",
            axis.title.y=element_blank(),
            axis.title=element_blank(),
            axis.text.y=element_blank(), 
            axis.text = element_text(color = "black", size = 7)
      )+
      # scale_fill_viridis(discrete = TRUE)+
      scale_fill_manual(values=c("#FDE725FF","#20A387FF", "#39568CFF", "#440154FF"))+      coord_flip()+
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())
    fig2C
    (fig2top<-(fig2A|fig2B|fig2C)) +plot_layout(ncol = 3)+plot_annotation(tag_levels = list(c("A", "B", "C")))

    ggsave("RDP1_Nform_Figure_2ABC.png", dpi=300, height=2, width=6, units="in")
    ggsave("RDP1_Nform_Figure_2ABC.pdf", dpi=300, height=2, width=6, units="in")
}  
{#OK*---- FIGURE 2DEF - Distribution of biomass by subpopulation ----
fig2D <-
  databycolumn%>%
    filter(!subpop == "admixed") %>% 
    # mutate(subpop = factor(subpop, levels=c("admixed","admixed-japonica","admixed-indica","aromatic","tropical-japonica","temperate-japonica","indica","aus"))) %>%
  ggplot(aes(x=factor(subpop, levels=c("admixed","admixed-japonica","aromatic","tropical-japonica","temperate-japonica","admixed-indica","indica","aus")), y=N_all_whole_exp))+ 
    # annotate("text", x="tropical-japonica", y=0.23, label= "B",fontface="bold") +
    xlab("Subpopulation")+
    ylab("Whole plant biomass (g)")+
    scale_y_continuous(breaks = seq(0, 0.25, by=0.05), 
                       # expand = c(0,0),
                       limits=c(0,0.25))+
    geom_violin(aes(fill=varietal),
                    alpha=0.8,
                    size=0.3)+ 
    geom_boxplot(
                 width=0.3, color="black", 
                 #outlier.shape = NA,
                 outlier.shape = 4,
                 alpha=0.5) +
    #geom_jitter(height = 0, width = 0.15, alpha=0.5, color="black", size=0.2)+
    theme(panel.background = element_rect(fill = "white", colour = "black",linewidth = .5, linetype = 1),
          panel.border = element_rect(fill = 'transparent', color = "black",linewidth = .5, linetype = 1),
          panel.grid.major.x =element_line(color = "gray",
                                           linewidth = 0.25,
                                           linetype = 2),
          legend.position="none",
          axis.title.y=element_blank(),
          axis.title=element_text(size=9),
          axis.text = element_text(color = "black", size = 7)
          )+
    scale_fill_viridis(discrete = TRUE) +
    annotate("text", x="admixed-japonica", y=0.23, label= "ab")+
    annotate("text", x="aromatic", y=0.23, label= "a")+
    annotate("text", x="tropical-japonica", y=0.23, label= "ab")+
    annotate("text", x="temperate-japonica", y=0.23, label= "bc")+
    annotate("text", x="admixed-indica", y=0.23, label= "de")+
    annotate("text", x="aus", y=0.23, label= "e")+
    annotate("text", x="indica", y=0.23, label= "cd")+
    
    coord_flip()+
  
    annotate("text", x="tropical-japonica", y=0.01, label= "(90)", size = 2) +
    annotate("text", x="temperate-japonica", y=0.01, label= "(95)", size = 2) +
    annotate("text", x="indica", y=0.01, label= "(84)", size = 2) +
    annotate("text", x="aus", y=0.01, label= "(60)", size = 2) +
    annotate("text", x="aromatic", y=0.01, label= "(14)", size = 2) +
    annotate("text", x="admixed-japonica", y=0.01, label= "(31)", size = 2) +
    annotate("text", x="admixed-indica", y=0.01, label= "(7)", size = 2)
    #annotate("text", x="admixed", y=0.01, label= "(9)", size = 2.5)
  fig2D
  
fig2E <- 
  databycolumn%>%
  filter(!subpop == "admixed") %>% 
  # mutate(subpop = factor(subpop, levels=c("admixed","admixed-japonica","admixed-indica","aromatic","tropical-japonica","temperate-japonica","indica","aus"))) %>%
  ggplot(aes(x=factor(subpop, levels=c("admixed","admixed-japonica","aromatic","tropical-japonica","temperate-japonica","admixed-indica","indica","aus")), y=N_all_shoot_exp))+
    # annotate("text", x="tropical-japonica", y=0.23, label= "B",fontface="bold") +
    xlab("Subpopulation")+
    ylab("Shoot biomass (g)")+
    scale_y_continuous(breaks = seq(0, 0.25, by=0.05),
                       # expand = c(0,0),
                       limits=c(0,0.25))+
  geom_violin(aes(fill=varietal),
              alpha=0.8,
              size=0.3)+ 
  geom_boxplot(
    width=0.3, color="black", 
    #outlier.shape = NA,
    outlier.shape = 4,
    alpha=0.5) +
    #geom_jitter(height = 0, width = 0.15, alpha=0.5, color="black", size=0.2)+
    theme(panel.background = element_rect(fill = "white", colour = "black",linewidth = .5, linetype = 1),
          panel.border = element_rect(fill = 'transparent', color = "black",linewidth = .5, linetype = 1),
          panel.grid.major.x =element_line(color = "gray",
                                           linewidth = 0.25,
                                           linetype = 2),
          legend.position="none",
          axis.title.y=element_blank(),
          axis.title=element_text(size=9),
          axis.text = element_text(color = "black", size = 7)
          )+
    scale_fill_viridis(discrete = TRUE) +
  annotate("text", x="admixed-japonica", y=0.23, label= "a")+
  annotate("text", x="aromatic", y=0.23, label= "a")+
  annotate("text", x="tropical-japonica", y=0.23, label= "a")+
  annotate("text", x="temperate-japonica", y=0.23, label= "b")+
  annotate("text", x="admixed-indica", y=0.23, label= "bc")+
  annotate("text", x="aus", y=0.23, label= "d")+
  annotate("text", x="indica", y=0.23, label= "c")+
    coord_flip()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

fig2F <- 
  databycolumn%>%
  filter(!subpop == "admixed") %>% 
  
  # mutate(subpop = factor(subpop, levels=c("admixed","admixed-japonica","admixed-indica","aromatic","tropical-japonica","temperate-japonica","indica","aus"))) %>%
  ggplot(aes(x=factor(subpop, levels=c("admixed","admixed-japonica","aromatic","tropical-japonica","temperate-japonica","admixed-indica","indica","aus")), y=N_all_root_exp))+
    # annotate("text", x="tropical-japonica", y=0.092, label= "C",fontface="bold") +
    xlab("Subpopulation")+
    ylab("Root biomass (g)")+
    scale_y_continuous(breaks = seq(0, 0.10, by=0.02),
                       #expand = c(0,0),
                       limits=c(0,0.10),)+
  geom_violin(aes(fill=varietal),
              alpha=0.8,
              size=0.3)+ 
  geom_boxplot(
    width=0.3, color="black", 
    #outlier.shape = NA,
    outlier.shape = 4,
    alpha=0.5) +
    #geom_jitter(height = 0, width = 0.15, alpha=0.5, color="black", size=0.2)+
      theme(panel.background = element_rect(fill = "white", colour = "black",linewidth = .5, linetype = 1),
            panel.border = element_rect(fill = 'transparent', color = "black",linewidth = .5, linetype = 1),
            panel.grid.major.x =element_line(color = "gray",
                                           linewidth = 0.25,
                                           linetype = 2),
          legend.position="none",
          axis.title.y=element_blank(),
          axis.title=element_text(size=9),
          axis.text = element_text(color = "black", size = 7)
          )+
    scale_fill_viridis(discrete = TRUE)+ 
    coord_flip()+
  annotate("text", x="admixed-japonica", y=0.09, label= "bc")+
  annotate("text", x="aromatic", y=0.09, label= "ab")+
  annotate("text", x="tropical-japonica", y=0.09, label= "cd")+
  annotate("text", x="temperate-japonica", y=0.09, label= "c")+
  annotate("text", x="admixed-indica", y=0.09, label= "d")+
  annotate("text", x="aus", y=0.09, label= "bc")+
  annotate("text", x="indica", y=0.09, label= "a")+
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

(fig2bottom<-(fig2D)|(fig2E)|(fig2F))+plot_layout(ncol = 3)+plot_annotation(tag_levels = list(c("D", "E", "F")))

ggsave("RDP1_Nform_Figure_2DEF.png", dpi=300, height=2, width=6, units="in")
ggsave("RDP1_Nform_Figure_2DEF.pdf", dpi=300, height=2, width=6, units="in")
}
  # fig2x<-grid.arrange(fig2A, fig2B, fig2C, fig2D,fig2E,fig2F, ncol=3) 
  # can't handle unequal plot size, https://cran.r-project.org/web/packages/gridExtra/vignettes/arrangeGrob.html 
  # https://ggplot2-book.org/arranging-plots.html?q=grid#laying-out-plots-side-by-side
  
#This configuration align the grids, BUT the bottom panel is made the same height as the top, which is not ideal... I went for this one for the manuscript
    (fig2<-(fig2A/fig2D)|(fig2B/fig2E)|(fig2C/fig2F))+plot_layout(ncol = 3)+plot_annotation(tag_levels = list(c("A","D","B","E", "C", "F")))

    ggsave("RDP1_Nform_Figure_2.png", dpi=300, height=4, width=6, units="in")
    ggsave("RDP1_Nform_Figure_2.pdf", dpi=300, height=4, width=6, units="in")
}
{#OK*---- FIGURE 3 - Fraction biomass partitioned to root by genetic x N ----
{#OK*---- FIGURE 3ABCD - Fraction biomass partitioned to root by varietal group x N ----
    f_labels <- data.frame(trt2 = c("NH4+ 0.3mM", "NH4+ 10mM", "NH4+ 3mM", "NO3- 3mM"), label = c("A", "B", "C", "D"))
    f_labelsdif <- data.frame(trt2 = c("NH4+ 0.3mM", "NH4+ 10mM", "NH4+ 3mM", "NO3- 3mM"), label= c("p<0.05", "p<0.05", "p<0.05", "p<0.05"))

    fig3top<-
  wholeroot%>%
      as.data.frame()%>%
      mutate(trt2 = str_replace_all(trt, paste0("_"), " "))%>%
      
      mutate(Varietal_Group = toupper(varietal))%>%
      mutate(Varietal_Group = factor(Varietal_Group, levels=c("INDICA","JAPONICA","ADMIXED"))) %>%
      rename_with(~ gsub("_", " ", .x, fixed = TRUE))%>%
      
      filter(!subpop=="admixed") %>%
      mutate(subpop = factor(subpop, levels=c("aus","indica","admixed-indica","temperate-japonica","tropical-japonica","aromatic","admixed-japonica","admixed"))) %>%
      ggplot(aes(x=`Varietal Group`, y=froot)) +
      
      geom_jitter(
        # aes(fill=`Varietal Group`),
        height = 0, 
        width = 0.1, 
        alpha=0.2, 
        color="black",
        size=0.2)+
      geom_violin(aes(fill=`Varietal Group`),
                  alpha=0.8,
                  size=0.3,
                  show.legend = TRUE)+
      
      facet_grid(.~trt2, 
                 #labeller = labeller(organ = organ.labs),
                 scales="free_y")+
      
      geom_boxplot(
        # aes(fill=`Varietal Group`),
        width=0.3, 
        color="black",
        alpha=0.5,
        outlier.shape = NA,
        # outlier.shape = 4,
        show.legend = TRUE
      ) +
      ylab("Fraction biomass partitioned to root")+
      xlab("Population")+  
      scale_y_continuous(breaks = seq(0, 0.8, by=0.1), limits=c(0,0.9))+
      theme(panel.background = element_rect(fill = "white", color = "black",linewidth = .5, linetype = 1),
            panel.border = element_rect(fill = 'transparent', color = "black",linewidth = .5, linetype = 1),
            panel.grid.major.y =element_line(color = "gray",
                                             linewidth = 0.25,
                                             linetype = 2),
            
            legend.text = element_text(size=7),
            legend.title = element_text(
              size=7,
              face = "bold"
              ),
            legend.key = element_rect(fill = "white", color = "black"),
           
            legend.position="top",
            # legend.box.background = element_rect(fill = "white", colour = "black"),
            # legend.position = c(.93, .83),
            # legend.justification = c("right", "top"),
            # legend.box.just = "right",
            # legend.margin = margin(3, 3, 3, 3),
          
            axis.title.y=element_text(size=7),   
            axis.title.x=element_blank(),
            axis.text.x = element_text(angle = 30, vjust = 1, hjust=1,size = 7),
            axis.text.y = element_text(size = 8),
            axis.text = element_text(color = "black"),
            
            strip.text.y = element_text(size=9,face="bold",angle = 270),
            # strip.text.x = element_text(size=9,face="bold"),
            strip.background = element_rect(colour="NA", fill="white")
      )+
      # guides(fill=
      #          guide_legend(direction = "horizontal",
      #                legend.key.height = unit(1, 'cm'), #change legend key height
      #                legend.key.width = unit(2, 'cm') #change legend key width
      #                )
      # )+
      # coord_flip()+
      scale_fill_viridis(discrete = TRUE)+      
      scale_color_viridis(discrete = TRUE)+
      stat_summary(fun=mean, color="blue",
                   geom="point",position=position_dodge(width=1), size = 1.5)+ 
      stat_summary(fun=mean, color="blue", aes(group=1),
                   geom="line", lwd=0.5, lty=1)+
      geom_text(x = "JAPONICA", y = 0.88, aes(label = label), data = f_labels )+
      geom_text(x = "JAPONICA",y = 0.75,aes(label = label), data = f_labelsdif,fontface="italic", size=3)
fig3top

    ggsave("RDP1_Nform_Figure_3ABCD.png", dpi=300, height=3, width=6, units="in")
    ggsave("RDP1_Nform_Figure_3ABCD.pdf", dpi=300, height=3, width=6, units="in")
} 
{#OK*---- FIGURE 3EFGH - Fraction biomass partitioned to root by subpop x N ----
    f_labels <- data.frame(trt2 = c("NH4+ 0.3mM", "NH4+ 10mM", "NH4+ 3mM", "NO3- 3mM"), label = c("E", "F", "G", "H"))
    
    fig3bottom<-
      wholeroot%>%
      
      as.data.frame()%>%
      mutate(trt2 = str_replace_all(trt, paste0("_"), " "))%>%
      
      mutate(Varietal_Group = toupper(varietal))%>%
      mutate(Varietal_Group = factor(Varietal_Group, levels=c("INDICA","JAPONICA","ADMIXED"))) %>%
      rename_with(~ gsub("_", " ", .x, fixed = TRUE))%>%
      
      filter(!subpop=="admixed") %>%
      mutate(subpop = factor(subpop, levels=c("aus","indica","admixed-indica","temperate-japonica","tropical-japonica","aromatic","admixed-japonica","admixed"))) %>%
      ggplot(aes(x=subpop, y=froot)) +
      
      geom_jitter(
        # aes(fill=`Varietal Group`),
        height = 0, 
        width = 0.1, 
        alpha=0.2, 
        color="black",
        size=0.2)+
      geom_violin(aes(fill=`Varietal Group`),
                  alpha=0.8,
                  size=0.3,
                  show.legend = TRUE)+
      
      facet_grid(.~trt2, 
                 #labeller = labeller(organ = organ.labs),
                 scales="free_y")+
      
      geom_boxplot(
        # aes(fill=`Varietal Group`),
        width=0.3, 
        color="black",
        alpha=0.5,
        outlier.shape = NA,
        # outlier.shape = 4,
        show.legend = TRUE
      ) +
      ylab("Fraction biomass partitioned to root")+
      xlab("Population")+  
      scale_y_continuous(breaks = seq(0, 0.8, by=0.1), limits=c(0,0.8))+
      theme(panel.background = element_rect(fill = "white", color = "black",linewidth = .5, linetype = 1),
            panel.border = element_rect(fill = 'transparent', color = "black",linewidth = .5, linetype = 1),
            panel.grid.major.y =element_line(color = "gray",
                                             linewidth = 0.25,
                                             linetype = 2),
            
            legend.text = element_text(size=7),
            legend.title = element_text(
              size=7,
              face = "bold"
              ),
            legend.key = element_rect(fill = "white", color = "black"),
           
            legend.position="none",
            # legend.box.background = element_rect(fill = "white", colour = "black"),
            # legend.position = c(.93, .83),
            # legend.justification = c("right", "top"),
            # legend.box.just = "right",
            # legend.margin = margin(3, 3, 3, 3),
          
            axis.title.y=element_text(size=7),   
            axis.title.x=element_blank(),
            axis.text.x = element_text(angle = 30, vjust = 1, hjust=1,size = 7),
            axis.text.y = element_text(size = 8),
            axis.text = element_text(color = "black"),
            
            strip.text.y = element_text(size=9,face="bold",angle = 270),
            strip.text.x = element_blank(),
            # strip.text.x = element_text(size=9,face="bold"),
            
            strip.background = element_rect(colour="NA", fill="white")
      )+
      # guides(fill=
      #          guide_legend(direction = "horizontal",
      #                legend.key.height = unit(1, 'cm'), #change legend key height
      #                legend.key.width = unit(2, 'cm') #change legend key width
      #                )
      # )+
      # coord_flip()+
      scale_fill_viridis(discrete = TRUE)+      
      scale_color_viridis(discrete = TRUE)+
      stat_summary(fun=mean, color="blue",
                   geom="point",position=position_dodge(width=1), size = 0.5)+ 
      stat_summary(fun=mean, color="blue", aes(group=1),
                   geom="line", lwd=0.5, lty=1)+
  
      geom_text(x = "admixed-japonica", y = 0.78, aes(label = label), data = f_labels)
    
    
    fig3bottom
    ggsave("RDP1_Nform_Figure_3EFGH.png", dpi=300, height=3, width=6, units="in")
    ggsave("RDP1_Nform_Figure_3EFGH.pdf", dpi=300, height=3, width=6, units="in")
}
  #This configuration align the grids, BUT the bottom panel is made the same height as the top, which is not ideal... I went for this one for the manuscript
  fig3<-(fig3top/fig3bottom)
  ggsave("RDP1_Nform_Figure_3.png", dpi=300, height=5, width=6, units="in")
  ggsave("RDP1_Nform_Figure_3.pdf", dpi=300, height=5, width=6, units="in")
}
{#OK*---- FIGURE 4 Biomass ratio to 3 mM NH4+ separated by organ and N ----
#Need to show 2 things: Differences betwen N trt AND shoot ratio of indica vs japonica at 3 mM NO3- because that is the only one with interaction (varietal x N)!
{#OK*---- FIGURE 4ABC Biomass ratio to 3 mM NH4+ separated by organ and N----
  f_labels <- data.frame(organ = c( "whole", "shoot","root"), label = c("A", "B", "C"))
  f_labelssigno3 <- data.frame(organ = c( "whole", "shoot","root"), label = c("a", "", "a"))
  f_labelssignh410 <- data.frame(organ = c( "whole", "shoot","root"), label = c("b", "", "b"))
  f_labelssignh40.3 <- data.frame(organ = c( "whole", "shoot","root"), label = c("a", "", "a"))
  
  fig4ABC<-
    data_ratio%>%
    # dplyr::filter ( exp == ("1")|exp == ("2"))%>%
    dplyr::filter ( exp == ("exp"))%>%
    dplyr::filter(!subpop=="admixed") %>%    

    mutate(trt_ratio =`to_NH4+_3mM`)%>%
    mutate(trt_ratio = str_replace_all(trt_ratio, paste0("_"), " "))%>%
    
    ggplot(aes(factor(trt_ratio),
               y=ratio))+
    
    xlab("")+
    ylab("Ratio of biomass (g/g biomass in NH4+ 3 mM)")+
    scale_y_continuous(breaks = seq(0, 4, by=0.5), limits=c(0,4.2))+
    geom_jitter(aes(fill="black",alpha=0.5),
                height = 0, 
                width = 0.1, 
                alpha=0.5, 
                color="black",
                size=0.1)+
    geom_violin(aes(fill=trt_ratio, alpha=1),
                alpha=0.8,
                size=0.3,
                color="black",
                draw_quantiles = NULL,
                trim = TRUE,
                show.legend = TRUE)+

    geom_boxplot(
      # aes(fill=trt_ratio, alpha=0.5),
      width=0.3, 
      color="black",
      alpha=0.5,
      outlier.shape = NA,
      # outlier.shape = 4,
      show.legend = TRUE)+
    geom_hline(yintercept=1,linetype=2, color = "red",size=0.5 )+
    
    # scale_fill_viridis(discrete = TRUE)+
    scale_fill_manual(values=c("white","#39568CFF","#440154FF", "#FDE725FF"))+
    
    theme(panel.background = element_rect(fill = "white", colour = "black"),legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    facet_grid(factor(organ, levels=c( "whole","shoot","root"))~.)+ #this arrangement works for setting the order of factors!!!! # https://www.statology.org/ggplot-facet-order/ #https://ggplot2.tidyverse.org/reference/facet_grid.html
    coord_flip()+
    theme(panel.background = element_rect(fill = "white", color = "black",linewidth = .5, linetype = 1),
          panel.border = element_rect(fill = 'transparent', color = "black",linewidth = .5, linetype = 1),
          panel.grid.major.x =element_blank(),
          
          legend.text = element_text(size=7),
          legend.title = element_text(
            size=7,
            face = "bold"
          ),
          legend.key = element_rect(fill = "white", color = "black"),
          
          legend.position="none",
          # legend.box.background = element_rect(fill = "white", colour = "black"),
          # legend.position = c(.93, .83),
          # legend.justification = c("right", "top"),
          # legend.box.just = "right",
          # legend.margin = margin(3, 3, 3, 3),
          
          axis.title.y=element_text(size=9),   
          axis.title.x=element_text(size=9),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5,size = 7),
          axis.text.y = element_text(size = 9),
          axis.text = element_text(color = "black"),
          
          strip.text.y = element_blank(),
          strip.text.x = element_blank(),
          # strip.text.x = element_text(size=9,face="bold"),
          
          strip.background = element_rect(colour="NA", fill="white")
    )+
    # theme(strip.text.x = element_text(size=9,face="bold"),
    #       strip.background = element_rect(colour="NA", fill="white"))+
    geom_text(x = "NO3- 3mM", y = 4.25, aes(label = label), data = f_labels) +
    geom_text(x = "NO3- 3mM", y = 0, aes(label = label), data = f_labelssigno3) +
    geom_text(x = "NH4+ 10mM", y = 0, aes(label = label), data = f_labelssignh410)+
    geom_text(x = "NH4+ 0.3mM", y = 0, aes(label = label), data = f_labelssignh40.3) 
  fig4ABC
  
  ggsave("RDP1_Nform_Figure_4ABC.png", dpi=300, height=6, width=3, units="in")
  ggsave("RDP1_Nform_Figure_4ABC.pdf", dpi=300, height=6, width=3, units="in")
}
{#OK*---- FIGURE 4D Inset plot showing differences between shoot ratio of INDICA vs JAPONICA from 3 mM NO3- ----
  # This one gets the one between 3 mM NO3- and shoot
  fig4D<-data_ratio%>%
    mutate(Varietal_Group = toupper(varietal))%>%
    mutate(Varietal_Group = factor(Varietal_Group, levels=c("INDICA","JAPONICA","ADMIXED"))) %>%    
    
  # dplyr::filter ( exp == ("1")|exp == ("2"))%>%
  dplyr::filter (exp == ("exp"))%>%
    
  dplyr::filter(!varietal=="admixed") %>%

  mutate(trt_ratio =`to_NH4+_3mM`)%>%
      dplyr::filter(trt_ratio=="NO3-_3mM") %>%
      mutate(trt_ratio = str_replace_all(trt_ratio, paste0("_"), " "))%>%      
    
      dplyr::filter(organ=="shoot") %>%
  
    mutate(subpop = factor(subpop, levels=c("aus","indica","admixed-indica","temperate-japonica","tropical-japonica","aromatic","admixed-japonica","admixed"))) %>%
  
  ggplot(aes(x=Varietal_Group,
         y=ratio))+
    
  xlab("")+
  ylab("NO3- 3 mM Shoot ratio")+
  scale_y_continuous(breaks = seq(0, 2.0, by=1.0), limits=c(0,2.0))+
    
    geom_violin(aes(alpha=1),
                fill="yellow",
                alpha=0.5,
                size=0.1,
                color="black",
                draw_quantiles = NULL,
                trim = TRUE,
                show.legend = TRUE)+
    # geom_jitter(
    #   height = 0, 
    #   width = 0.1, 
    #   alpha=0.2, 
    #   color="black",
    #   size=0.2)+
    geom_boxplot(
      aes(alpha=1),
      fill="yellow",
      width=0.4,
      size=0.2,
      color="black",
      alpha=1,
      outlier.shape = NA,
      # outlier.shape = 4,
      show.legend = TRUE)+
    geom_hline(yintercept=1,linetype=3, color = "red",size=0.5 )+
    
  scale_fill_viridis(discrete = TRUE)+
    theme(panel.background = element_rect(fill = "white", colour = "black"),legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
   
    theme(panel.background = element_rect(color = "black",linewidth = .5, linetype = 1),
          panel.border = element_rect(fill = 'transparent', color = "black",linewidth = .5, linetype = 1),
          panel.grid.major.x =element_blank(),
          
          legend.text = element_text(size=7),
          legend.title = element_text(
            size=7,
            face = "bold"
          ),
          legend.key = element_rect(fill = "white", color = "black"),
          
          legend.position="none",
          # legend.box.background = element_rect(fill = "white", colour = "black"),
          # legend.position = c(.93, .83),
          # legend.justification = c("right", "top"),
          # legend.box.just = "right",
          # legend.margin = margin(3, 3, 3, 3),
          
          axis.ticks.y = element_blank(),
          
          axis.title.y=element_text(size=7),   
          axis.title.x=element_text(size=7),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5,size = 7),
          axis.text.y = element_text(size = 8),
          axis.text = element_text(color = "black"),
          
          strip.text.y = element_blank(),
          strip.text.x = element_blank(),
          # strip.text.x = element_text(size=9,face="bold"),
          
          strip.background = element_rect(colour="NA", fill="white")
    )+
    # theme(panel.background = element_rect(fill = "white", colour = "black"),legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    facet_grid(trt_ratio~organ)+
    coord_flip()+
    # theme(strip.text.x = element_text(size=9,face="bold"),
    # strip.background = element_rect(colour="NA", fill="white"))+
    annotate("text", x=" ", y=0.7, label= "P = 0.0395", size = 2.5, fontface = 'italic')+
    annotate("text", x=" ", y=1.8, label= "D", size = 4)+
    scale_x_discrete(limits=rev)
    
ggsave("RDP1_Nform_Figure_4D.png", dpi=300, height=2, width=2, units="in")
ggsave("RDP1_Nform_Figure_4D.pdf", dpi=300, height=6, width=6, units="in")

# Put all components of figure 4 together
  ggdraw() +
  draw_plot(fig4ABC) +
  draw_plot(fig4D, x = 0.56, y = 0.436, width = 0.33, height = 0.22)

ggsave("RDP1_Nform_Figure_4.png", dpi=300, height=6, width=5, units="in")
ggsave("RDP1_Nform_Figure_4.pdf", dpi=300, height=6, width=5, units="in")
}
}
{#OK*----FIGURE S2 Connected dw by organ and population----
#Include this figure to dissect the variations between NO3- vs NH4+ growth at 3 mM between INDIVA vs JAPONICA (Follow up from Figure 4D). *Shoot at NO3- vs NH4+ shows interactions! 
  
{#OK*----Subpopulation----
  organ.labs  <- as_labeller(c(`whole` = "Whole plant", `shoot` = "Shoot", `root` = "Root")) #https://ggplot2.tidyverse.org/reference/as_labeller.html new codes from ggplot
  
  figsxb_comparenform_subpop<- data_mean %>%  
    # dplyr::filter ( exp == ("1")|exp == ("2"))%>%
    dplyr::filter (exp == ("exp"))%>%
    
  filter(trt== "NH4+_3mM" | trt== "NO3-_3mM") %>%
  filter(!grepl("admixed", varietal)) %>%
    mutate(trt2 = str_replace_all(trt, paste0("_"), " "))%>%
    
  ggplot(aes(trt2, dw)) +
  geom_violin(aes(fill=trt),alpha=0.8,  size=0.2)+
  geom_boxplot(
    # aes(fill=trt, alpha=1),
    alpha=0.5,width=0.3,outlier.shape = NA)+
    geom_line(aes(group=GWAS_ID), alpha = 0.1)+
  
    facet_grid(organ~factor(subpop, levels = c("aus","indica","admixed-indica","temperate-japonica","tropical-japonica","aromatic", "admixed-japonica")),labeller = labeller(organ =organ.labs),scales="free_y")+ #https://ggplot2.tidyverse.org/reference/labellers.html
  ylab("Biomass (g)")+
  xlab("")+
  
    theme(axis.text = element_text(size = 9),
          panel.border = element_rect(fill = 'transparent', color = "black",linewidth = .5, linetype = 1),
        axis.title.x = element_text(size=9),
        axis.title.y = element_text(size=9),
        legend.text = element_blank(),
        legend.title = element_blank(),
        legend.position="none",
        axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))+
  
    scale_fill_viridis(discrete = TRUE)+
    scale_color_viridis(discrete = TRUE)+
  
    stat_summary(fun=mean, colour="blue",
               geom="point",position=position_dodge(width=0.5))+ 
    stat_summary(fun=mean, colour="blue", aes(group=1),
               geom="line", lwd=1, lty=1)+
    theme(strip.text.y = element_text(size=9,face="bold"),strip.text.x = element_text(size=5,face="bold"),strip.background = element_rect(colour="NA", fill="white"))+
    theme(panel.background = element_rect(fill = "white", colour = "black"))
  figsxb_comparenform_subpop
  
ggsave("RDP1_Nform_Figure_S2.png", dpi=300, height=5, width=6, units="in")
ggsave("RDP1_Nform_Figure_S2.pdf", dpi=300, height=5, width=6, units="in")
}
# {#OK - Not as clear as by subpopulation*----Varietal group ----
#   data_mean %>%  
#     # dplyr::filter ( exp == ("1")|exp == ("2"))%>%
#     dplyr::filter (exp == ("exp"))%>%
#     
#     mutate(Varietal_Group = toupper(varietal))%>%
#     mutate(Varietal_Group = factor(Varietal_Group, levels=c("INDICA","JAPONICA","ADMIXED"))) %>%    
#     
#     filter(trt== "NH4+_3mM" | trt== "NO3-_3mM") %>%
#     filter(!grepl("admixed", subpop)) %>%
#     
#     mutate(trt2 = str_replace_all(trt, paste0("_"), " "))%>%
#     
#     ggplot(aes(trt2, dw)) +
#     geom_violin(aes(fill=trt),alpha=.5)+
#     geom_boxplot(aes(fill=trt, alpha=1),width=0.3,outlier.shape = NA)+
#     geom_line(aes(group=GWAS_ID), alpha = 0.1)+
#     
#     facet_grid(factor(organ, levels = c("whole","shoot","root"))~factor(Varietal_Group, levels=c("INDICA","JAPONICA")), labeller = labeller(organ = organ.labs),scales="free_y")+
#     ylab("Biomass (g)")+
#     xlab("")+
#     
#     theme(axis.text = element_text(size = 9),
#           axis.title.x = element_text(size=9),
#           axis.title.y = element_text(size=9),
#           legend.text = element_blank(),
#           legend.title = element_blank(),
#           legend.position="none",
#           axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))+
#     
#     scale_fill_viridis(discrete = TRUE)+
#     scale_color_viridis(discrete = TRUE)+
#     
#     stat_summary(fun=mean, colour="white",
#                  geom="point",position=position_dodge(width=0.75))+ 
#     stat_summary(fun=mean, colour="blue", aes(group=1),
#                  geom="line", lwd=1, lty=1)+
#     theme(strip.text.y = element_text(size=9,face="bold"),strip.text.x = element_text(size=5,face="bold"),strip.background = element_rect(colour="NA", fill="white"))+
#     theme(panel.background = element_rect(fill = "white", colour = "black"))
#   
#   # plot_grid(figx_comparenform , labels = "A")
#   ggsave("figsxa_comparenform_varietal.png", dpi=300, height=5, width=3, units="in")
#   ggsave("figsxa_comparenform_varietal.pdf", dpi=300, height=5, width=6, units="in")
# }
}
}