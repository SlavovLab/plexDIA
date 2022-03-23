
source("functions_parameters.R")
library(reticulate)
source_python("py_isoEnv.py")
library(dplyr)

'%!in%' <- function(x,y)!('%in%'(x,y))


#Quick data functions:
pD_channel <- function(df){
  
  df <- df %>% dplyr::mutate("channel_name" = ifelse(grepl("-0|mTRAQ0", Modified.Sequence), "mTRAQ0",
                                              ifelse(grepl("-4|mTRAQ4", Modified.Sequence), "mTRAQ4", "mTRAQ8")))
  
  return(df)
}

pD_seqcharge <- function(df){
  df$seqcharge <- paste0(df$Modified.Sequence, df$Precursor.Charge)
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ0\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ-K-0\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ-n-0\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ4\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ-K-4\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ-n-4\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ8\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ-K-8\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ-n-8\\)")
  return(df)
}

pD_seqcharge2 <- function(df){
  df$seqcharge <- paste0(df$Precursor.Id)
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ0\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ-K-0\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ-n-0\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ4\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ-K-4\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ-n-4\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ8\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ-K-8\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ-n-8\\)")
  return(df)
}

pD_rmMixSpec <- function(df){
  df$HY<-F
  df$HY[grepl("HUMAN",df$Protein.Names) & grepl("YEAST",df$Protein.Names) & grepl("ECOLI",df$Protein.Names)] <-T
  df$HY[grepl("HUMAN",df$Protein.Names) & grepl("ECOLI",df$Protein.Names)] <-T
  df$HY[grepl("HUMAN",df$Protein.Names) & grepl("YEAST",df$Protein.Names)] <-T
  df$HY[grepl("ECOLI",df$Protein.Names) & grepl("YEAST",df$Protein.Names)] <-T
  table(df$HY)
  df <- df[grepl("FALSE", df$HY),] %>% dplyr::select(-"HY")
  
  return(df)
}

pD_species <- function(df){
  
  df <- df%>% dplyr::mutate("species" = ifelse(grepl("ECOLI", Protein.Names), "E. coli",
                                        ifelse(grepl("HUMAN", Protein.Names), "H. sapiens",
                                               ifelse(grepl("YEAST", Protein.Names), "S. cerevisiae", "remove"))))
  
  df <- df[!grepl("remove", df$species),]
  
  return(df)
}

pD_Cterm <- function(df){
  df$Cterm <- str_sub(df$Stripped.Sequence, -1, -1)
  df <- df %>% dplyr::mutate("Cterm" = ifelse(grepl("K", Cterm), "K", "R"))
}

pD_channel_LF <- function(df){
  
  df <- df %>% dplyr::mutate("channel_name" = ifelse(grepl("Melanoma", Celltype), "mTRAQ0",
                                              ifelse(grepl("PDAC", Celltype), "mTRAQ4", "mTRAQ8")))
  
  return(df)
}

#Quant functions
pD_isotopicCO <- function(df){
  #df <- b7_1
  
  cnames <- colnames(df)
  df <- pD_seqcharge(df)
  df$seqcharge_run <- paste0(df$seqcharge, df$Run)
  
  mTRAQ0 <- list(C=7, H=12, N=2, O=1)
  mTRAQ4 <- list(C=4, H=12, N=1, O=1)
  mTRAQ8 <- list(C=1, H=12, N=0, O=1)
  
  
  # count number of each mTRAQ label type per precursor:
  df$d0 <- str_count(df$Modified.Sequence, "mTRAQ0|mTRAQ-K-0|mTRAQ-n-0")
  df$d4 <- str_count(df$Modified.Sequence, "mTRAQ4|mTRAQ-K-4|mTRAQ-n-4")
  df$d8 <- str_count(df$Modified.Sequence, "mTRAQ8|mTRAQ-K-8|mTRAQ-n-8")
  
  uni_prec <- df %>% dplyr::distinct(Modified.Sequence, .keep_all=T)
  
  deamid <- list(N=-1, H=-1, O=1)
  ox <- list(O=1)
  carba <- list(C=2, H=3, N=1, O=1)
  
  precs <- vector(mode = "list", length = nrow(uni_prec))
  # get chemical formula for each peptide
  for(i in 1:nrow(uni_prec)){
    tempseq <- paste0(uni_prec$Stripped.Sequence[i])
    modseq <- paste0(uni_prec$Modified.Sequence[i])
    numLab_d0 <- uni_prec$d0[i]
    numLab_d4 <- uni_prec$d4[i]
    numLab_d8 <- uni_prec$d8[i]
    
    ## other modifications
    d <- str_count(uni_prec$Precursor.Id[i], "UniMod:7")
    o <- str_count(uni_prec$Precursor.Id[i], "UniMod:35")
    c <- str_count(uni_prec$Precursor.Id[i], "UniMod:4")
    
    ch <- as.numeric(paste0(uni_prec$Precursor.Charge[i]))
    el_temp <- ConvertPeptide(tempseq)
    el_temp$C <- el_temp$C+(mTRAQ0$C*numLab_d0)+(mTRAQ4$C*numLab_d4)+(mTRAQ8$C*numLab_d8)+(carba$C*c)
    el_temp$H <- el_temp$H+(mTRAQ0$H*numLab_d0)+(mTRAQ4$H*numLab_d4)+(mTRAQ8$H*numLab_d8)+(carba$H*c)
    el_temp$N <- el_temp$N+(mTRAQ0$N*numLab_d0)+(mTRAQ4$N*numLab_d4)+(mTRAQ8$N*numLab_d8)+(carba$N*c)
    el_temp$O <- el_temp$O+(mTRAQ0$O*numLab_d0)+(mTRAQ4$O*numLab_d4)+(mTRAQ8$O*numLab_d8)+(carba$O*c)
    el_temp_df <- data.frame(el_temp)
    rownames(el_temp_df) <- modseq
    precs[[i]] <- el_temp_df
  }
  
  iso_env_list <- vector(mode = "list", length = length(precs))
  # calculate distribution of isotopic envelope
  for(i in 1:length(precs)){
    modseq <- rownames(data.frame(precs[i]))
    temp <- as.numeric(data.frame(precs[i]))
    iso <- iso_distr(temp)
    iso_df <- data.frame(t(iso))
    rownames(iso_df) <- modseq
    iso_env_list[[i]] <- iso
  }
  
  #unlist into a dataframe
  indx <- sapply(iso_env_list, length)
  iso_env_df <- as.data.frame(do.call(rbind,lapply(iso_env_list, `length<-`,
                                                   max(indx))))
  colnames(iso_env_df) <- c(paste0("X", seq(1,ncol(iso_env_df), by =1)))
  rnames <- uni_prec$Modified.Sequence
  rownames(iso_env_df) <- rnames
  iso_env_df[is.na(iso_env_df)] <- 0
  all_carryOver <- iso_env_df %>% rowSums()
  iso_env_df <- iso_env_df/all_carryOver  #get fraction of each peak 
  iso_env_df$modseq <- rnames
  
  # Now join with the quantitative data, and correct for isotopic carryover
  dat <- df %>% left_join(iso_env_df, by =c("Modified.Sequence" = "modseq"))
  dat$num_labs <- dat$d0+dat$d4+dat$d8
  dat <- pD_channel(dat)
  dat <- dat %>% mutate("y" = (X5+X6)/(X1+X2+X5+X6))
  d0_dat <- dat[grepl("mTRAQ0|mTRAQ-K-0|mTRAQ-n-0", dat$Modified.Sequence),] %>%
    dplyr::select("seqcharge_run", "num_labs", "Ms1.Area", "y") %>% dplyr::rename("d0"="Ms1.Area", "y1" = "y")
  d4_dat <- dat[grepl("mTRAQ4|mTRAQ-K-4|mTRAQ-n-4", dat$Modified.Sequence),] %>%
    dplyr::select("seqcharge_run", "Ms1.Area", "y") %>% dplyr::rename("d4"="Ms1.Area", "y2" = "y")
  d8_dat <- dat[grepl("mTRAQ8|mTRAQ-K-8|mTRAQ-n-8", dat$Modified.Sequence),] %>%
    dplyr::select("seqcharge_run", "Ms1.Area", "y") %>% dplyr::rename("d8"="Ms1.Area", "y3" = "y")
  all_dat <- d0_dat %>% full_join(d4_dat, by =c("seqcharge_run" = "seqcharge_run")) %>%
    full_join(d8_dat, by =c("seqcharge_run" = "seqcharge_run"))
  all_dat$d0_iso <- all_dat$d0+all_dat$d0*all_dat$y1
  all_dat$d4_iso <- all_dat$d4 - all_dat$d0*all_dat$y1 + all_dat$d4*all_dat$y2
  all_dat$d8_iso <- all_dat$d8 - all_dat$d4*all_dat$y2 + all_dat$d8*all_dat$y3
  #only do the isotopic correction for precursors with < 2 labels.. otherwise, the difference in mass is enough to avoid impactful isotopic carryover (such is the case with K-containing peptides).
  all_dat <- all_dat %>% dplyr::mutate(d0_iso = ifelse(num_labs<2, d0_iso, d0)) 
  all_dat <- all_dat %>% dplyr::mutate(d4_iso = ifelse(num_labs<2, d4_iso, d4)) 
  all_dat <- all_dat %>% dplyr::mutate(d8_iso = ifelse(num_labs<2, d8_iso, d8)) 
  all_dat$d0_iso[all_dat$d0_iso<0]<-0 #replace negatives with 0
  all_dat$d4_iso[all_dat$d4_iso<0]<-0 #replace negatives with 0
  all_dat$d8_iso[all_dat$d8_iso<0]<-0 #replace negatives with 0
  all_dat_d0 <- all_dat %>% dplyr::select("seqcharge_run", "d0_iso") %>% dplyr::rename("Ms1.Area_iso" = "d0_iso") %>% na.omit()
  all_dat_d4 <- all_dat %>% dplyr::select("seqcharge_run", "d4_iso")%>% dplyr::rename("Ms1.Area_iso" = "d4_iso") %>% na.omit()
  all_dat_d8 <- all_dat %>% dplyr::select("seqcharge_run", "d8_iso")%>% dplyr::rename("Ms1.Area_iso" = "d8_iso") %>% na.omit()
  d0_dat <- dat[grepl("mTRAQ0|mTRAQ-K-0|mTRAQ-n-0", dat$Modified.Sequence),]
  d4_dat <- dat[grepl("mTRAQ4|mTRAQ-K-4|mTRAQ-n-4", dat$Modified.Sequence),]
  d8_dat <- dat[grepl("mTRAQ8|mTRAQ-K-8|mTRAQ-n-8", dat$Modified.Sequence),]
  
  d0_dat <- d0_dat %>% left_join(all_dat_d0, by =c("seqcharge_run" = "seqcharge_run"))
  d4_dat <- d4_dat %>% left_join(all_dat_d4, by =c("seqcharge_run" = "seqcharge_run"))
  d8_dat <- d8_dat %>% left_join(all_dat_d8, by =c("seqcharge_run" = "seqcharge_run"))
  dat <- rbind(d0_dat, d4_dat, d8_dat)
  dat$Ms1.Area_iso[is.na(dat$Ms1.Area_iso)] <-0
  dat <- dat %>% dplyr::select(all_of(cnames), "Ms1.Area_iso")
  #b7_2 <- dat
  return(dat)
}

pD_isoCO_Ms1Ext <- function(df, Ms1Ext_df){
  #df <- QE1
  # Ms1Ext_df <- QE_Ms1ext
  uni_runs <- unique(df$Run)
  df <- pD_seqcharge(df)
  df$seqcharge_run <- paste0(df$seqcharge, df$Run)
  Ms1Ext_df <- pD_seqcharge(Ms1Ext_df)
  Ms1Ext_df <- Ms1Ext_df %>% dplyr::select(-contains("Quality")) %>%
    dplyr::select(-contains("decoy")) %>% dplyr::select(-contains("QValue"))
  #Ms1Ext_df <- Ms1Ext_df[,!grepl("decoy", colnames(Ms1Ext_df))]
  #Ms1Ext_df <- Ms1Ext_df[,!grepl("QValue", colnames(Ms1Ext_df))] #remove qvalue columns
  #Ms1Ext_df <- Ms1Ext_df[,!grepl("Quality", colnames(Ms1Ext_df))] #remove qvalue columns
  final_all <- data.frame(matrix(ncol = 6, nrow = 0))
  colnames(final_all) <- c("Protein.Group", "Protein.Names", "Ms1.Extracted", "Precursor.Charge","seqcharge", "Run")
  for(i in 1:length(uni_runs)){
    temp <- Ms1Ext_df[,grepl(paste0(uni_runs[i],"|","Modified.Sequence","|","Stripped.Sequence","|","Protein.Group|seqcharge|Protein.Names|Precursor.Charge"), colnames(Ms1Ext_df))]
    temp <- temp %>% dplyr::select(-contains("Quality"))
    x <- substr(colnames(temp), nchar(colnames(temp)), nchar(colnames(temp)))
    #d0
    temp0 <- temp[,grepl("Precursor.Charge|Modified.Sequence|Stripped.Sequence|Protein.Group|seqcharge|Protein.Names|raw.0", colnames(temp))]
    colnames(temp0) <- c("Protein.Group","Protein.Names", "Stripped.Sequence","Modified.Sequence", "Precursor.Charge","Ms1.Extracted", "seqcharge")
    temp0$channel_name <- "mTRAQ0"
    temp0$Modified.Sequence <- str_replace_all(temp0$Modified.Sequence, "\\(mTRAQ\\)","\\(mTRAQ0\\)")
    temp0 <- na.omit(temp0)
    #d4
    temp4 <- temp[,grepl("Precursor.Charge|Modified.Sequence|Stripped.Sequence|Protein.Group|seqcharge|Protein.Names|raw.4", colnames(temp))]
    colnames(temp4) <- c("Protein.Group","Protein.Names","Stripped.Sequence","Modified.Sequence","Precursor.Charge", "Ms1.Extracted", "seqcharge")
    temp4$channel_name <- "mTRAQ4"
    temp4$Modified.Sequence <- str_replace_all(temp4$Modified.Sequence, "\\(mTRAQ\\)","\\(mTRAQ4\\)")
    temp4 <- na.omit(temp4)
    #d8
    temp8 <- temp[,grepl("Precursor.Charge|Modified.Sequence|Stripped.Sequence|Protein.Group|seqcharge|Protein.Names|raw.8", colnames(temp))]
    colnames(temp8) <-  c("Protein.Group","Protein.Names","Stripped.Sequence","Modified.Sequence", "Precursor.Charge","Ms1.Extracted", "seqcharge")
    temp8$channel_name <- "mTRAQ8"
    temp8$Modified.Sequence <- str_replace_all(temp8$Modified.Sequence, "\\(mTRAQ\\)","\\(mTRAQ8\\)")
    temp8 <- na.omit(temp8)
    #combine
    all_temp <- rbind(temp0,temp4,temp8)
    all_temp$Run <- paste0(uni_runs[i])
    final_all <- rbind(final_all, all_temp)
  }
  final_all$seqcharge_run <- paste0(final_all$seqcharge, final_all$Run)
  cnames <- colnames(final_all)
  mTRAQ0 <- list(C=7, H=12, N=2, O=1)
  mTRAQ4 <- list(C=4, H=12, N=1, O=1)
  mTRAQ8 <- list(C=1, H=12, N=0, O=1)
  
  
  # count number of each mTRAQ label type per precursor:
  final_all$d0 <- str_count(final_all$Modified.Sequence, "mTRAQ0|mTRAQ-K-0|mTRAQ-n-0")
  final_all$d4 <- str_count(final_all$Modified.Sequence, "mTRAQ4|mTRAQ-K-4|mTRAQ-n-4")
  final_all$d8 <- str_count(final_all$Modified.Sequence, "mTRAQ8|mTRAQ-K-8|mTRAQ-n-8")
  
  uni_prec <- final_all %>% dplyr::distinct(Modified.Sequence, .keep_all=T)
  
  deamid <- list(N=-1, H=-1, O=1)
  ox <- list(O=1)
  carba <- list(C=2, H=3, N=1, O=1)
  
  precs <- vector(mode = "list", length = nrow(uni_prec))
  # get chemical formula for each peptide
  for(i in 1:nrow(uni_prec)){
    tempseq <- paste0(uni_prec$Stripped.Sequence[i])
    modseq <- paste0(uni_prec$Modified.Sequence[i])
    numLab_d0 <- uni_prec$d0[i]
    numLab_d4 <- uni_prec$d4[i]
    numLab_d8 <- uni_prec$d8[i]
    
    ## other modifications
    d <- str_count(uni_prec$Modified.Sequence[i], "UniMod:7")
    o <- str_count(uni_prec$Modified.Sequence[i], "UniMod:35")
    c <- str_count(uni_prec$Modified.Sequence[i], "UniMod:4")
    
    # ch <- as.numeric(paste0(uni_prec$Precursor.Charge[i]))
    el_temp <- ConvertPeptide(tempseq)
    el_temp$C <- el_temp$C+(mTRAQ0$C*numLab_d0)+(mTRAQ4$C*numLab_d4)+(mTRAQ8$C*numLab_d8)+(carba$C*c)
    el_temp$H <- el_temp$H+(mTRAQ0$H*numLab_d0)+(mTRAQ4$H*numLab_d4)+(mTRAQ8$H*numLab_d8)+(carba$H*c)
    el_temp$N <- el_temp$N+(mTRAQ0$N*numLab_d0)+(mTRAQ4$N*numLab_d4)+(mTRAQ8$N*numLab_d8)+(carba$N*c)
    el_temp$O <- el_temp$O+(mTRAQ0$O*numLab_d0)+(mTRAQ4$O*numLab_d4)+(mTRAQ8$O*numLab_d8)+(carba$O*c)
    el_temp_df <- data.frame(el_temp)
    rownames(el_temp_df) <- modseq
    precs[[i]] <- el_temp_df
  }
  
  iso_env_list <- vector(mode = "list", length = length(precs))
  # calculate distribution of isotopic envelope
  for(i in 1:length(precs)){
    modseq <- rownames(data.frame(precs[i]))
    temp <- as.numeric(data.frame(precs[i]))
    iso <- iso_distr(temp)
    iso_df <- data.frame(t(iso))
    rownames(iso_df) <- modseq
    iso_env_list[[i]] <- iso
  }
  
  #unlist into a dataframe
  indx <- sapply(iso_env_list, length)
  iso_env_df <- as.data.frame(do.call(rbind,lapply(iso_env_list, `length<-`,
                                                   max(indx))))
  colnames(iso_env_df) <- c(paste0("X", seq(1,ncol(iso_env_df), by =1)))
  rnames <- uni_prec$Modified.Sequence
  rownames(iso_env_df) <- rnames
  iso_env_df[is.na(iso_env_df)] <- 0
  all_carryOver <- iso_env_df %>% rowSums()
  iso_env_df <- iso_env_df/all_carryOver  #get fraction of each peak 
  iso_env_df$modseq <- rnames
  
  # Now join with the quantitative data, and correct for isotopic carryover
  dat <- final_all %>% left_join(iso_env_df, by =c("Modified.Sequence" = "modseq"))
  dat$num_labs <- dat$d0+dat$d4+dat$d8
  dat <- pD_channel(dat)
  dat <- dat %>% dplyr::mutate("y" = (X5+X6)/(X1+X2+X5+X6))
  d0_dat <- dat[grepl("mTRAQ0|mTRAQ-K-0|mTRAQ-n-0", dat$Modified.Sequence),] %>%
    dplyr::select("seqcharge_run", "num_labs", "Ms1.Extracted", "y") %>% dplyr::rename("d0"="Ms1.Extracted", "y1" = "y")
  d4_dat <- dat[grepl("mTRAQ4|mTRAQ-K-4|mTRAQ-n-4", dat$Modified.Sequence),] %>%
    dplyr::select("seqcharge_run", "Ms1.Extracted", "y") %>% dplyr::rename("d4"="Ms1.Extracted", "y2" = "y")
  d8_dat <- dat[grepl("mTRAQ8|mTRAQ-K-8|mTRAQ-n-8", dat$Modified.Sequence),] %>%
    dplyr::select("seqcharge_run", "Ms1.Extracted", "y") %>% dplyr::rename("d8"="Ms1.Extracted", "y3" = "y")
  all_dat <- d0_dat %>% full_join(d4_dat, by =c("seqcharge_run" = "seqcharge_run")) %>%
    full_join(d8_dat, by =c("seqcharge_run" = "seqcharge_run"))
  all_dat$d0_iso <- all_dat$d0+all_dat$d0*all_dat$y1
  all_dat$d4_iso <- all_dat$d4 - all_dat$d0*all_dat$y1 + all_dat$d4*all_dat$y2
  all_dat$d8_iso <- all_dat$d8 - all_dat$d4*all_dat$y2 + all_dat$d8*all_dat$y3
  #only do the isotopic correction for precursors with < 2 labels.. otherwise, the difference in mass is enough to avoid impactful isotopic carryover (such is the case with K-containing peptides).
  all_dat <- all_dat %>% dplyr::mutate(d0_iso = ifelse(num_labs<2, d0_iso, d0)) 
  all_dat <- all_dat %>% dplyr::mutate(d4_iso = ifelse(num_labs<2, d4_iso, d4)) 
  all_dat <- all_dat %>% dplyr::mutate(d8_iso = ifelse(num_labs<2, d8_iso, d8)) 
  all_dat$d0_iso[all_dat$d0_iso<0]<-0 #replace negatives with 0
  all_dat$d4_iso[all_dat$d4_iso<0]<-0 #replace negatives with 0
  all_dat$d8_iso[all_dat$d8_iso<0]<-0 #replace negatives with 0
  all_dat_d0 <- all_dat %>% dplyr::select("seqcharge_run", "d0_iso") %>% dplyr::rename("Ms1.Extracted_iso" = "d0_iso") %>% na.omit()
  all_dat_d4 <- all_dat %>% dplyr::select("seqcharge_run", "d4_iso")%>% dplyr::rename("Ms1.Extracted_iso" = "d4_iso") %>% na.omit()
  all_dat_d8 <- all_dat %>% dplyr::select("seqcharge_run", "d8_iso")%>% dplyr::rename("Ms1.Extracted_iso" = "d8_iso") %>% na.omit()
  d0_dat <- dat[grepl("mTRAQ0|mTRAQ-K-0|mTRAQ-n-0", dat$Modified.Sequence),]
  d4_dat <- dat[grepl("mTRAQ4|mTRAQ-K-4|mTRAQ-n-4", dat$Modified.Sequence),]
  d8_dat <- dat[grepl("mTRAQ8|mTRAQ-K-8|mTRAQ-n-8", dat$Modified.Sequence),]
  
  d0_dat <- d0_dat %>% left_join(all_dat_d0, by =c("seqcharge_run" = "seqcharge_run"))
  d4_dat <- d4_dat %>% left_join(all_dat_d4, by =c("seqcharge_run" = "seqcharge_run"))
  d8_dat <- d8_dat %>% left_join(all_dat_d8, by =c("seqcharge_run" = "seqcharge_run"))
  dat <- rbind(d0_dat, d4_dat, d8_dat)
  dat$Ms1.Extracted_iso[is.na(dat$Ms1.Extracted_iso)] <-0
  dat <- dat %>% dplyr::select(all_of(cnames), "Ms1.Extracted_iso")
  
  return(dat)
}

pD_isoCO_Ms1Ext_tims <- function(df, Ms1Ext_df){
  #df <- QE1
  # Ms1Ext_df <- QE_Ms1ext
  uni_runs <- unique(df$Run)
  df <- pD_seqcharge(df)
  df$seqcharge_run <- paste0(df$seqcharge, df$Run)
  Ms1Ext_df <- pD_seqcharge(Ms1Ext_df)
  Ms1Ext_df <- Ms1Ext_df %>% dplyr::select(-contains("Quality")) %>%
    dplyr::select(-contains("decoy")) %>% dplyr::select(-contains("QValue"))
  #Ms1Ext_df <- Ms1Ext_df[,!grepl("decoy", colnames(Ms1Ext_df))]
  #Ms1Ext_df <- Ms1Ext_df[,!grepl("QValue", colnames(Ms1Ext_df))] #remove qvalue columns
  #Ms1Ext_df <- Ms1Ext_df[,!grepl("Quality", colnames(Ms1Ext_df))] #remove qvalue columns
  final_all <- data.frame(matrix(ncol = 6, nrow = 0))
  colnames(final_all) <- c("Protein.Group", "Protein.Names", "Ms1.Extracted", "Precursor.Charge","seqcharge", "Run")
  for(i in 1:length(uni_runs)){
    temp <- Ms1Ext_df[,grepl(paste0(uni_runs[i],"|","Modified.Sequence","|","Stripped.Sequence","|","Protein.Group|seqcharge|Protein.Names|Precursor.Charge"), colnames(Ms1Ext_df))]
    temp <- temp %>% dplyr::select(-contains("Quality"))
    x <- substr(colnames(temp), nchar(colnames(temp)), nchar(colnames(temp)))
    #d0
    temp0 <- temp[,grepl("Precursor.Charge|Modified.Sequence|Stripped.Sequence|Protein.Group|seqcharge|Protein.Names|d.0", colnames(temp))]
    colnames(temp0) <- c("Protein.Group","Protein.Names", "Stripped.Sequence","Modified.Sequence", "Precursor.Charge","Ms1.Extracted", "seqcharge")
    temp0$channel_name <- "mTRAQ0"
    temp0$Modified.Sequence <- str_replace_all(temp0$Modified.Sequence, "\\(mTRAQ\\)","\\(mTRAQ0\\)")
    temp0 <- na.omit(temp0)
    #d4
    temp4 <- temp[,grepl("Precursor.Charge|Modified.Sequence|Stripped.Sequence|Protein.Group|seqcharge|Protein.Names|d.4", colnames(temp))]
    colnames(temp4) <- c("Protein.Group","Protein.Names","Stripped.Sequence","Modified.Sequence","Precursor.Charge", "Ms1.Extracted", "seqcharge")
    temp4$channel_name <- "mTRAQ4"
    temp4$Modified.Sequence <- str_replace_all(temp4$Modified.Sequence, "\\(mTRAQ\\)","\\(mTRAQ4\\)")
    temp4 <- na.omit(temp4)
    #d8
    temp8 <- temp[,grepl("Precursor.Charge|Modified.Sequence|Stripped.Sequence|Protein.Group|seqcharge|Protein.Names|d.8", colnames(temp))]
    colnames(temp8) <-  c("Protein.Group","Protein.Names","Stripped.Sequence","Modified.Sequence", "Precursor.Charge","Ms1.Extracted", "seqcharge")
    temp8$channel_name <- "mTRAQ8"
    temp8$Modified.Sequence <- str_replace_all(temp8$Modified.Sequence, "\\(mTRAQ\\)","\\(mTRAQ8\\)")
    temp8 <- na.omit(temp8)
    #combine
    all_temp <- rbind(temp0,temp4,temp8)
    all_temp$Run <- paste0(uni_runs[i])
    final_all <- rbind(final_all, all_temp)
  }
  final_all$seqcharge_run <- paste0(final_all$seqcharge, final_all$Run)
  cnames <- colnames(final_all)
  mTRAQ0 <- list(C=7, H=12, N=2, O=1)
  mTRAQ4 <- list(C=4, H=12, N=1, O=1)
  mTRAQ8 <- list(C=1, H=12, N=0, O=1)
  
  
  # count number of each mTRAQ label type per precursor:
  final_all$d0 <- str_count(final_all$Modified.Sequence, "mTRAQ0|mTRAQ-K-0|mTRAQ-n-0")
  final_all$d4 <- str_count(final_all$Modified.Sequence, "mTRAQ4|mTRAQ-K-4|mTRAQ-n-4")
  final_all$d8 <- str_count(final_all$Modified.Sequence, "mTRAQ8|mTRAQ-K-8|mTRAQ-n-8")
  
  uni_prec <- final_all %>% distinct(Modified.Sequence, .keep_all=T)
  
  deamid <- list(N=-1, H=-1, O=1)
  ox <- list(O=1)
  carba <- list(C=2, H=3, N=1, O=1)
  
  precs <- vector(mode = "list", length = nrow(uni_prec))
  # get chemical formula for each peptide
  for(i in 1:nrow(uni_prec)){
    tempseq <- paste0(uni_prec$Stripped.Sequence[i])
    modseq <- paste0(uni_prec$Modified.Sequence[i])
    numLab_d0 <- uni_prec$d0[i]
    numLab_d4 <- uni_prec$d4[i]
    numLab_d8 <- uni_prec$d8[i]
    
    ## other modifications
    d <- str_count(uni_prec$Modified.Sequence[i], "UniMod:7")
    o <- str_count(uni_prec$Modified.Sequence[i], "UniMod:35")
    c <- str_count(uni_prec$Modified.Sequence[i], "UniMod:4")
    
    # ch <- as.numeric(paste0(uni_prec$Precursor.Charge[i]))
    el_temp <- ConvertPeptide(tempseq)
    el_temp$C <- el_temp$C+(mTRAQ0$C*numLab_d0)+(mTRAQ4$C*numLab_d4)+(mTRAQ8$C*numLab_d8)+(carba$C*c)
    el_temp$H <- el_temp$H+(mTRAQ0$H*numLab_d0)+(mTRAQ4$H*numLab_d4)+(mTRAQ8$H*numLab_d8)+(carba$H*c)
    el_temp$N <- el_temp$N+(mTRAQ0$N*numLab_d0)+(mTRAQ4$N*numLab_d4)+(mTRAQ8$N*numLab_d8)+(carba$N*c)
    el_temp$O <- el_temp$O+(mTRAQ0$O*numLab_d0)+(mTRAQ4$O*numLab_d4)+(mTRAQ8$O*numLab_d8)+(carba$O*c)
    el_temp_df <- data.frame(el_temp)
    rownames(el_temp_df) <- modseq
    precs[[i]] <- el_temp_df
  }
  
  iso_env_list <- vector(mode = "list", length = length(precs))
  # calculate distribution of isotopic envelope
  for(i in 1:length(precs)){
    modseq <- rownames(data.frame(precs[i]))
    temp <- as.numeric(data.frame(precs[i]))
    iso <- iso_distr(temp)
    iso_df <- data.frame(t(iso))
    rownames(iso_df) <- modseq
    iso_env_list[[i]] <- iso
  }
  
  #unlist into a dataframe
  indx <- sapply(iso_env_list, length)
  iso_env_df <- as.data.frame(do.call(rbind,lapply(iso_env_list, `length<-`,
                                                   max(indx))))
  colnames(iso_env_df) <- c(paste0("X", seq(1,ncol(iso_env_df), by =1)))
  rnames <- uni_prec$Modified.Sequence
  rownames(iso_env_df) <- rnames
  iso_env_df[is.na(iso_env_df)] <- 0
  all_carryOver <- iso_env_df %>% rowSums()
  iso_env_df <- iso_env_df/all_carryOver  #get fraction of each peak 
  iso_env_df$modseq <- rnames
  
  # Now join with the quantitative data, and correct for isotopic carryover
  dat <- final_all %>% left_join(iso_env_df, by =c("Modified.Sequence" = "modseq"))
  dat$num_labs <- dat$d0+dat$d4+dat$d8
  dat <- pD_channel(dat)
  dat <- dat %>% dplyr::mutate("y" = (X5+X6)/(X1+X2+X5+X6))
  d0_dat <- dat[grepl("mTRAQ0|mTRAQ-K-0|mTRAQ-n-0", dat$Modified.Sequence),] %>%
    dplyr::select("seqcharge_run", "num_labs", "Ms1.Extracted", "y") %>% dplyr::rename("d0"="Ms1.Extracted", "y1" = "y")
  d4_dat <- dat[grepl("mTRAQ4|mTRAQ-K-4|mTRAQ-n-4", dat$Modified.Sequence),] %>%
    dplyr::select("seqcharge_run", "Ms1.Extracted", "y") %>% dplyr::rename("d4"="Ms1.Extracted", "y2" = "y")
  d8_dat <- dat[grepl("mTRAQ8|mTRAQ-K-8|mTRAQ-n-8", dat$Modified.Sequence),] %>%
    dplyr::select("seqcharge_run", "Ms1.Extracted", "y") %>% dplyr::rename("d8"="Ms1.Extracted", "y3" = "y")
  all_dat <- d0_dat %>% full_join(d4_dat, by =c("seqcharge_run" = "seqcharge_run")) %>%
    full_join(d8_dat, by =c("seqcharge_run" = "seqcharge_run"))
  all_dat$d0_iso <- all_dat$d0+all_dat$d0*all_dat$y1
  all_dat$d4_iso <- all_dat$d4 - all_dat$d0*all_dat$y1 + all_dat$d4*all_dat$y2
  all_dat$d8_iso <- all_dat$d8 - all_dat$d4*all_dat$y2 + all_dat$d8*all_dat$y3
  #only do the isotopic correction for precursors with < 2 labels.. otherwise, the difference in mass is enough to avoid impactful isotopic carryover (such is the case with K-containing peptides).
  all_dat <- all_dat %>% dplyr::mutate(d0_iso = ifelse(num_labs<2, d0_iso, d0)) 
  all_dat <- all_dat %>% dplyr::mutate(d4_iso = ifelse(num_labs<2, d4_iso, d4)) 
  all_dat <- all_dat %>% dplyr::mutate(d8_iso = ifelse(num_labs<2, d8_iso, d8)) 
  all_dat$d0_iso[all_dat$d0_iso<0]<-0 #replace negatives with 0
  all_dat$d4_iso[all_dat$d4_iso<0]<-0 #replace negatives with 0
  all_dat$d8_iso[all_dat$d8_iso<0]<-0 #replace negatives with 0
  all_dat_d0 <- all_dat %>% dplyr::select("seqcharge_run", "d0_iso") %>% dplyr::rename("Ms1.Extracted_iso" = "d0_iso") %>% na.omit()
  all_dat_d4 <- all_dat %>% dplyr::select("seqcharge_run", "d4_iso")%>% dplyr::rename("Ms1.Extracted_iso" = "d4_iso") %>% na.omit()
  all_dat_d8 <- all_dat %>% dplyr::select("seqcharge_run", "d8_iso")%>% dplyr::rename("Ms1.Extracted_iso" = "d8_iso") %>% na.omit()
  d0_dat <- dat[grepl("mTRAQ0|mTRAQ-K-0|mTRAQ-n-0", dat$Modified.Sequence),]
  d4_dat <- dat[grepl("mTRAQ4|mTRAQ-K-4|mTRAQ-n-4", dat$Modified.Sequence),]
  d8_dat <- dat[grepl("mTRAQ8|mTRAQ-K-8|mTRAQ-n-8", dat$Modified.Sequence),]
  
  d0_dat <- d0_dat %>% left_join(all_dat_d0, by =c("seqcharge_run" = "seqcharge_run"))
  d4_dat <- d4_dat %>% left_join(all_dat_d4, by =c("seqcharge_run" = "seqcharge_run"))
  d8_dat <- d8_dat %>% left_join(all_dat_d8, by =c("seqcharge_run" = "seqcharge_run"))
  dat <- rbind(d0_dat, d4_dat, d8_dat)
  dat$Ms1.Extracted_iso[is.na(dat$Ms1.Extracted_iso)] <-0
  dat <- dat %>% dplyr::select(all_of(cnames), "Ms1.Extracted_iso")
  
  return(dat)
}

pD_diann_quant <- function(df, group.header="Protein.Group", id.header = "Unlabelled.Precursor", quantity.header = "Ms1.Area"){
  protein.groups <- diann_maxlfq(df[df$Lib.PG.Q.Value <= 0.01,], group.header=paste0(group.header), id.header = paste0(id.header), quantity.header = paste0(quantity.header))
  protein.groups_n <- protein.groups                                      
  PGs <- rownames(protein.groups)
  protein.groups_n <- cbind(protein.groups_n,PGs)
  df_names <- df %>% dplyr::select("Protein.Group", "Protein.Names") %>% dplyr::distinct(Protein.Group, .keep_all=T)
  
  protein.groups_n <- data.frame(protein.groups_n)
  protein.groups_n <- protein.groups_n %>% left_join(df_names, by =c("PGs"="Protein.Group"))
  
  protein.groups_n$HY<-F
  protein.groups_n$HY[grepl("HUMAN",protein.groups_n$Protein.Names) & grepl("YEAST",protein.groups_n$Protein.Names) & grepl("ECOLI",protein.groups_n$Protein.Names)] <-T
  protein.groups_n$HY[grepl("HUMAN",protein.groups_n$Protein.Names) & grepl("ECOLI",protein.groups_n$Protein.Names)] <-T
  protein.groups_n$HY[grepl("HUMAN",protein.groups_n$Protein.Names) & grepl("YEAST",protein.groups_n$Protein.Names)] <-T
  protein.groups_n$HY[grepl("ECOLI",protein.groups_n$Protein.Names) & grepl("YEAST",protein.groups_n$Protein.Names)] <-T
  table(protein.groups_n$HY)
  protein.groups_n <- protein.groups_n[grepl("FALSE", protein.groups_n$HY),]
  
  protein.groups_n <- protein.groups_n%>% mutate("species" = ifelse(grepl("ECOLI", Protein.Names), "E. coli",
                                                                    ifelse(grepl("HUMAN", Protein.Names), "H. sapiens",
                                                                           ifelse(grepl("YEAST", Protein.Names), "S. cerevisiae", "remove"))))
  
  protein.groups_n <- protein.groups_n[!grepl("remove", protein.groups_n$species),]
  return(protein.groups_n)
}

pD_protQuant <- function(df, runHeader =  "Run", quantHeader = "Ms1.Area", 
                         groupHeader = "Protein.Group", sampleHeader = "channel_name",
                         reqIntPrecs_withinRun = FALSE, reqIntPrecs_wholedataset= FALSE){
  #df <- b7_1
  df[[quantHeader]] <- as.numeric(df[[quantHeader]])
  df[[groupHeader]] <- as.character(df[[groupHeader]])
  df[[sampleHeader]] <- as.character(df[[sampleHeader]])
  df[[runHeader]] <- as.character(df[[runHeader]])
  df$Run_sample <- paste0(df[[runHeader]],"_",df[[sampleHeader]])
  if(reqIntPrecs_withinRun==TRUE){
    df <- pD_seqcharge(df)
    df <- df %>% dplyr::group_by(seqcharge, Run) %>% dplyr::add_count() %>% dplyr::filter(n==3) %>% dplyr::ungroup()}
  if(reqIntPrecs_wholedataset==TRUE){
    df <- pD_seqcharge(df)
    df <- df %>% dplyr::group_by(seqcharge) %>% dplyr::add_count() %>% dplyr::filter(n==3*length(unique(.data[[runHeader]]))) %>% dplyr::ungroup()}
  df <- df %>% dplyr::group_by(.data[[groupHeader]], Run_sample) %>% 
    dplyr::mutate(Prot_sum  = sum(.data[[quantHeader]])) %>% dplyr::ungroup() %>%
    dplyr::select("Run_sample", "Protein.Names", paste0(groupHeader), "Prot_sum") %>% 
    dplyr::distinct(.data[[groupHeader]], Run_sample, .keep_all=T)
  df <- reshape2::dcast(df, df[[groupHeader]]+Protein.Names~Run_sample, value.var = "Prot_sum")
  colnames(df)[1] <- paste0(groupHeader)
  
  return(df)
  
  
}

#Plotting functions

pD_prepPlot <- function(prot_df, Run = "Run name"){
  #Run <- "wJD804"
  #prot_df <- QE7_2
  prot_df <- data.frame(prot_df)
  df <- pD_species(prot_df)
  pDIAa <- paste0(Run,"_mTRAQ0")
  pDIAb <- paste0(Run,"_mTRAQ4")
  pDIAc <- paste0(Run,"_mTRAQ8")
  
  df1 <- df
  df1$d0_d4 <- as.numeric(df1[[pDIAa]])/as.numeric(df1[[pDIAb]])
  df1 <- df1[!is.na(df1$d0_d4), ]
  df1 <- df1[!is.infinite(df1$d0_d4), ]
  df1 <- df1[which(df1$d0_d4>0), ]
  med_human <- df1[grepl("H. sapiens", df1$species),]
  med_human <- median(med_human$d0_d4, na.rm=T)
  scalar <- 1/med_human
  df1$d0_d4 <- df1$d0_d4*scalar   #normalize human median on 1:1 ratio
  df1 <- df1 %>% dplyr::select("Protein.Group", "species", all_of(pDIAb), "d0_d4")
  colnames(df1)[3] <- "sample2_int"
  df1$sample2_int <- as.numeric(df1$sample2_int)
  
  df2 <- df
  df2$d0_d8 <- as.numeric(df2[[pDIAa]])/as.numeric(df2[[pDIAc]])
  df2 <- df2[!is.na(df2$d0_d8), ]
  df2 <- df2[!is.infinite(df2$d0_d8), ]
  df2 <- df2[which(df2$d0_d8>0), ]
  med_human <- df2[grepl("H. sapiens", df2$species),]
  med_human <- median(med_human$d0_d8, na.rm=T)
  scalar <- 1/med_human
  df2$d0_d8 <- df2$d0_d8*scalar   #normalize human median on 1:1 ratio
  df2 <- df2 %>% dplyr::select("Protein.Group", "species", all_of(pDIAc), "d0_d8")
  colnames(df2)[3] <- "sample2_int"
  df2$sample2_int <- as.numeric(df2$sample2_int)
  
  df3 <- df
  df3$d4_d8 <- as.numeric(df3[[pDIAb]])/as.numeric(df3[[pDIAc]])
  df3 <- df3[!is.na(df3$d4_d8), ]
  df3 <- df3[!is.infinite(df3$d4_d8), ]
  df3 <- df3[which(df3$d4_d8>0), ]
  med_human <- df3[grepl("H. sapiens", df3$species),]
  med_human <- median(med_human$d4_d8, na.rm=T)
  scalar <- 1/med_human
  df3$d4_d8 <- df3$d4_d8*scalar   #normalize human median on 1:1 ratio
  df3 <- df3 %>% dplyr::select("Protein.Group", "species", all_of(pDIAc), "d4_d8")
  colnames(df3)[3] <- "sample2_int"
  df3$sample2_int <- as.numeric(df3$sample2_int)
  mylist <- list(df1,df2,df3)
  return(mylist)
}

pD_plotLFQBench <- function(pD_list, LF_list, pD_list_name = "plexDIA", LF_list_name = "LF-DIA"){
  
  
  #pD_list <- pD_list
  #LF_list <- LF_list
  #intersect LF and multiDIA
  protein.groups_n_omit <- pD_list[[1]]
  protein.groups_n_omit_LF_a <- LF_list[[1]]
  protein.groups_n_omit_noint <- pD_list[[1]]
  protein.groups_n_omit_LF_noint <- LF_list[[1]]
  
  protein.groups_n_omit_a <- protein.groups_n_omit[protein.groups_n_omit$Protein.Group%in%protein.groups_n_omit_LF_a$Protein.Group,]
  protein.groups_n_omit_LF_a <- protein.groups_n_omit_LF_a[protein.groups_n_omit_LF_a$Protein.Group%in%protein.groups_n_omit$Protein.Group,]
  
  
  LF_mDIA_d0d4 <- rbind(protein.groups_n_omit_a, protein.groups_n_omit_LF_a)
  
  ################################# plotting #########################
  
  x_d0d4 <- log2(LF_mDIA_d0d4$sample2_int) %>% sort()
  x_d0d4_val <- as.numeric(quantile(x_d0d4,probs=c(.0025,.9975)))
  y_d0d4 <- log2(LF_mDIA_d0d4$d0_d4) %>% sort()
  y_d0d4_val <- as.numeric(quantile(y_d0d4,probs=c(.0025,.9975)))
  ####################################
  #d0d4
  a2 <- ggplot(protein.groups_n_omit_LF_a, aes(x=log2(sample2_int), y=log2(d0_d4), color = species, alpha = species)) +
    geom_point() + 
    scale_color_manual(values = c("#619CFF","#00BA38", "#F8766D")) +
    geom_hline(yintercept=2, linetype="dashed", 
               color = "#619CFF", size=1)+
    geom_hline(yintercept=0, linetype="dashed", 
               color = "#00BA38", size=1)+
    geom_hline(yintercept=-1, linetype="dashed", 
               color = "#F8766D", size=1) +
    scale_alpha_manual(name = "category", values = c(0.3, 0.08, 0.22),guide = "none") +
    labs(x = expression(paste(Log["2"],", B")), y="Log2, delta0:delta4",title = paste0(LF_list_name))+
    #subtitle = paste0(nrow(LF_PGs_lim)," Protein groups"))+ 
    geom_smooth(method = "loess", method.args = list(family = "symmetric")) + 
    coord_cartesian(ylim = c(y_d0d4_val[1],y_d0d4_val[2]),
                    xlim = c(x_d0d4_val[1],x_d0d4_val[2])) +
    theme_classic()+
    theme(legend.position="none", 
          plot.title = element_text(hjust = 0.5, size=18),
          axis.ticks.y = element_blank(),
          plot.subtitle = element_text(size =12, hjust = 0.5),
          axis.title.y=element_blank(),
          axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=14),
          axis.text.y = element_blank())
  
  
  
  a4 <- ggplot(protein.groups_n_omit_LF_a, aes(x=species, y=log2(d0_d4), color = species))+
    geom_boxplot() +   
    scale_color_manual(values = c("#619CFF","#00BA38", "#F8766D")) +
    geom_hline(yintercept=2, linetype="dashed", 
               color = "#619CFF", size=1)+
    geom_hline(yintercept=0, linetype="dashed", 
               color = "#00BA38", size=1)+
    geom_hline(yintercept=-1, linetype="dashed", 
               color = "#F8766D", size=1) +
    scale_x_discrete(labels=c("H. sapiens" = "", "E. coli" = "",
                              "S. cerevisiae" = "")) +
    labs(x="",title = paste0(LF_list_name))+ #subtitle="")+
    theme(axis.ticks = element_blank(),
          axis.title.y=element_blank(),
          axis.text.y = element_blank(),
          panel.background=element_blank(), panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_text(size=32),
          plot.subtitle = element_text(size =12, hjust = 0.5),
          plot.title = element_text(size =18, hjust = 0.5)) +
    coord_cartesian(ylim = c(y_d0d4_val[1],y_d0d4_val[2]))
  ####################################
  #d0d4
  a1 <- ggplot(protein.groups_n_omit_a, aes(x=log2(sample2_int), y=log2(d0_d4), color = species, alpha = species)) +
    geom_point() + 
    scale_color_manual(values = c("#619CFF","#00BA38", "#F8766D")) +
    geom_hline(yintercept=2, linetype="dashed", 
               color = "#619CFF", size=1)+
    geom_hline(yintercept=0, linetype="dashed", 
               color = "#00BA38", size=1)+
    geom_hline(yintercept=-1, linetype="dashed", 
               color = "#F8766D", size=1) +
    scale_alpha_manual(name = "category", values = c(0.3, 0.08, 0.22),guide = "none") +
    labs(x = expression(paste(Log["2"],", B")), y=expression(paste(Log["2"],", A/B")),title = paste0(pD_list_name))+
    geom_smooth(method = "loess", method.args = list(family = "symmetric")) + 
    coord_cartesian(ylim = c(y_d0d4_val[1],y_d0d4_val[2]),
                    xlim = c(x_d0d4_val[1],x_d0d4_val[2])) +
    theme_classic()+
    theme(legend.position="none", 
          plot.title = element_text(size =18, hjust = 0.5),
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          plot.subtitle = element_text(size =12, hjust = 0.5),
          axis.title.x = element_text(size=14),
          axis.title.y = element_text(size=14),
          legend.text = element_text(face = "italic"))
  
  
  a3 <- ggplot(protein.groups_n_omit_a, aes(x=species, y=log2(d0_d4), color = species))+
    geom_boxplot() +   
    scale_color_manual(values = c("#619CFF","#00BA38", "#F8766D")) +
    geom_hline(yintercept=2, linetype="dashed", 
               color = "#619CFF", size=1)+
    geom_hline(yintercept=0, linetype="dashed", 
               color = "#00BA38", size=1)+
    geom_hline(yintercept=-1, linetype="dashed", 
               color = "#F8766D", size=1) +
    scale_x_discrete(labels=c("H. sapiens" = "", "E. coli" = "",
                              "S. cerevisiae" = "")) +
    labs(x="",title = paste0(pD_list_name), color = "Species:")+
    coord_cartesian(ylim = c(y_d0d4_val[1],y_d0d4_val[2])) +
    theme(axis.ticks = element_blank(),
          axis.title.y=element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          panel.background=element_blank(), panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.title.x = element_text(size=32),
          plot.background=element_blank(),
          plot.subtitle = element_text(size =12, hjust = 0.5),
          plot.title = element_text(size =18, hjust = 0.5),
          legend.position = "top",
          legend.text = element_text(size=14, face= "italic"),
          legend.title = element_text(size=14),
          legend.spacing.x = unit(1.0, 'line'))
  
  ################################## coverage of protein group ratios:
  joined <- protein.groups_n_omit[protein.groups_n_omit$Protein.Group%in%protein.groups_n_omit_LF_a$Protein.Group,]
  
  types <- c("plexDIA", "LF-DIA", "Intersected")
  numbers <- c(nrow(protein.groups_n_omit_noint), nrow(protein.groups_n_omit_LF_noint), nrow(joined))
  df_cover <- data.frame("Cond" = types, "IDs" = numbers)
  
  df_cover$Cond <- factor(df_cover$Cond, levels=c("plexDIA", "LF-DIA", "Intersected"))
  
  a_cover <- ggplot(df_cover, aes(x=as.factor(Cond), y=IDs, fill=Cond)) + geom_bar(stat="identity", colour="black", alpha =0.7) +
    geom_text(data=df_cover, aes(x=as.factor(Cond), y=IDs+300, label=IDs), col='black', size=4.5) +
    scale_fill_manual(values = c("grey15", "grey50", "grey85")) +
    theme_classic() + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=14),
          axis.ticks.x = element_blank(),
          # plot.title = element_text(size=14, hjust = 0.5),
          plot.subtitle = element_text(size = 12, hjust = 0.5),
          plot.title = element_text(size = 13.5, hjust = 0.5, face="bold"),
          legend.position = "none",
          strip.text.x = element_text(size=12, face = "italic"),
          legend.text = element_text(size=12))+
    labs(x = "", y = "Protein Ratios (A/B)", fill = "")
  
  
  ####################
  ####################
  ####################
  
  #intersect LF and multiDIA
  protein.groups_n_omit_b <- pD_list[[2]]
  protein.groups_n_omit_LF_b <- LF_list[[2]]
  
  protein.groups_n_omit_LF_b <- protein.groups_n_omit_LF_b[!grepl("H.", protein.groups_n_omit_LF_b$species),]
  protein.groups_n_omit_b <- protein.groups_n_omit_b[!grepl("H.", protein.groups_n_omit_b$species),]
  
  protein.groups_n_omit_LF_b_noInt <- protein.groups_n_omit_LF_b
  protein.groups_n_omit_b_noInt <- protein.groups_n_omit_b
  
  protein.groups_n_omit_b <- protein.groups_n_omit_b[protein.groups_n_omit_b$Protein.Group%in%protein.groups_n_omit_LF_b$Protein.Group,]
  protein.groups_n_omit_LF_b <- protein.groups_n_omit_LF_b[protein.groups_n_omit_LF_b$Protein.Group%in%protein.groups_n_omit_b$Protein.Group,]
  
  LF_mDIA_d0d4 <- rbind(protein.groups_n_omit_b, protein.groups_n_omit_LF_b)
  
  ################################# plotting #########################
  
  
  x_d0d4 <- log2(LF_mDIA_d0d4$sample2_int) %>% sort()
  x_d0d4_val <- as.numeric(quantile(x_d0d4,probs=c(.0025,.9975)))
  y_d0d4 <- log2(LF_mDIA_d0d4$d0_d8) %>% sort()
  y_d0d4_val <- as.numeric(quantile(y_d0d4,probs=c(.0025,.9975)))
  ####################################
  #d0d4
  b2 <- ggplot(protein.groups_n_omit_LF_b, aes(x=log2(sample2_int), y=log2(d0_d8), color = species, alpha = species)) +
    geom_point() + 
    scale_color_manual(values = c("#619CFF", "#F8766D")) +
    geom_hline(yintercept=log2(2/3), linetype="dashed", 
               color = "#619CFF", size=1)+
    geom_hline(yintercept=log2(3), linetype="dashed", 
               color = "#F8766D", size=1) +
    scale_alpha_manual(name = "category", values = c(0.3, 0.08, 0.26),guide = "none") +
    labs(x = expression(paste(Log["2"],", C")), y="Log2, delta0:delta8",title = paste0(LF_list_name))+
    #subtitle = paste0(nrow(LF_PGs_lim)," Protein groups"))+ 
    geom_smooth(method = "loess", method.args = list(family = "symmetric")) +
    coord_cartesian(ylim = c(y_d0d4_val[1],y_d0d4_val[2]),
                    xlim = c(x_d0d4_val[1],x_d0d4_val[2])) +
    theme_classic()+
    theme(legend.position="none", 
          plot.title = element_text(hjust = 0.5, size=18),
          axis.ticks.y = element_blank(),
          plot.subtitle = element_text(size =12, hjust = 0.5),
          axis.title.y=element_blank(),
          axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=14),
          axis.text.y = element_blank())
  
  
  b4 <- ggplot(protein.groups_n_omit_LF_b, aes(x=species, y=log2(d0_d8), color = species))+
    geom_boxplot() +   
    scale_color_manual(values = c("#619CFF", "#F8766D")) +
    geom_hline(yintercept=log2(2/3), linetype="dashed", 
               color = "#619CFF", size=1)+
    geom_hline(yintercept=log2(3), linetype="dashed", 
               color = "#F8766D", size=1) +
    scale_x_discrete(labels=c("H. sapiens" = "", "E. coli" = "",
                              "S. cerevisiae" = "")) +
    labs(x="",title = paste0(LF_list_name))+ #subtitle="")+
    coord_cartesian(ylim = c(y_d0d4_val[1],y_d0d4_val[2])) +
    theme(axis.ticks = element_blank(),
          axis.title.y=element_blank(),
          axis.text.y = element_blank(),
          panel.background=element_blank(), panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_text(size=32),
          plot.subtitle = element_text(size =12, hjust = 0.5),
          plot.title = element_text(size =18, hjust = 0.5))
  
  
  ####################################
  #d0d4
  b1 <- ggplot(protein.groups_n_omit_b, aes(x=log2(sample2_int), y=log2(d0_d8), color = species, alpha = species)) +
    geom_point() + #scale_bolor_brewer(palette = "PuOr")+
    scale_color_manual(values = c("#619CFF", "#F8766D")) +
    geom_hline(yintercept=log2(2/3), linetype="dashed", 
               color = "#619CFF", size=1)+
    geom_hline(yintercept=log2(3), linetype="dashed", 
               color = "#F8766D", size=1) +
    scale_alpha_manual(name = "category", values = c(0.3, 0.08, 0.26),guide = "none") +
    labs(x = expression(paste(Log["2"],", C")), y=expression(paste(Log["2"],", A/C")),title = paste0(pD_list_name))+
    #subtitle = paste0(nrow(d0_d4_pt)," Protein groups"))+ 
    geom_smooth(method = "loess", method.args = list(family = "symmetric")) + 
    coord_cartesian(ylim = c(y_d0d4_val[1],y_d0d4_val[2]),
                    xlim = c(x_d0d4_val[1],x_d0d4_val[2])) +
    theme_classic()+
    theme(legend.position="none", 
          plot.title = element_text(size =18, hjust = 0.5),
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          plot.subtitle = element_text(size =12, hjust = 0.5),
          axis.title.x = element_text(size=14),
          axis.title.y = element_text(size=14),
          legend.text = element_text(face = "italic"))
  
  
  b3 <- ggplot(protein.groups_n_omit_b, aes(x=species, y=log2(d0_d8), color = species))+
    geom_boxplot() +   
    scale_color_manual(values = c("#619CFF", "#F8766D")) +
    geom_hline(yintercept=log2(2/3), linetype="dashed", 
               color = "#619CFF", size=1)+
    geom_hline(yintercept=log2(3), linetype="dashed", 
               color = "#F8766D", size=1) +
    scale_x_discrete(labels=c("H. sapiens" = "", "E. coli" = "",
                              "S. cerevisiae" = "")) +
    labs(x="",title = paste0(pD_list_name), color = "Species:")+
    coord_cartesian(ylim = c(y_d0d4_val[1],y_d0d4_val[2])) +# subtitle = "")+
    theme(axis.ticks = element_blank(),
          axis.title.y=element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          panel.background=element_blank(), panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.title.x = element_text(size=32),
          plot.background=element_blank(),
          plot.subtitle = element_text(size =12, hjust = 0.5),
          plot.title = element_text(size =18, hjust = 0.5),
          legend.position = "top",
          legend.text = element_text(size=14, face= "italic"),
          legend.title = element_text(size=14),
          legend.spacing.x = unit(1.0, 'line'))
  
  
  
  ################################## coverage of protein group ratios:
  joined <- protein.groups_n_omit_b[protein.groups_n_omit_b$Protein.Group%in%protein.groups_n_omit_LF_b$Protein.Group,]
  
  types <- c("MultiDIA", "Label-Free", "Intersected")
  numbers <- c(nrow(protein.groups_n_omit_b_noInt), nrow(protein.groups_n_omit_LF_b_noInt), nrow(joined))
  df_cover <- data.frame("Cond" = types, "IDs" = numbers)
  
  df_cover$Cond <- factor(df_cover$Cond, levels=c("MultiDIA", "Label-Free", "Intersected"))
  
  b_cover <- ggplot(df_cover, aes(x=as.factor(Cond), y=IDs, fill=Cond)) + geom_bar(stat="identity", colour="black", alpha =0.7) +
    geom_text(data=df_cover, aes(x=as.factor(Cond), y=IDs+100, label=IDs), col='black', size=4.5) +
    scale_fill_manual(values = c("grey15", "grey50", "grey85")) +
    #scale_fill_manual(values = c("lightblue2", "lightpink2", "plum2")) +
    theme_classic() + #ylim(0,6000)+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=14),
          axis.ticks.x = element_blank(),
          # plot.title = element_text(size=14, hjust = 0.5),
          plot.subtitle = element_text(size = 12, hjust = 0.5),
          plot.title = element_text(size = 13.5, hjust = 0.5, face="bold"),
          legend.position = "none",
          strip.text.x = element_text(size=12, face = "italic"),
          legend.text = element_text(size=12))+
    #scale_fill_manual(values = c("lightslategrey", "darkorange", "mediumpurple1")) +
    labs(x = "", y = "Protein Ratios (A/C)", fill = "", title = expression(paste(italic(H.~sapiens),' excluded')))
  
  ################
  #################
  #################
  protein.groups_n_omit_c <- pD_list[[3]]
  protein.groups_n_omit_LF_c <- LF_list[[3]]
  
  protein.groups_n_omit_LF_c <- protein.groups_n_omit_LF_c[!grepl("H.", protein.groups_n_omit_LF_c$species),]
  protein.groups_n_omit_c <- protein.groups_n_omit_c[!grepl("H.", protein.groups_n_omit_c$species),]
  
  protein.groups_n_omit_noint <- protein.groups_n_omit_c
  protein.groups_n_omit_LF_noint <- protein.groups_n_omit_LF_c
  
  protein.groups_n_omit_c <- protein.groups_n_omit_c[protein.groups_n_omit_c$Protein.Group%in%protein.groups_n_omit_LF_c$Protein.Group,]
  protein.groups_n_omit_LF_c <- protein.groups_n_omit_LF_c[protein.groups_n_omit_LF_c$Protein.Group%in%protein.groups_n_omit_c$Protein.Group,]
  
  LF_mDIA_d0d4 <- rbind(protein.groups_n_omit_c, protein.groups_n_omit_LF_c)
  
  ################################# plotting #########################
  
  x_d0d4 <- log2(LF_mDIA_d0d4$sample2_int) %>% sort()
  x_d0d4_val <- as.numeric(quantile(x_d0d4,probs=c(.0025,.9975)))
  y_d0d4 <- log2(LF_mDIA_d0d4$d4_d8) %>% sort()
  y_d0d4_val <- as.numeric(quantile(y_d0d4,probs=c(.0025,.9975)))
  ####################################
  #d0d4
  c2 <- ggplot(protein.groups_n_omit_LF_c, aes(x=log2(sample2_int), y=log2(d4_d8), color = species, alpha = species)) +
    geom_point() + 
    scale_color_manual(values = c("#619CFF", "#F8766D")) +
    geom_hline(yintercept=log2(1/6), linetype="dashed", 
               color = "#619CFF", size=1)+
    geom_hline(yintercept=log2(6), linetype="dashed", 
               color = "#F8766D", size=1) +
    scale_alpha_manual(name = "category", values = c(0.3, 0.08, 0.26),guide = "none") +
    labs(x = expression(paste(Log["2"],", C")), y="Log2, delta4:delta8",title = paste0(LF_list_name))+
    #subtitle = paste0(nrow(LF_PGs_lim)," Protein groups"))+ 
    geom_smooth(method = "loess", method.args = list(family = "symmetric")) + 
    coord_cartesian(ylim = c(y_d0d4_val[1],y_d0d4_val[2]),
                    xlim = c(x_d0d4_val[1],x_d0d4_val[2])) +
    theme_classic()+ 
    theme(legend.position="none", 
          plot.title = element_text(hjust = 0.5, size=18),
          axis.ticks.y = element_blank(),
          plot.subtitle = element_text(size =12, hjust = 0.5),
          axis.title.y=element_blank(),
          axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=14),
          axis.text.y = element_blank())
  
  
  c4 <- ggplot(protein.groups_n_omit_LF_c, aes(x=species, y=log2(d4_d8), color = species))+
    geom_boxplot() +   
    scale_color_manual(values = c("#619CFF", "#F8766D")) +
    geom_hline(yintercept=log2(1/6), linetype="dashed", 
               color = "#619CFF", size=1)+
    geom_hline(yintercept=log2(6), linetype="dashed", 
               color = "#F8766D", size=1) +
    scale_x_discrete(labels=c("H. sapiens" = "", "E. coli" = "",
                              "S. cerevisiae" = "")) +
    labs(x="",title = paste0(LF_list_name))+
    coord_cartesian(ylim = c(y_d0d4_val[1],y_d0d4_val[2])) +#subtitle="")+
    theme(axis.ticks = element_blank(),
          axis.title.y=element_blank(),
          axis.text.y = element_blank(),
          panel.background=element_blank(), panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_text(size=32),
          plot.subtitle = element_text(size =12, hjust = 0.5),
          plot.title = element_text(size =18, hjust = 0.5))
  
  
  ####################################
  #d0d4
  c1 <- ggplot(protein.groups_n_omit_c, aes(x=log2(sample2_int), y=log2(d4_d8), color = species, alpha = species)) +
    geom_point() + #scale_color_brewer(palette = "PuOr")+
    scale_color_manual(values = c("#619CFF", "#F8766D")) +
    geom_hline(yintercept=log2(1/6), linetype="dashed", 
               color = "#619CFF", size=1)+
    geom_hline(yintercept=log2(6), linetype="dashed", 
               color = "#F8766D", size=1) +
    scale_alpha_manual(name = "category", values = c(0.3, 0.08, 0.26),guide = "none") +
    labs(x = expression(paste(Log["2"],", C")), y=expression(paste(Log["2"],", B/C")),title = paste0(pD_list_name))+
    #subtitle = paste0(nrow(d0_d4_pt)," Protein groups"))+ 
    geom_smooth(method = "loess", method.args = list(family = "symmetric")) +
    coord_cartesian(ylim = c(y_d0d4_val[1],y_d0d4_val[2]),
                    xlim = c(x_d0d4_val[1],x_d0d4_val[2])) +
    theme_classic()+
    theme(legend.position="none", 
          plot.title = element_text(size =18, hjust = 0.5),
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          plot.subtitle = element_text(size =12, hjust = 0.5),
          axis.title.x = element_text(size=14),
          axis.title.y = element_text(size=14),
          legend.text = element_text(face = "italic"))
  
  
  c3 <- ggplot(protein.groups_n_omit_c, aes(x=species, y=log2(d4_d8), color = species))+
    geom_boxplot() +   
    scale_color_manual(values = c("#619CFF", "#F8766D")) +
    geom_hline(yintercept=log2(1/6), linetype="dashed", 
               color = "#619CFF", size=1)+
    geom_hline(yintercept=log2(6), linetype="dashed", 
               color = "#F8766D", size=1) +
    scale_x_discrete(labels=c("H. sapiens" = "", "E. coli" = "",
                              "S. cerevisiae" = "")) +
    labs(x="",title = paste0(pD_list_name), color = "Species:")+
    coord_cartesian(ylim = c(y_d0d4_val[1],y_d0d4_val[2])) +# subtitle = "")+
    theme(axis.ticks = element_blank(),
          axis.title.y=element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          panel.background=element_blank(), panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.title.x = element_text(size=32),
          plot.background=element_blank(),
          plot.subtitle = element_text(size =12, hjust = 0.5),
          plot.title = element_text(size =18, hjust = 0.5),
          legend.position = "top",
          legend.text = element_text(size=14, face= "italic"),
          legend.title = element_text(size=14),
          legend.spacing.x = unit(1.0, 'line'))
  
  
  ################################## coverage of protein group ratios:
  joined <- protein.groups_n_omit_c[protein.groups_n_omit_c$Protein.Group%in%protein.groups_n_omit_LF_c$Protein.Group,]
  
  types <- c("plexDIA", "LF-DIA", "Intersected")
  numbers <- c(nrow(protein.groups_n_omit_noint), nrow(protein.groups_n_omit_LF_noint), nrow(joined))
  df_cover <- data.frame("Cond" = types, "IDs" = numbers)
  
  df_cover$Cond <- factor(df_cover$Cond, levels=c("plexDIA", "LF-DIA", "Intersected"))
  
  c_cover <- ggplot(df_cover, aes(x=as.factor(Cond), y=IDs, fill=Cond)) + geom_bar(stat="identity", colour="black", alpha =0.7) +
    geom_text(data=df_cover, aes(x=as.factor(Cond), y=IDs+100, label=IDs), col='black', size=4.5) +
    scale_fill_manual(values = c("grey15", "grey50", "grey85")) +
    #scale_fill_manual(values = c("lightblue2", "lightpink2", "plum2")) +
    theme_classic() + #ylim(0,6000)+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=14),
          axis.ticks.x = element_blank(),
          # plot.title = element_text(size=14, hjust = 0.5),
          plot.subtitle = element_text(size = 12, hjust = 0.5),
          plot.title = element_text(size = 13.5, hjust = 0.5, face="bold"),
          legend.position = "top",
          legend.direction = "vertical",
          strip.text.x = element_text(size=12, face = "italic"),
          legend.text = element_text(size=12))+
    #scale_fill_manual(values = c("lightslategrey", "darkorange", "mediumpurple1")) +
    labs(x = "", y = "Protein Ratios (B/C)", fill = "", title = expression(paste(italic(H.~sapiens)," excluded")))
  
  ############## combine
  g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}
  
  mylegend<-g_legend(a3)
  mylegend1<-g_legend(c_cover)
  
  
  a_fin <- grid.arrange(mylegend, nrow=2,heights=c(1, 10),
                        arrangeGrob(a_cover,
                                    a1 + theme(legend.position="none"),
                                    a2 + theme(legend.position="none"), 
                                    a3 + theme(legend.position="none"),
                                    a4 + theme(legend.position="none"), nrow = 1, widths=c(1.68,2.25,2,0.7,0.7)))
  
  c_fin <- grid.arrange(nrow=1,
                        arrangeGrob(c_cover + theme(legend.position="none"),
                                    c1 + theme(legend.position="none"),
                                    c2 + theme(legend.position="none"), 
                                    c3 + theme(legend.position="none"),
                                    c4 + theme(legend.position="none"), nrow = 1, widths=c(1.68,2.25,2,0.7,0.7)))
  
  b_fin <- grid.arrange(nrow=1,
                        arrangeGrob(b_cover + theme(legend.position="none"),
                                    b1 + theme(legend.position="none"),
                                    b2 + theme(legend.position="none"), 
                                    b3 + theme(legend.position="none"),
                                    b4 + theme(legend.position="none"), nrow = 1, widths=c(1.68,2.25,2,0.7,0.7)))
  
  
  figure_fin <- ggarrange(a_fin, NULL, b_fin, NULL, c_fin,
                          ncol = 1, nrow = 5,
                          heights=c(1.11, 0.1, 1, 0.1, 1))
  
  #ggsave("Figure3abc.pdf", figure_fin, width=11,height=9, dpi=400)
  
  return(figure_fin)
  
  
}

plot_XIC <- function(report, XIC, sample1 = "PDAC", 
                     sample2 = "U-937", Run = "wJD1166", Protein_x = "GLO1",
                     Modified.Sequence_1 = "\\(mTRAQ-n-4\\)SLDFYTR2",
                     Modified.Sequence_2 = "\\(mTRAQ-n-8\\)SLDFYTR2",
                     Stripped.Sequence = "SLDFYTR",
                     RT_start = 33.25,
                     RT_end = 33.5,
                     Remove_yions=F){
  
  XIC_Seq <- XIC[grepl(paste0(Run), XIC$File.Name),]
  #XIC_Seq <- pD_Cterm(XIC_Seq)
  if(Remove_yions==TRUE){
    XIC_Seq <- XIC_Seq[!grepl("y", XIC_Seq$FragmentType),] #remove y ions
  }
  report1 <- report[grepl(paste0(Run), report$Run),]
  #############
  Seq_frag_4 <- report1[grepl(paste0(Modified.Sequence_1), report1$Precursor.Id),]
  Seq_frag_4_cor <- unlist(strsplit(Seq_frag_4$Fragment.Correlations, ";"))
  Seq_frag_4_info <- unlist(strsplit(Seq_frag_4$Fragment.Info, ";"))
  Seq_frag_4_info <- sub("\\/.*","\\",Seq_frag_4_info)
  Seq_frag_4_info <- str_replace(Seq_frag_4_info, "\\^", "_")
  
  Seq_frag_8 <- report1[grepl(paste0(Modified.Sequence_2), report1$Precursor.Id),]
  Seq_frag_8_cor <- unlist(strsplit(Seq_frag_8$Fragment.Correlations, ";"))
  Seq_frag_8_info <- unlist(strsplit(Seq_frag_8$Fragment.Info, ";"))
  Seq_frag_8_info <- sub("\\/.*","\\",Seq_frag_8_info)
  Seq_frag_8_info <- str_replace(Seq_frag_8_info, "\\^", "_")
  
  
  Seq_cor <- data.frame(c(Seq_frag_4_cor, Seq_frag_8_cor))
  colnames(Seq_cor) <- "Cor"
  Seq_info <- data.frame(c(Seq_frag_4_info, Seq_frag_8_info))
  colnames(Seq_info) <- "ion"
  Seq_merge <- cbind(Seq_info, Seq_cor)
  if(Remove_yions==TRUE){
    Seq_merge <- Seq_merge[!grepl("y", Seq_merge$ion),]  #remove y ions if plexDIA superposition
  }
  Seq_merge <- Seq_merge %>% dplyr::group_by(ion) %>% dplyr::mutate("mean_cor" = mean(as.numeric(Cor))) %>% ungroup()
  Seq_merge <- Seq_merge[grepl("b|y", Seq_merge$ion),]
  #Seq_keep_ions <- Seq_merge %>% filter(mean_cor > 0.5)
  Seq_keep_ions <- Seq_merge %>% distinct(ion, .keep_all=T) %>% arrange(desc(mean_cor))  %>% top_n(4) 
  #Seq_keep_ions[6,] <- Seq_merge[3,]
  
  ############# d4
  Seq <- XIC_Seq[grepl(paste0(Modified.Sequence_1), XIC_Seq$Precursor.Id),]
  Seq$ion <- paste0(Seq$FragmentType,Seq$FragmentSeriesNumber,"_",Seq$FragmentCharge)
  #Seq_lim <- Seq[,grep("^ion|^X", colnames(Seq))]
  Seq_lim <- Seq %>% dplyr::select(ion,X0,X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,X16,X17,X18,X19,X20,X21,X22,X23,X24)
  
  mylist <- list() #create an empty list
  
  for(i in 1:12){  
    temp <- Seq_lim[(2*i+1):(2*i+2),]
    ion_temp <- paste0(temp[2,1])
    temp <- temp[,-1] #remove "ion"
    temp_t <- t(as.matrix(temp))
    colnames(temp_t) <- c("RT", "Intensity")
    temp_t <- data.frame(temp_t)
    temp_t$Ion <- as.character(paste0(ion_temp))
    temp_t$Quant <- "MS2"
    mylist[[i]] <- temp_t
  }
  
  df_Seq <- do.call("rbind",mylist) #combine all vectors into a matrix
  df_Seq$ion_type <- substr(df_Seq$Ion,1,1)
  df_Seq <- df_Seq[df_Seq$Ion%in%Seq_keep_ions$ion,]
  
  ### MS1
  Seq_MS1 <- Seq_lim[1:2,]
  Seq_MS1 <- Seq_MS1[,-1] #remove "ion"
  Seq_MS1_t <- t(as.matrix(Seq_MS1))
  colnames(Seq_MS1_t) <- c("RT", "Intensity")
  Seq_MS1_t <- data.frame(Seq_MS1_t)
  Seq_MS1_t$Ion <- "Precursor"
  Seq_MS1_t$Quant <- "MS1"
  Seq_MS1_t$ion_type <- "Precursor"
  
  Seq_all_d4 <- rbind(df_Seq, Seq_MS1_t)
  
  ############# d8
  Seq <- XIC_Seq[grepl(paste0(Modified.Sequence_2), XIC_Seq$Precursor.Id),]
  Seq$ion <- paste0(Seq$FragmentType,Seq$FragmentSeriesNumber,"_",Seq$FragmentCharge)
  #Seq_lim <- Seq[,grep("^ion|^X", colnames(Seq))]
  Seq_lim <- Seq %>% dplyr::select(ion,X0,X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,X16,X17,X18,X19,X20,X21,X22,X23,X24)
  
  mylist <- list() #create an empty list
  
  for(i in 1:12){  
    temp <- Seq_lim[(2*i+1):(2*i+2),]
    ion_temp <- paste0(temp[2,1])
    temp <- temp[,-1] #remove "ion"
    temp_t <- t(as.matrix(temp))
    colnames(temp_t) <- c("RT", "Intensity")
    temp_t <- data.frame(temp_t)
    temp_t$Ion <- as.character(paste0(ion_temp))
    temp_t$Quant <- "MS2"
    mylist[[i]] <- temp_t
  }
  
  df_Seq <- do.call("rbind",mylist) #combine all vectors into a matrix
  df_Seq$ion_type <- substr(df_Seq$Ion,1,1)
  df_Seq <- df_Seq[df_Seq$Ion%in%Seq_keep_ions$ion,]
  
  ### MS1
  Seq_MS1 <- Seq_lim[1:2,]
  Seq_MS1 <- Seq_MS1[,-1] #remove "ion"
  Seq_MS1_t <- t(as.matrix(Seq_MS1))
  colnames(Seq_MS1_t) <- c("RT", "Intensity")
  Seq_MS1_t <- data.frame(Seq_MS1_t)
  Seq_MS1_t$Ion <- "Precursor"
  Seq_MS1_t$Quant <- "MS1"
  Seq_MS1_t$ion_type <- "Precursor"
  
  Seq_all_d8 <- rbind(df_Seq, Seq_MS1_t)
  
  ######### combine:
  Seq_all_d4$Channel <- paste0(sample1)
  Seq_all_d8$Channel <- paste0(sample2)
  Seq_fin <- rbind(Seq_all_d4, Seq_all_d8)
  Seq_fin$Protein <- paste0(Protein_x)
  Seq_fin$Peptide <- paste0(Stripped.Sequence)
  
  bp <- Seq_fin
  bp$Channel<- factor(bp$Channel, level = c(paste0(sample1), paste0(sample2)))
  #bp$Intensity <- as.numeric(bp$Intensity)
  #bp$RT <- as.numeric(bp$RT)
  
  #############
  
  bp <- bp %>% dplyr::mutate(int_mirror = ifelse(grepl(paste0(sample1), Channel), -Intensity, Intensity))
  bp$ion_chan <- paste0(bp$Ion,bp$Channel)
  bp2 <- bp #%>% left_join(bp1, by =c("ion_chan" = "ion_chan"))
  MS1_int <- bp[grepl("Precursor", bp$ion_type),] %>% filter(RT>RT_start & RT<RT_end)
  MS1_int <- max(abs(MS1_int$Intensity))
  
  MS2_int <- bp[!grepl("Precursor", bp$ion_type),] %>% filter(RT>RT_start & RT<RT_end)
  MS2_int <- max(abs(MS2_int$Intensity))
  
  #scales_x <- list(
  #Protein_x = scale_x_continuous(limits = c(RT_start, RT_end), breaks = seq(RT_start, RT_end, length.out=3)))
  
  bp$ion_chan_quant <- paste0(bp$ion_chan, bp$Quant)
  bp <- bp %>% group_by(Protein, Quant) %>% mutate("norm_int" = Intensity/max(Intensity)) %>% ungroup()
  
  bp1 <- bp
  
  bp_int <- expand.grid(
    RT            = seq(from=min(bp$RT), max(bp$RT), by=.0005), #Decide resolution here.
    ion_chan              = unique(bp$ion_chan),
    stringsAsFactors = FALSE
  )
  
  bp1 <- bp %>% dplyr::full_join(bp_int, by=c("ion_chan", "RT")) %>%
    #dplyr::right_join(bp_int, by=c("RT" = "RT"))# %>% 
    dplyr::group_by(ion_chan) %>% 
    dplyr::mutate(y_interpolated   = spline(x=RT, y=int_mirror  , xout=RT)$y)%>% 
    dplyr::ungroup()
  bp1 <- bp1 %>% dplyr::mutate("Quant" = ifelse(grepl("b|y", ion_chan), "MS2", "MS1"))
  bp1$Protein <- paste0(Protein_x)
  bp1_PDAC <- bp1[grepl(paste0(sample1), bp1$ion_chan),]
  bp1_PDAC$y_interpolated[which(bp1_PDAC$y_interpolated>0)] <- 0
  bp1_U937 <- bp1[grepl(paste0(sample2), bp1$ion_chan),]
  bp1_U937$y_interpolated[which(bp1_U937$y_interpolated<0)] <- 0
  bp1 <- rbind(bp1_PDAC, bp1_U937)
  bp1_alpha <- bp1 %>%na.omit() %>% filter(abs(int_mirror)>0) %>%
    dplyr::mutate(alpha1=(abs(int_mirror))/(max(abs(int_mirror)))) %>% 
    ungroup() %>% distinct(ion_chan, RT, .keep_all=T) %>% 
    dplyr::mutate("RT_ion_chan" = paste0(RT,ion_chan)) %>% 
    dplyr::select("RT_ion_chan", "alpha1")
  bp1$RT_ion_chan <- paste0(bp1$RT,bp1$ion_chan)
  bp1 <- bp1 %>%left_join(bp1_alpha, by =c("RT_ion_chan"="RT_ion_chan"))
  #bp1$alpha2 <- (bp1$alpha1*0.2) + 0.8 #set new alpha range between 0.2-1.0
  bp1_MS1 <- bp1[grepl("MS1", bp1$Quant),] %>% filter(RT>RT_start & RT<RT_end)
  bp1_MS2 <- bp1[grepl("MS2", bp1$Quant),] %>% filter(RT>RT_start & RT<RT_end)
  
  ####################
  MS1plot <- ggplot(bp1_MS1, aes(x=RT, y=-y_interpolated, group=ion_chan, color = norm_int))  + 
    geom_point(data = bp1_MS1 %>% dplyr::filter(abs(int_mirror)>0), aes(y=-int_mirror, alpha=alpha1), size=2.8) +
    scale_alpha(range = c(0.4, 1))+
    #geom_line(aes(color=abs(..y..)/(log2(abs(..y..)))), size=1.5, alpha=0.7) + 
    geom_line(aes(color=abs(..y..)),size=1.5, alpha=0.9) + 
    geom_hline(yintercept = 0, size=2)+
    facet_grid_sc(cols = vars(Protein), rows = vars(Quant), 
                  labeller = label_wrap_gen(width = 1, multi_line = TRUE))+
    theme_bw()  + ylim(-MS1_int-MS1_int*0.05,MS1_int+MS1_int*0.05) + 
    #xlim(RT_start, RT_end) +
    scale_x_continuous(limits = c(RT_start, RT_end), breaks = round(seq(RT_start, RT_end, length.out=3),2))+
    labs(x = "Retention time (min)", y = "Intensity", fill = "Log2, Median Abundance") + 
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 14),
          strip.text.y = element_text(size = 14),
          strip.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          panel.spacing.x = unit(1.3, "lines"),
          legend.position = "none") + scale_color_viridis(option="turbo")
  
  MS2plot <- ggplot(bp1_MS2, aes(x=RT, y=-y_interpolated, group=ion_chan, color = norm_int))  + 
    geom_point(data = bp1_MS2 %>% dplyr::filter(abs(int_mirror)>0), aes(y=-int_mirror, alpha=alpha1), size=2.8) +
    scale_alpha(range = c(0.4, 1))+
    geom_line(aes(color=abs(..y..)),size=1.5, alpha=0.9) +
    geom_hline(yintercept = 0, size=2)+
    facet_grid_sc(cols = vars(Protein), rows = vars(Quant), 
                  labeller = label_wrap_gen(width = 1, multi_line = TRUE))+
    theme_bw()   + ylim(-MS2_int-MS2_int*0.05,MS2_int+MS2_int*0.05) + 
    #xlim(RT_start, RT_end) +
    scale_x_continuous(limits = c(RT_start, RT_end), breaks = round(seq(RT_start, RT_end, length.out=3),2))+
    labs(x = "Retention time (min)", y = "Intensity", fill = "Log2, Median Abundance") + 
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 14),
          strip.text.y = element_text(size = 14),
          strip.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          panel.spacing.x = unit(1.3, "lines"),
          legend.position = "none") + scale_color_viridis(option="turbo")
  
  MS1MS2plot <- ggarrange(MS1plot, MS2plot, ncol=1)
  return(MS1MS2plot)
}

#GGabby function credit largely to: https://pascal-martin.netlify.app/post/nicer-scatterplot-in-gggally/
GGscatterPlot <- function(data, mapping, ..., 
                          method = "spearman") {
  
  q_colors <- 15 # for no particular reason
  v_colors <-  viridis(q_colors)
  #Get correlation coefficient
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)
  
  cor <- cor(x, y, method = method)
  #Assemble data frame
  df <- data.frame(x = x, y = y)
  # PCA
  nonNull <- x!=0 & y!=0
  dfpc <- prcomp(~x+y, df[nonNull,])
  df$cols <- predict(dfpc, df)[,1]
  
  coeff.tls <- COEF::tls(y ~ x +0, df)
  coeff.tls <- coeff.tls[[1]]
  
  # Define the direction of color range based on PC1 orientation:
  dfsum <- x+y
  colDirection <- ifelse(dfsum[which.max(df$cols)] < 
                           dfsum[which.min(df$cols)],
                         1,
                         -1)
  #Get 2D density for alpha
  dens2D <- MASS::kde2d(df$x, df$y)
  df$density <- fields::interp.surface(dens2D , 
                                       df[,c("x", "y")])
  
  if (any(df$density==0)) {
    mini2D = min(df$density[df$density!=0]) #smallest non zero value
    df$density[df$density==0] <- mini2D
  }
  #Prepare plot
  x_d0d4 <- df$x %>% sort()
  x_d0d4_val <- as.numeric(quantile(x_d0d4,probs=c(.005,.995)))
  y_d0d4 <- df$y %>% sort()
  y_d0d4_val <- as.numeric(quantile(y_d0d4,probs=c(.005,.995)))
  
  pp <- ggplot(df, aes(x=x, y=y, color = density)) +
    ggplot2::geom_point(shape=16, show.legend = FALSE, size=0.75) +
    #ggplot2::scale_color_viridis() +
    scale_color_gradientn(colours = v_colors) +
    ggplot2::scale_alpha(range = c(.05, .6)) +
    #ggplot2::geom_smooth(method = COEF::tls)+
    scale_fill_gradient2(low = "red", mid = "white", high="green", limits = c(0,0.8),
                         midpoint=0.40)+
    #geom_abline(intercept = coeff.tls[1], slope = coeff.tls[2], col="black",linetype =1, size=0.8 ) +
    ggplot2::geom_label(
      data = data.frame(
        xlabel = min(x_d0d4_val[1], na.rm = TRUE),
        ylabel = max(y_d0d4_val[2], na.rm = TRUE),
        lab = paste0("\u03c1", " = ",round(cor, digits = 2))),
      mapping = ggplot2::aes(x = xlabel, 
                             y = ylabel, 
                             label = lab, fill=cor),
      hjust = 0, vjust = 1,
      size = 3, fontface = "bold", 
      inherit.aes = FALSE # do not inherit anything from the ...
    ) + xlim(x_d0d4_val[1],x_d0d4_val[2])+ 
    ylim(y_d0d4_val[1],y_d0d4_val[2])+
    theme_classic()
  
  return(pp)
}









######## other misc (probably wont use)

MS1_mTRAQ_d0d4 <- function(ev2_0,ev2_4){
  
  ev2_04 <- ev2_0 %>% inner_join(ev2_4, by =c("seqcharge_file"="seqcharge_file"))
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm" = Ms1.Area/median(Ms1.Area)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm" = Ms1.Area_d4/median(Ms1.Area_d4)) %>% ungroup()
  ev2_04 <- ev2_04[which(ev2_04$Ms1.Area & ev2_04$Ms1.Area_d4 > 0),] #keep only non-zero values
  
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d4_avg_int" = mean((Ms1.Area +Ms1.Area_d4)/2)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_norm_avg" = mean(d0_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_norm_avg" = mean(d4_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d0_d4_rat_individual" = d0_norm/d4_norm)
  ev2_04 <- ev2_04 %>% mutate("d0_d4_rat" = d0_norm_avg/d4_norm_avg)
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("pep_rat_CV_d0d4" = sd(d0_d4_rat_individual, na.rm=T)/mean(d0_d4_rat_individual)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  
  ev2_04 <- ev2_04%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04$Protein.Names), "HUMAN", 
                                                ifelse(grepl("YEAST", ev2_04$Protein.Names), "YEAST", 
                                                       ifelse(grepl("ECOLI", ev2_04$Protein.Names), "ECOLI", "unknown"))))
  
  ev2_04 <- ev2_04[!grepl("unknown", ev2_04$species),]
  ev2_04$Quant <- "MS1"
  ev2_04$Label <- "mTRAQ"
  
  return(ev2_04)
  
}

pD_PrecRatios <- function(df, quant.header = "Ms1.Area_iso"){
  #df <- ev2_lim
  df <- ev
  df <- pD_seqcharge(df)
  df <- pD_channel(df)
  df <- pD_rmMixSpec(df)
  df$run_chan <- paste0(df$Run, "_", df$channel_name)
  df$seqcharge_file <- paste0(df$Run, "_", df$seqcharge)
  ev2_04 <- reshape2::dcast(df, seqcharge_file+Stripped.Sequence+seqcharge+Run+Protein.Names~channel_name, value.var = quant.header)
  ev2_04[(ev2_04==0)] <- NA
  ev2_04$d0_d4 <- ev2_04$mTRAQ0/ev2_04$mTRAQ4
  ev2_04$d0_d8 <- ev2_04$mTRAQ0/ev2_04$mTRAQ8
  ev2_04$d4_d8 <- ev2_04$mTRAQ4/ev2_04$mTRAQ8
  med_summary <- ev2_04 %>% dplyr::filter(grepl("HUMAN", Protein.Names)) %>%
    dplyr::group_by(Run) %>% 
    dplyr::summarise_at(c("d0_d4", "d0_d8", "d4_d8"), median, na.rm = TRUE)%>%
    reshape2::melt(.data)
  med_summary <- med_summary[-1,]
  colnames(med_summary) <- c("labs", "med")
  med_summary$scalar <- 1/as.numeric(med_summary$med)
  
  ev2_04_m <- reshape2::melt(ev2_04)
  ev2_04_m <- ev2_04_m[!grepl("mTRAQ", ev2_04_m$variable),]
  ev2_04_m <- ev2_04_m %>% left_join(med_summary, by =c("variable" = "labs"))
  ev2_04_m <- ev2_04_m[!is.na(ev2_04_m$value),]
  ev2_04_m$rat_norm <- ev2_04_m$value*ev2_04_m$scalar
  
  return(ev2_04_m)
}


########### functions for plotting

reorder_within <- function(x, by, within, fun = max, sep = "___", ...) {
  if (!is.list(within)) {
    within <- list(within)
  }
  new_x <- do.call(paste, c(list(x, sep = sep), within))
  stats::reorder(new_x, by, FUN = fun)
}


scale_x_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_x_discrete(labels = function(x) gsub(reg, "", x), ...)
}


scale_y_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_y_discrete(labels = function(x) gsub(reg, "", x), ...)
}





# Translates new (post 1.8.1 b12) channels to the old format
# Old Channels were denoted like (mTRAQ0) new ones are denoted like (mTRAQ-K-0)
translate_diann_channel_format <- function(.input_df, columns = c("Precursor.Id","Modified.Sequence")){
  if (length(columns) < 1){
    print('translate_diann_channel_format, no columns specified')
    return(.input_df)
  }
  
  if (nrow(.input_df) < 1){
    print('translate_diann_channel_format, dataframe is empty')
    return(.input_df)
  }
  
  # check if channel is in old format
  test_precursor <- .input_df[[columns[1]]][[1]]
  label_occurences <- str_count(test_precursor, 'mTRAQ-[a-zA-Z]-')
  if(label_occurences == 0){
    return(.input_df)
  }
  
  
  for (column in columns) {
    .input_df[[column]] = sapply(.input_df[[column]], .update_channel)
  }
  return(.input_df)
  
}

.update_channel <- function(sequence){
  groups <- str_match_all(sequence, "mTRAQ-([a-zA-Z])-([0-9]+)")
  
  if (length(groups) > 0 ){
    groups <- groups[[1]]
    
    for(i in 1:nrow(groups)){
      sequence <- str_replace_all(sequence, groups[i,1], paste0('mTRAQ',groups[i,3]))
    }
  }
  
  return(sequence)
}


.get_channel <- function(sequence){
  label = ''
  
  for (channel in c('0','4','8')) {
    mod = channel
    if (grepl( mod, sequence, fixed = TRUE)){
      label <- channel
    }
  }
  
  return(label)
  
}

.lf_channel <- function(sequence){
  sequence <- str_replace_all(sequence, paste0('\\(mTRAQ[0-9]\\)'), '')
  return(sequence)
}

.get_sequence <- function(sequence){
  sequence <- str_replace_all(sequence, paste0('[0-9]'), '')
  return(sequence)
}


clust <- function(x) {
  umelt <- tidyr::spread(x, prot_info, Prot_norm_int)
  rownames(umelt) <- umelt$C_Q
  umelt1 <- umelt
  rownames(umelt1) <- umelt1$C_Q
  #2) calculate euclidian dist
  t_umelt <- t(as.matrix(umelt1))[-c(1:3),]
  ord <- hclust( dist(t_umelt, method = "euclidean"), method = "ward.D" )$order
  umelt <- melt(t_umelt)
  umelt$Var1 <- factor(umelt$Var1, levels= rownames(t_umelt)[ord])
  return(umelt)
}


