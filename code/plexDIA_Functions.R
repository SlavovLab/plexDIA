

#write functions for MS1 quantitation, MS2 quantitation, and averaged MS1 MS2 qaunt


######## MS1.Area:

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

######


######## MS1.Area:

MS1_mTRAQ_d0d8 <- function(ev2_0,ev2_4){
  
  ev2_04 <- ev2_0 %>% inner_join(ev2_4, by =c("seqcharge_file"="seqcharge_file"))
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm" = Ms1.Area/median(Ms1.Area)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm" = Ms1.Area_d8/median(Ms1.Area_d8)) %>% ungroup()
  ev2_04 <- ev2_04[which(ev2_04$Ms1.Area & ev2_04$Ms1.Area_d8 > 0),] #keep only non-zero values
  
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d4_avg_int" = mean((Ms1.Area +Ms1.Area_d8)/2)) %>% ungroup()
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

###########


######## MS1.Area:

MS1_mTRAQ_d4d8 <- function(ev2_0,ev2_4){
  
  ev2_04 <- ev2_0 %>% inner_join(ev2_4, by =c("seqcharge_file"="seqcharge_file"))
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm" = Ms1.Area_d4/median(Ms1.Area_d4)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm" = Ms1.Area_d8/median(Ms1.Area_d8)) %>% ungroup()
  ev2_04 <- ev2_04[which(ev2_04$Ms1.Area_d4 & ev2_04$Ms1.Area_d8 > 0),] #keep only non-zero values
  
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d4_avg_int" = mean((Ms1.Area_d4 +Ms1.Area_d8)/2)) %>% ungroup()
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
################################################ MS2 quant.
MS2_t_mTRAQ_d0d4 <- function(ev2_0,ev2_4){
  
  ev2_04 <- ev2_0 %>% inner_join(ev2_4, by =c("seqcharge_file"="seqcharge_file"))
  #ev2_04 <- ev2_04[which(ev2_04$Precursor.Translated & ev2_04$Precursor.Translated_d4 > 0),] #keep only non-zero values
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm" = Precursor.Translated/median(Precursor.Translated)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm" = Precursor.Translated_d4/median(Precursor.Translated_d4)) %>% ungroup()
  ev2_04 <- ev2_04[which(ev2_04$Precursor.Translated & ev2_04$Precursor.Translated_d4 > 0),] #keep only non-zero values
  
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d4_avg_int" = mean((Precursor.Translated +Precursor.Translated_d4)/2)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_norm_avg" = mean(d0_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_norm_avg" = mean(d4_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d0_d4_rat_individual" = d0_norm/d4_norm)
  ev2_04 <- ev2_04 %>% mutate("d0_d4_rat" = d0_norm_avg/d4_norm_avg)
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("pep_rat_CV_d0d4" = sd(d0_d4_rat_individual, na.rm=T)/mean(d0_d4_rat_individual)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  
  
  ev2_04 <- ev2_04%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04$Protein.Names), "HUMAN", 
                                                ifelse(grepl("YEAST", ev2_04$Protein.Names), "YEAST", 
                                                       ifelse(grepl("ECOLI", ev2_04$Protein.Names), "ECOLI", "unknown"))))
  
  ev2_04 <- ev2_04[!grepl("unknown", ev2_04$species),]
  ev2_04$Quant <- "MS2_translated"
  ev2_04$Label <- "mTRAQ"
  
  return(ev2_04)
  
}

#
MS2_t_mTRAQ_d0d8 <- function(ev2_0,ev2_4){
  
  ev2_04 <- ev2_0 %>% inner_join(ev2_8, by =c("seqcharge_file"="seqcharge_file"))
  #ev2_04 <- ev2_04[which(ev2_04$Precursor.Translated & ev2_04$Precursor.Translated_d8 > 0),] #keep only non-zero values
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm" = Precursor.Translated/median(Precursor.Translated)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm" = Precursor.Translated_d8/median(Precursor.Translated_d8)) %>% ungroup()
  ev2_04 <- ev2_04[which(ev2_04$Precursor.Translated & ev2_04$Precursor.Translated_d8 > 0),] #keep only non-zero values
  
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d4_avg_int" = mean((Precursor.Translated +Precursor.Translated_d8)/2)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_norm_avg" = mean(d0_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_norm_avg" = mean(d4_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d0_d4_rat_individual" = d0_norm/d4_norm)
  ev2_04 <- ev2_04 %>% mutate("d0_d4_rat" = d0_norm_avg/d4_norm_avg)
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("pep_rat_CV_d0d4" = sd(d0_d4_rat_individual, na.rm=T)/mean(d0_d4_rat_individual)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  
  
  ev2_04 <- ev2_04%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04$Protein.Names), "HUMAN", 
                                                ifelse(grepl("YEAST", ev2_04$Protein.Names), "YEAST", 
                                                       ifelse(grepl("ECOLI", ev2_04$Protein.Names), "ECOLI", "unknown"))))
  
  ev2_04 <- ev2_04[!grepl("unknown", ev2_04$species),]
  ev2_04$Quant <- "MS2_translated"
  ev2_04$Label <- "mTRAQ"
  
  return(ev2_04)
  
}

#
MS2_t_mTRAQ_d4d8 <- function(ev2_0,ev2_4){
  
  ev2_4 <- ev2_lim[grepl("mTRAQ4",ev2_lim$Precursor.Id),]
  ev2_04 <- ev2_4 %>% inner_join(ev2_8, by =c("seqcharge_file"="seqcharge_file"))
 # ev2_04 <- ev2_04[which(ev2_04$Precursor.Translated & ev2_04$Precursor.Translated_d8 > 0),] #keep only non-zero values
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm" = Precursor.Translated/median(Precursor.Translated)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm" = Precursor.Translated_d8/median(Precursor.Translated_d8)) %>% ungroup()
  ev2_04 <- ev2_04[which(ev2_04$Precursor.Translated & ev2_04$Precursor.Translated_d8 > 0),] #keep only non-zero values
  
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d4_avg_int" = mean((Precursor.Translated +Precursor.Translated_d8)/2)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_norm_avg" = mean(d0_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_norm_avg" = mean(d4_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d0_d4_rat_individual" = d0_norm/d4_norm)
  ev2_04 <- ev2_04 %>% mutate("d0_d4_rat" = d0_norm_avg/d4_norm_avg)
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("pep_rat_CV_d0d4" = sd(d0_d4_rat_individual, na.rm=T)/mean(d0_d4_rat_individual)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  
  
  ev2_04 <- ev2_04%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04$Protein.Names), "HUMAN", 
                                                ifelse(grepl("YEAST", ev2_04$Protein.Names), "YEAST", 
                                                       ifelse(grepl("ECOLI", ev2_04$Protein.Names), "ECOLI", "unknown"))))
  
  ev2_04 <- ev2_04[!grepl("unknown", ev2_04$species),]
  ev2_04$Quant <- "MS2_translated"
  ev2_04$Label <- "mTRAQ"
  
  return(ev2_04)
  
}

################################################ MS1 and MS2 quant. averaged:

MS1MS2_mTRAQ_d0d4 <- function(ev2_0,ev2_4){
  
  ev2_04a <- ev2_0 %>% inner_join(ev2_4, by =c("seqcharge_file"="seqcharge_file"))
  #ev2_04 <- ev2_04[which((ev2_04$Precursor.Translated & ev2_04$Precursor.Translated_d4 > 0) | (ev2_04$Ms1.Translated & ev2_04$Ms1.Translated_d4 > 0)),] #keep only non-zero values
  ev2_04 <- ev2_04a[which((ev2_04a$Precursor.Translated & ev2_04a$Precursor.Translated_d4 > 0)),] #keep only non-zero values
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm_MS2" = Precursor.Translated/median(Precursor.Translated)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm_MS2" = Precursor.Translated_d4/median(Precursor.Translated_d4)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d4_avg_int_MS2" = mean((Precursor.Translated +Precursor.Translated_d4)/2)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_norm_avg_MS2" = mean(d0_norm_MS2)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_norm_avg_MS2" = mean(d4_norm_MS2)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d0_d4_rat_MS2" = d0_norm_avg_MS2/d4_norm_avg_MS2)
  ev2_04_1 <- ev2_04a[which((ev2_04a$Ms1.Translated & ev2_04a$Ms1.Translated_d4 > 0)),] #keep only non-zero values
  ev2_04_1 <- ev2_04_1 %>% group_by(Run) %>% mutate("d0_norm_MS1" = Ms1.Translated/median(Ms1.Translated)) %>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% group_by(Run) %>% mutate("d4_norm_MS1" = Ms1.Translated_d4/median(Ms1.Translated_d4)) %>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% group_by(seqcharge) %>% mutate("d0_d4_avg_int_MS1" = mean((Ms1.Translated +Ms1.Translated_d4)/2)) %>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% group_by(seqcharge) %>% mutate("d0_norm_avg_MS1" = mean(d0_norm_MS1)) %>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% group_by(seqcharge) %>% mutate("d4_norm_avg_MS1" = mean(d4_norm_MS1)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% mutate("d0_d4_rat_MS1" = d0_norm_avg_MS1/d4_norm_avg_MS1)
  ev2_04_1 <- ev2_04_1 %>% dplyr::select(Protein.Names,seqcharge, d0_d4_avg_int_MS1, d0_d4_rat_MS1)
  
  ev2_04_both <- ev2_04 %>% full_join(ev2_04_1, by =c("seqcharge" = "seqcharge"))
  ev2_04_both$d0_d4_rat <- rowMeans(subset(ev2_04_both, select = c(d0_d4_rat_MS2, d0_d4_rat_MS1)), na.rm = TRUE)
  ev2_04_both$d0_d4_avg_int <- rowMeans(subset(ev2_04_both, select = c(d0_d4_avg_int_MS2, d0_d4_avg_int_MS1)), na.rm = TRUE)
  
  
  ev2_04_both <- ev2_04_both%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04_both$Protein.Names.x) | grepl("HUMAN", ev2_04_both$Protein.Names.y), "HUMAN", 
                                                ifelse(grepl("YEAST", ev2_04_both$Protein.Names.x) | grepl("YEAST",ev2_04_both$Protein.Names.y), "YEAST", 
                                                       ifelse(grepl("ECOLI", ev2_04_both$Protein.Names.x) | grepl("ECOLI", ev2_04_both$Protein.Names.y), "ECOLI", "unknown"))))
  
  ev2_04_both <- ev2_04_both[!grepl("unknown", ev2_04_both$species),]
  ev2_04_both$Quant <- "MS2_MS1"
  ev2_04_both$Label <- "mTRAQ"
  
  return(ev2_04_both)
  
}

#d0_d8
MS1MS2_mTRAQ_d0d8 <- function(ev2_0,ev2_4){
  
  ev2_04a <- ev2_0 %>% inner_join(ev2_4, by =c("seqcharge_file"="seqcharge_file"))
  #ev2_04 <- ev2_04[which((ev2_04$Precursor.Translated & ev2_04$Precursor.Translated_d8 > 0) | (ev2_04$Ms1.Translated & ev2_04$Ms1.Translated_d8 > 0)),] #keep only non-zero values
  ev2_04 <- ev2_04a[which((ev2_04a$Precursor.Translated & ev2_04a$Precursor.Translated_d8 > 0)),] #keep only non-zero values
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm_MS2" = Precursor.Translated/median(Precursor.Translated)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d8_norm_MS2" = Precursor.Translated_d8/median(Precursor.Translated_d8)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d8_avg_int_MS2" = mean((Precursor.Translated +Precursor.Translated_d8)/2)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_norm_avg_MS2" = mean(d0_norm_MS2)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d8_norm_avg_MS2" = mean(d8_norm_MS2)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d0_d8_rat_MS2" = d0_norm_avg_MS2/d8_norm_avg_MS2)
  ev2_04_1 <- ev2_04a[which((ev2_04a$Ms1.Translated & ev2_04a$Ms1.Translated_d8 > 0)),] #keep only non-zero values
  ev2_04_1 <- ev2_04_1 %>% group_by(Run) %>% mutate("d0_norm_MS1" = Ms1.Translated/median(Ms1.Translated)) %>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% group_by(Run) %>% mutate("d8_norm_MS1" = Ms1.Translated_d8/median(Ms1.Translated_d8)) %>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% group_by(seqcharge) %>% mutate("d0_d8_avg_int_MS1" = mean((Ms1.Translated +Ms1.Translated_d8)/2)) %>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% group_by(seqcharge) %>% mutate("d0_norm_avg_MS1" = mean(d0_norm_MS1)) %>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% group_by(seqcharge) %>% mutate("d8_norm_avg_MS1" = mean(d8_norm_MS1)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% mutate("d0_d8_rat_MS1" = d0_norm_avg_MS1/d8_norm_avg_MS1)
  ev2_04_1 <- ev2_04_1 %>% dplyr::select(Protein.Names,seqcharge, d0_d8_avg_int_MS1, d0_d8_rat_MS1)
  
  ev2_04_both <- ev2_04 %>% full_join(ev2_04_1, by =c("seqcharge" = "seqcharge"))
  ev2_04_both$d0_d8_rat <- rowMeans(subset(ev2_04_both, select = c(d0_d8_rat_MS2, d0_d8_rat_MS1)), na.rm = TRUE)
  ev2_04_both$d0_d8_avg_int <- rowMeans(subset(ev2_04_both, select = c(d0_d8_avg_int_MS2, d0_d8_avg_int_MS1)), na.rm = TRUE)
  
  
  ev2_04_both <- ev2_04_both%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04_both$Protein.Names.x) | grepl("HUMAN", ev2_04_both$Protein.Names.y), "HUMAN", 
                                                          ifelse(grepl("YEAST", ev2_04_both$Protein.Names.x) | grepl("YEAST",ev2_04_both$Protein.Names.y), "YEAST", 
                                                                 ifelse(grepl("ECOLI", ev2_04_both$Protein.Names.x) | grepl("ECOLI", ev2_04_both$Protein.Names.y), "ECOLI", "unknown"))))
  
  ev2_04_both <- ev2_04_both[!grepl("unknown", ev2_04_both$species),]
  ev2_04_both$Quant <- "MS2_MS1"
  ev2_04_both$Label <- "mTRAQ"
  
  return(ev2_04_both)

}
#d4_d8

MS1MS2_mTRAQ_d4d8 <- function(ev2_0,ev2_4){
  
  ev2_4 <- ev2_lim[grepl("Nter4|Lys4",ev2_lim$Precursor.Id),] 
  
  ev2_04 <- ev2_4 %>% inner_join(ev2_8, by =c("seqcharge_file"="seqcharge_file"))
  ev2_04 <- ev2_04[which((ev2_04$Precursor.Translated & ev2_04$Precursor.Translated_d8 > 0) | (ev2_04$Ms1.Translated & ev2_04$Ms1.Translated_d8 > 0)),] #keep only non-zero values
  ev2_04 <- ev2_04[which((ev2_04$Precursor.Translated & ev2_04$Precursor.Translated_d8 > 0)),] #keep only non-zero values
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm_MS2" = Precursor.Translated/median(Precursor.Translated)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d8_norm_MS2" = Precursor.Translated_d8/median(Precursor.Translated_d8)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_d8_avg_int_MS2" = mean((Precursor.Translated +Precursor.Translated_d8)/2)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_norm_avg_MS2" = mean(d4_norm_MS2)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d8_norm_avg_MS2" = mean(d8_norm_MS2)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d4_d8_rat_MS2" = d4_norm_avg_MS2/d8_norm_avg_MS2)
  #ev2_04_1 <- ev2_04a[which((ev2_04a$Ms1.Translated & ev2_04a$Ms1.Translated_d8 > 0)),] #keep only non-zero values
  ev2_04_1 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm_MS1" = Ms1.Translated/median(Ms1.Translated)) %>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% group_by(Run) %>% mutate("d8_norm_MS1" = Ms1.Translated_d8/median(Ms1.Translated_d8)) %>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% group_by(seqcharge) %>% mutate("d4_d8_avg_int_MS1" = mean((Ms1.Translated +Ms1.Translated_d8)/2)) %>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% group_by(seqcharge) %>% mutate("d4_norm_avg_MS1" = mean(d4_norm_MS1)) %>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% group_by(seqcharge) %>% mutate("d8_norm_avg_MS1" = mean(d8_norm_MS1)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  ev2_04_both <- ev2_04_1 %>% mutate("d4_d8_rat_MS1" = d4_norm_avg_MS1/d8_norm_avg_MS1)
  #ev2_04_1 <- ev2_04_1 %>% dplyr::select(Protein.Names,seqcharge, d4_d8_avg_int_MS1, d4_d8_rat_MS1)
  
 # ev2_04_both <- ev2_04 %>% full_join(ev2_04_1, by =c("seqcharge" = "seqcharge"))
  ev2_04_both$d4_d8_rat <- rowMeans(subset(ev2_04_both, select = c(d4_d8_rat_MS2, d4_d8_rat_MS1)), na.rm = TRUE)
  ev2_04_both$d4_d8_avg_int <- rowMeans(subset(ev2_04_both, select = c(d4_d8_avg_int_MS2, d4_d8_avg_int_MS1)), na.rm = TRUE)
  
  ev2_04_both <- ev2_04_both%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04_both$Protein.Names), "HUMAN", 
                                                          ifelse(grepl("YEAST", ev2_04_both$Protein.Names), "YEAST", 
                                                                 ifelse(grepl("ECOLI", ev2_04_both$Protein.Names), "ECOLI", "unknown"))))
  
  # 
  # ev2_04_both <- ev2_04_both%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04_both$Protein.Names.x) | grepl("HUMAN", ev2_04_both$Protein.Names.y), "HUMAN", 
  #                                                         ifelse(grepl("YEAST", ev2_04_both$Protein.Names.x) | grepl("YEAST",ev2_04_both$Protein.Names.y), "YEAST", 
  #                                                                ifelse(grepl("ECOLI", ev2_04_both$Protein.Names.x) | grepl("ECOLI", ev2_04_both$Protein.Names.y), "ECOLI", "unknown"))))
  # 
  ev2_04_both <- ev2_04_both[!grepl("unknown", ev2_04_both$species),]
  ev2_04_both$Quant <- "MS2_MS1"
  ev2_04_both$Label <- "mTRAQ"
  
  return(ev2_04_both)

}


###################################################################
#####################################################################
####################################################################  LABEL FREE

################################################ MS2 quant.
MS2_LF_d0d4 <- function(ev2_0,ev2_4){
  
  ev2_04 <- ev2_0 %>% inner_join(ev2_4, by =c("seqcharge_file"="seqcharge_file"))
  ev2_04 <- ev2_04[which(ev2_04$Precursor.Quantity & ev2_04$Precursor.Quantity_d4 > 0),] #keep only non-zero values
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm" = Precursor.Quantity/median(Precursor.Quantity)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm" = Precursor.Quantity_d4/median(Precursor.Quantity_d4)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d4_avg_int" = mean((Precursor.Quantity +Precursor.Quantity_d4)/2)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_norm_avg" = mean(d0_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_norm_avg" = mean(d4_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d0_d4_rat_individual" = d0_norm/d4_norm)
  ev2_04 <- ev2_04 %>% mutate("d0_d4_rat" = d0_norm_avg/d4_norm_avg)
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("pep_rat_CV_d0d4" = sd(d0_d4_rat_individual, na.rm=T)/mean(d0_d4_rat_individual)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  
  ev2_04 <- ev2_04%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04$Protein.Names), "HUMAN", 
                                                ifelse(grepl("YEAST", ev2_04$Protein.Names), "YEAST", 
                                                       ifelse(grepl("ECOLI", ev2_04$Protein.Names), "ECOLI", "unknown"))))
  
  ev2_04 <- ev2_04[!grepl("unknown", ev2_04$species),]
  ev2_04$Quant <- "MS2"
  ev2_04$Label <- "Unlabeled"
  
  return(ev2_04)
  
}

#
MS2_LF_d0d8 <- function(ev2_0,ev2_4){
  
  ev2_04 <- ev2_0 %>% inner_join(ev2_8, by =c("seqcharge_file"="seqcharge_file"))
  ev2_04 <- ev2_04[which(ev2_04$Precursor.Quantity & ev2_04$Precursor.Quantity_d8 > 0),] #keep only non-zero values
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm" = Precursor.Quantity/median(Precursor.Quantity)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm" = Precursor.Quantity_d8/median(Precursor.Quantity_d8)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d8_avg_int" = mean((Precursor.Quantity +Precursor.Quantity_d8)/2)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_norm_avg" = mean(d0_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_norm_avg" = mean(d4_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d0_d8_rat_individual" = d0_norm/d4_norm)
  ev2_04 <- ev2_04 %>% mutate("d0_d8_rat" = d0_norm_avg/d4_norm_avg)
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("pep_rat_CV_d0d8" = sd(d0_d8_rat_individual, na.rm=T)/mean(d0_d8_rat_individual)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  
  ev2_04 <- ev2_04%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04$Protein.Names), "HUMAN", 
                                                ifelse(grepl("YEAST", ev2_04$Protein.Names), "YEAST", 
                                                       ifelse(grepl("ECOLI", ev2_04$Protein.Names), "ECOLI", "unknown"))))
  
  ev2_04 <- ev2_04[!grepl("unknown", ev2_04$species),]
  ev2_04$Quant <- "MS2"
  ev2_04$Label <- "Unlabeled"
  
  return(ev2_04)
  
}

#
MS2_LF_d4d8 <- function(ev2_0,ev2_4){
  
  ev2_4 <- ev2_lim[grepl(paste0(RF1[2]),ev2_lim$Run),] #assumes the 2nd run is delta4 equivalent
  ev2_04 <- ev2_4 %>% inner_join(ev2_8, by =c("seqcharge_file"="seqcharge_file_d8"))
  ev2_04 <- ev2_04[which(ev2_04$Precursor.Quantity & ev2_04$Precursor.Quantity_d8 > 0),] #keep only non-zero values
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm" = Precursor.Quantity/median(Precursor.Quantity)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm" = Precursor.Quantity_d8/median(Precursor.Quantity_d8)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_d8_avg_int" = mean((Precursor.Quantity +Precursor.Quantity_d8)/2)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_norm_avg" = mean(d0_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_norm_avg" = mean(d4_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d4_d8_rat_individual" = d0_norm/d4_norm)
  ev2_04 <- ev2_04 %>% mutate("d4_d8_rat" = d0_norm_avg/d4_norm_avg)
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("pep_rat_CV_d4d8" = sd(d4_d8_rat_individual, na.rm=T)/mean(d4_d8_rat_individual)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  
  
  ev2_04 <- ev2_04%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04$Protein.Names), "HUMAN", 
                                                ifelse(grepl("YEAST", ev2_04$Protein.Names), "YEAST", 
                                                       ifelse(grepl("ECOLI", ev2_04$Protein.Names), "ECOLI", "unknown"))))
  
  ev2_04 <- ev2_04[!grepl("unknown", ev2_04$species),]
  ev2_04$Quant <- "MS2"
  ev2_04$Label <- "Unlabeled"
  
  return(ev2_04)
  
}


################################################ MS2 quant.
MS2_t_LF_d0d4 <- function(ev2_0,ev2_4){
  
  ev2_04 <- ev2_0 %>% inner_join(ev2_4, by =c("seqcharge"="seqcharge_d4"))
  #ev2_04 <- ev2_04[which(ev2_04$Precursor.Translated & ev2_04$Precursor.Translated_d4 > 0),] #keep only non-zero values
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm" = Precursor.Translated/median(Precursor.Translated)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm" = Precursor.Translated_d4/median(Precursor.Translated_d4)) %>% ungroup()
  ev2_04 <- ev2_04[which(ev2_04$Precursor.Translated & ev2_04$Precursor.Translated_d4 > 0),] #keep only non-zero values
  
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d4_avg_int" = mean((Precursor.Translated +Precursor.Translated_d4)/2)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_norm_avg" = mean(d0_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_norm_avg" = mean(d4_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d0_d4_rat_individual" = d0_norm/d4_norm)
  ev2_04 <- ev2_04 %>% mutate("d0_d4_rat" = d0_norm_avg/d4_norm_avg)
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("pep_rat_CV_d0d4" = sd(d0_d4_rat_individual, na.rm=T)/mean(d0_d4_rat_individual)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  
  
  ev2_04 <- ev2_04%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04$Protein.Names), "HUMAN", 
                                                ifelse(grepl("YEAST", ev2_04$Protein.Names), "YEAST", 
                                                       ifelse(grepl("ECOLI", ev2_04$Protein.Names), "ECOLI", "unknown"))))
  
  ev2_04 <- ev2_04[!grepl("unknown", ev2_04$species),]
  ev2_04$Quant <- "MS2_translated"
  ev2_04$Label <- "Unlabeled"
  
  return(ev2_04)
  
}

#
MS2_t_LF_d0d8 <- function(ev2_0,ev2_4){
  
  ev2_04 <- ev2_0 %>% inner_join(ev2_8, by =c("seqcharge"="seqcharge_d4"))
  #ev2_04 <- ev2_04[which(ev2_04$Precursor.Translated & ev2_04$Precursor.Translated_d8 > 0),] #keep only non-zero values
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm" = Precursor.Translated/median(Precursor.Translated)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm" = Precursor.Translated_d8/median(Precursor.Translated_d8)) %>% ungroup()
  ev2_04 <- ev2_04[which(ev2_04$Precursor.Translated & ev2_04$Precursor.Translated_d8 > 0),] #keep only non-zero values
  
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d8_avg_int" = mean((Precursor.Translated +Precursor.Translated_d8)/2)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_norm_avg" = mean(d0_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_norm_avg" = mean(d4_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d0_d8_rat_individual" = d0_norm/d4_norm)
  ev2_04 <- ev2_04 %>% mutate("d0_d8_rat" = d0_norm_avg/d4_norm_avg)
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("pep_rat_CV_d0d8" = sd(d0_d8_rat_individual, na.rm=T)/mean(d0_d8_rat_individual)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  
  
  ev2_04 <- ev2_04%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04$Protein.Names), "HUMAN", 
                                                ifelse(grepl("YEAST", ev2_04$Protein.Names), "YEAST", 
                                                       ifelse(grepl("ECOLI", ev2_04$Protein.Names), "ECOLI", "unknown"))))
  
  ev2_04 <- ev2_04[!grepl("unknown", ev2_04$species),]
  ev2_04$Quant <- "MS2_translated"
  ev2_04$Label <- "Unlabeled"
  
  return(ev2_04)
  
}

#
MS2_t_LF_d4d8 <- function(ev2_0,ev2_4){
  
  ev2_4 <- ev2_lim[grepl(paste0(RF1[2],"|",RF1[3]),ev2_lim$Run),] #assumes the 2nd run is delta4 equivalent
  ev2_04 <- ev2_4 %>% inner_join(ev2_8, by =c("seqcharge"="seqcharge_d8"))
  # ev2_04 <- ev2_04[which(ev2_04$Precursor.Translated & ev2_04$Precursor.Translated_d8 > 0),] #keep only non-zero values
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm" = Precursor.Translated/median(Precursor.Translated)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm" = Precursor.Translated_d8/median(Precursor.Translated_d8)) %>% ungroup()
  ev2_04 <- ev2_04[which(ev2_04$Precursor.Translated & ev2_04$Precursor.Translated_d8 > 0),] #keep only non-zero values
  
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_d8_avg_int" = mean((Precursor.Translated +Precursor.Translated_d8)/2)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_norm_avg" = mean(d0_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_norm_avg" = mean(d4_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d4_d8_rat_individual" = d0_norm/d4_norm)
  ev2_04 <- ev2_04 %>% mutate("d4_d8_rat" = d0_norm_avg/d4_norm_avg)
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("pep_rat_CV_d4d8" = sd(d4_d8_rat_individual, na.rm=T)/mean(d4_d8_rat_individual)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  
  
  ev2_04 <- ev2_04%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04$Protein.Names), "HUMAN", 
                                                ifelse(grepl("YEAST", ev2_04$Protein.Names), "YEAST", 
                                                       ifelse(grepl("ECOLI", ev2_04$Protein.Names), "ECOLI", "unknown"))))
  
  ev2_04 <- ev2_04[!grepl("unknown", ev2_04$species),]
  ev2_04$Quant <- "MS2_translated"
  ev2_04$Label <- "Unlabeled"
  
  return(ev2_04)
  
}

################################################ MS1 quant.
MS1_LF_d0d4 <- function(ev2_0,ev2_4){
  
  ev2_04 <- ev2_0 %>% inner_join(ev2_4, by =c("seqcharge_file"="seqcharge_file_d4"))
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm" = Ms1.Area/median(Ms1.Area)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm" = Ms1.Area_d4/median(Ms1.Area_d4)) %>% ungroup()
  ev2_04 <- ev2_04[which(ev2_04$Ms1.Area & ev2_04$Ms1.Area_d4 > 0),] #keep only non-zero values
  
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d4_avg_int" = mean((Ms1.Area +Ms1.Area_d4)/2)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_norm_avg" = mean(d0_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_norm_avg" = mean(d4_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d0_d4_rat_individual" = d0_norm/d4_norm)
  ev2_04 <- ev2_04 %>% mutate("d0_d4_rat" = d0_norm_avg/d4_norm_avg)
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("pep_rat_CV_d0d4" = sd(d0_d4_rat_individual, na.rm=T)/mean(d0_d4_rat_individual)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  
  
  #  ev2_04 <- ev2_04[which(ev2_04$Ms1.Area & ev2_04$Ms1.Area_d4 > 0),] #keep only non-zero values
  #  ev2_04$d0_norm <- ev2_04$Ms1.Area/median(ev2_04$Ms1.Area)
  #  ev2_04$d4_norm <- ev2_04$Ms1.Area_d4/median(ev2_04$Ms1.Area_d4)
  #  ev2_04$d0_d4_rat <- ev2_04$d0_norm/ev2_04$d4_norm
  #  ev2_04$d0_d4_avg_int <- (ev2_04$Ms1.Area+ev2_04$Ms1.Area_d4)/2
  
  ev2_04 <- ev2_04%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04$Protein.Names), "HUMAN", 
                                                ifelse(grepl("YEAST", ev2_04$Protein.Names), "YEAST", 
                                                       ifelse(grepl("ECOLI", ev2_04$Protein.Names), "ECOLI", "unknown"))))
  
  ev2_04 <- ev2_04[!grepl("unknown", ev2_04$species),]
  ev2_04$Quant <- "MS1"
  ev2_04$Label <- "Unlabeled"
  
  return(ev2_04)
  
}

#
MS1_LF_d0d8 <- function(ev2_0,ev2_4){
  
  ev2_04 <- ev2_0 %>% inner_join(ev2_4, by =c("seqcharge"="seqcharge_d8"))
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm" = Ms1.Area/median(Ms1.Area)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm" = Ms1.Area_d8/median(Ms1.Area_d8)) %>% ungroup()
  ev2_04 <- ev2_04[which(ev2_04$Ms1.Area & ev2_04$Ms1.Area_d8 > 0),] #keep only non-zero values
  
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d8_avg_int" = mean((Ms1.Area +Ms1.Area_d8)/2)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_norm_avg" = mean(d0_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_norm_avg" = mean(d4_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d0_d8_rat_individual" = d0_norm/d4_norm)
  ev2_04 <- ev2_04 %>% mutate("d0_d8_rat" = d0_norm_avg/d4_norm_avg)
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("pep_rat_CV_d0d8" = sd(d0_d8_rat_individual, na.rm=T)/mean(d0_d8_rat_individual)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  
  # ev2_04 <- ev2_0 %>% inner_join(ev2_4, by =c("seqcharge_file"="seqcharge_file"))
  # ev2_04 <- ev2_04[which(ev2_04$Ms1.Area & ev2_04$Ms1.Area_d8 > 0),] #keep only non-zero values
  # ev2_04$d0_norm <- ev2_04$Ms1.Area/median(ev2_04$Ms1.Area)
  # ev2_04$d8_norm <- ev2_04$Ms1.Area_d8/median(ev2_04$Ms1.Area_d8)
  # ev2_04$d0_d8_rat <- ev2_04$d0_norm/ev2_04$d8_norm
  # ev2_04$d0_d8_avg_int <- (ev2_04$Ms1.Area+ev2_04$Ms1.Area_d8)/2
  
  ev2_04 <- ev2_04%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04$Protein.Names), "HUMAN", 
                                                ifelse(grepl("YEAST", ev2_04$Protein.Names), "YEAST", 
                                                       ifelse(grepl("ECOLI", ev2_04$Protein.Names), "ECOLI", "unknown"))))
  
  ev2_04 <- ev2_04[!grepl("unknown", ev2_04$species),]
  ev2_04$Quant <- "MS1"
  ev2_04$Label <- "Unlabeled"
  
  return(ev2_04)
  
}


#
MS1_LF_d4d8 <- function(ev2_4,ev2_8){
  
  ev2_4 <- ev2_lim[grepl(paste0(RF2),ev2_lim$Run),] #assumes the 2nd run is delta4 equivalent
  ev2_04 <- ev2_4 %>% inner_join(ev2_8, by =c("seqcharge"="seqcharge_d8"))
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm" = Ms1.Area/median(Ms1.Area)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm" = Ms1.Area_d8/median(Ms1.Area_d8)) %>% ungroup()
  ev2_04 <- ev2_04[which(ev2_04$Ms1.Area & ev2_04$Ms1.Area_d8 > 0),] #keep only non-zero values
  
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_d8_avg_int" = mean((Ms1.Area +Ms1.Area_d8)/2)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_norm_avg" = mean(d0_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_norm_avg" = mean(d4_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d4_d8_rat_individual" = d0_norm/d4_norm)
  ev2_04 <- ev2_04 %>% mutate("d4_d8_rat" = d0_norm_avg/d4_norm_avg)
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("pep_rat_CV_d4d8" = sd(d4_d8_rat_individual, na.rm=T)/mean(d4_d8_rat_individual)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  
  # ev2_04 <- ev2_4 %>% inner_join(ev2_8, by =c("seqcharge_file"="seqcharge_file"))
  # ev2_04 <- ev2_04[which(ev2_04$Ms1.Area & ev2_04$Ms1.Area_d8 > 0),] #keep only non-zero values
  # ev2_04$d4_norm <- ev2_04$Ms1.Area/median(ev2_04$Ms1.Area)
  # ev2_04$d8_norm <- ev2_04$Ms1.Area_d8/median(ev2_04$Ms1.Area_d8)
  # ev2_04$d4_d8_rat <- ev2_04$d4_norm/ev2_04$d8_norm
  # ev2_04$d4_d8_avg_int <- (ev2_04$Ms1.Area+ev2_04$Ms1.Area_d8)/2
  
  ev2_04 <- ev2_04%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04$Protein.Names), "HUMAN", 
                                                ifelse(grepl("YEAST", ev2_04$Protein.Names), "YEAST", 
                                                       ifelse(grepl("ECOLI", ev2_04$Protein.Names), "ECOLI", "unknown"))))
  
  ev2_04 <- ev2_04[!grepl("unknown", ev2_04$species),]
  ev2_04$Quant <- "MS1"
  ev2_04$Label <- "Unlabeled"
  
  return(ev2_04)
  
}



############ LF MS1 quant output minimum peptide intensity

MS1_LF_d0d4_min <- function(ev2_0,ev2_4){
  
  ev2_04 <- ev2_0 %>% inner_join(ev2_4, by =c("seqcharge"="seqcharge_d4"))
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm" = Ms1.Area/median(Ms1.Area)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm" = Ms1.Area_d4/median(Ms1.Area_d4)) %>% ungroup()
  ev2_04 <- ev2_04[which(ev2_04$Ms1.Area & ev2_04$Ms1.Area_d4 > 0),] #keep only non-zero values
  
  ev2_04 <- ev2_04 %>% rowwise() %>% mutate("d0_d4_avg_int" = min(Ms1.Area, Ms1.Area_d4)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d4_avg_int" = median(d0_d4_avg_int)) %>% ungroup()
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
  ev2_04$Label <- "Unlabeled"
  
  return(ev2_04)
  
}

#
MS1_LF_d0d8_min <- function(ev2_0,ev2_4){
  
  ev2_04 <- ev2_0 %>% inner_join(ev2_4, by =c("seqcharge"="seqcharge_d8"))
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm" = Ms1.Area/median(Ms1.Area)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm" = Ms1.Area_d8/median(Ms1.Area_d8)) %>% ungroup()
  ev2_04 <- ev2_04[which(ev2_04$Ms1.Area & ev2_04$Ms1.Area_d8 > 0),] #keep only non-zero values
  
  ev2_04 <- ev2_04 %>% rowwise() %>% mutate("d0_d8_avg_int" = min(Ms1.Area, Ms1.Area_d8)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d8_avg_int" = median(d0_d8_avg_int)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_norm_avg" = mean(d0_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_norm_avg" = mean(d4_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d0_d8_rat_individual" = d0_norm/d4_norm)
  ev2_04 <- ev2_04 %>% mutate("d0_d8_rat" = d0_norm_avg/d4_norm_avg)
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("pep_rat_CV_d0d8" = sd(d0_d8_rat_individual, na.rm=T)/mean(d0_d8_rat_individual)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  

  ev2_04 <- ev2_04%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04$Protein.Names), "HUMAN", 
                                                ifelse(grepl("YEAST", ev2_04$Protein.Names), "YEAST", 
                                                       ifelse(grepl("ECOLI", ev2_04$Protein.Names), "ECOLI", "unknown"))))
  
  ev2_04 <- ev2_04[!grepl("unknown", ev2_04$species),]
  ev2_04$Quant <- "MS1"
  ev2_04$Label <- "Unlabeled"
  
  return(ev2_04)
  
}


#
MS1_LF_d4d8_min <- function(ev2_4,ev2_8){
  
  ev2_4 <- ev2_lim[grepl(paste0(RF2),ev2_lim$Run),] #assumes the 2nd run is delta4 equivalent
  ev2_04 <- ev2_4 %>% inner_join(ev2_8, by =c("seqcharge"="seqcharge_d8"))
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm" = Ms1.Area/median(Ms1.Area)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm" = Ms1.Area_d8/median(Ms1.Area_d8)) %>% ungroup()
  ev2_04 <- ev2_04[which(ev2_04$Ms1.Area & ev2_04$Ms1.Area_d8 > 0),] #keep only non-zero values
  
  ev2_04 <- ev2_04 %>% rowwise() %>% mutate("d4_d8_avg_int" = min(Ms1.Area, Ms1.Area_d8)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_d8_avg_int" = median(d4_d8_avg_int)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_norm_avg" = mean(d0_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_norm_avg" = mean(d4_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d4_d8_rat_individual" = d0_norm/d4_norm)
  ev2_04 <- ev2_04 %>% mutate("d4_d8_rat" = d0_norm_avg/d4_norm_avg)
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("pep_rat_CV_d4d8" = sd(d4_d8_rat_individual, na.rm=T)/mean(d4_d8_rat_individual)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()

  ev2_04 <- ev2_04%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04$Protein.Names), "HUMAN", 
                                                ifelse(grepl("YEAST", ev2_04$Protein.Names), "YEAST", 
                                                       ifelse(grepl("ECOLI", ev2_04$Protein.Names), "ECOLI", "unknown"))))
  
  ev2_04 <- ev2_04[!grepl("unknown", ev2_04$species),]
  ev2_04$Quant <- "MS1"
  ev2_04$Label <- "Unlabeled"
  
  return(ev2_04)
  
}



####################################









################################################ MS1 and MS2 quant. averaged:

MS1MS2_LF_d0d4 <- function(ev2_0,ev2_4){
  
  ev2_04a <- ev2_0 %>% inner_join(ev2_4, by =c("seqcharge_file"="seqcharge_file"))
  #ev2_04 <- ev2_04[which((ev2_04$Precursor.Translated & ev2_04$Precursor.Translated_d4 > 0) | (ev2_04$Ms1.Area & ev2_04$Ms1.Area_d4 > 0)),] #keep only non-zero values
  ev2_04 <- ev2_04a[which((ev2_04a$Precursor.Translated & ev2_04a$Precursor.Translated_d4 > 0)),] #keep only non-zero values
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm_MS2" = Precursor.Translated/median(Precursor.Translated)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm_MS2" = Precursor.Translated_d4/median(Precursor.Translated_d4)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d4_avg_int_MS2" = mean((Precursor.Translated +Precursor.Translated_d4)/2)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_norm_avg_MS2" = mean(d0_norm_MS2)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_norm_avg_MS2" = mean(d4_norm_MS2)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d0_d4_rat_MS2" = d0_norm_avg_MS2/d4_norm_avg_MS2)
  ev2_04_1 <- ev2_04a[which((ev2_04a$Ms1.Area & ev2_04a$Ms1.Area_d4 > 0)),] #keep only non-zero values
  ev2_04_1 <- ev2_04_1 %>% group_by(Run) %>% mutate("d0_norm_MS1" = Ms1.Area/median(Ms1.Area)) %>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% group_by(Run) %>% mutate("d4_norm_MS1" = Ms1.Area_d4/median(Ms1.Area_d4)) %>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% group_by(seqcharge) %>% mutate("d0_d4_avg_int_MS1" = mean((Ms1.Area +Ms1.Area_d4)/2)) %>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% group_by(seqcharge) %>% mutate("d0_norm_avg_MS1" = mean(d0_norm_MS1)) %>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% group_by(seqcharge) %>% mutate("d4_norm_avg_MS1" = mean(d4_norm_MS1)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% mutate("d0_d4_rat_MS1" = d0_norm_avg_MS1/d4_norm_avg_MS1)
  ev2_04_1 <- ev2_04_1 %>% dplyr::select(Protein.Names,seqcharge, d0_d4_avg_int_MS1, d0_d4_rat_MS1)
  
  ev2_04_both <- ev2_04 %>% full_join(ev2_04_1, by =c("seqcharge" = "seqcharge"))
  ev2_04_both$d0_d4_rat <- rowMeans(subset(ev2_04_both, select = c(d0_d4_rat_MS2, d0_d4_rat_MS1)), na.rm = TRUE)
  ev2_04_both$d0_d4_avg_int <- rowMeans(subset(ev2_04_both, select = c(d0_d4_avg_int_MS2, d0_d4_avg_int_MS1)), na.rm = TRUE)
  
  
  ev2_04_both <- ev2_04_both%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04_both$Protein.Names.x) | grepl("HUMAN", ev2_04_both$Protein.Names.y), "HUMAN", 
                                                          ifelse(grepl("YEAST", ev2_04_both$Protein.Names.x) | grepl("YEAST",ev2_04_both$Protein.Names.y), "YEAST", 
                                                                 ifelse(grepl("ECOLI", ev2_04_both$Protein.Names.x) | grepl("ECOLI", ev2_04_both$Protein.Names.y), "ECOLI", "unknown"))))
  
  ev2_04_both <- ev2_04_both[!grepl("unknown", ev2_04_both$species),]
  ev2_04_both$Quant <- "MS2_MS1"
  ev2_04_both$Label <- "Unlabeled"
  
  return(ev2_04_both)
  
}

#d0_d8
MS1MS2_LF_d0d8 <- function(ev2_0,ev2_4){
  
  ev2_04a <- ev2_0 %>% inner_join(ev2_4, by =c("seqcharge_file"="seqcharge_file"))
  #ev2_04 <- ev2_04[which((ev2_04$Precursor.Translated & ev2_04$Precursor.Translated_d8 > 0) | (ev2_04$Ms1.Area & ev2_04$Ms1.Area_d8 > 0)),] #keep only non-zero values
  ev2_04 <- ev2_04a[which((ev2_04a$Precursor.Translated & ev2_04a$Precursor.Translated_d8 > 0)),] #keep only non-zero values
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm_MS2" = Precursor.Translated/median(Precursor.Translated)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d8_norm_MS2" = Precursor.Translated_d8/median(Precursor.Translated_d8)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d8_avg_int_MS2" = mean((Precursor.Translated +Precursor.Translated_d8)/2)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_norm_avg_MS2" = mean(d0_norm_MS2)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d8_norm_avg_MS2" = mean(d8_norm_MS2)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d0_d8_rat_MS2" = d0_norm_avg_MS2/d8_norm_avg_MS2)
  ev2_04_1 <- ev2_04a[which((ev2_04a$Ms1.Area & ev2_04a$Ms1.Area_d8 > 0)),] #keep only non-zero values
  ev2_04_1 <- ev2_04_1 %>% group_by(Run) %>% mutate("d0_norm_MS1" = Ms1.Area/median(Ms1.Area)) %>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% group_by(Run) %>% mutate("d8_norm_MS1" = Ms1.Area_d8/median(Ms1.Area_d8)) %>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% group_by(seqcharge) %>% mutate("d0_d8_avg_int_MS1" = mean((Ms1.Area +Ms1.Area_d8)/2)) %>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% group_by(seqcharge) %>% mutate("d0_norm_avg_MS1" = mean(d0_norm_MS1)) %>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% group_by(seqcharge) %>% mutate("d8_norm_avg_MS1" = mean(d8_norm_MS1)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% mutate("d0_d8_rat_MS1" = d0_norm_avg_MS1/d8_norm_avg_MS1)
  ev2_04_1 <- ev2_04_1 %>% dplyr::select(Protein.Names,seqcharge, d0_d8_avg_int_MS1, d0_d8_rat_MS1)
  
  ev2_04_both <- ev2_04 %>% full_join(ev2_04_1, by =c("seqcharge" = "seqcharge"))
  ev2_04_both$d0_d8_rat <- rowMeans(subset(ev2_04_both, select = c(d0_d8_rat_MS2, d0_d8_rat_MS1)), na.rm = TRUE)
  ev2_04_both$d0_d8_avg_int <- rowMeans(subset(ev2_04_both, select = c(d0_d8_avg_int_MS2, d0_d8_avg_int_MS1)), na.rm = TRUE)
  
  
  ev2_04_both <- ev2_04_both%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04_both$Protein.Names.x) | grepl("HUMAN", ev2_04_both$Protein.Names.y), "HUMAN", 
                                                          ifelse(grepl("YEAST", ev2_04_both$Protein.Names.x) | grepl("YEAST",ev2_04_both$Protein.Names.y), "YEAST", 
                                                                 ifelse(grepl("ECOLI", ev2_04_both$Protein.Names.x) | grepl("ECOLI", ev2_04_both$Protein.Names.y), "ECOLI", "unknown"))))
  
  ev2_04_both <- ev2_04_both[!grepl("unknown", ev2_04_both$species),]
  ev2_04_both$Quant <- "MS2_MS1"
  ev2_04_both$Label <- "Unlabeled"
  
  return(ev2_04_both)
  
}
#d4_d8

MS1MS2_LF_d4d8 <- function(ev2_0,ev2_4){
  
  ev2_4 <- ev2_lim[grepl(paste0(RF1[2]),ev2_lim$Run),] #assumes the 2nd run is delta4 equivalent
  
  ev2_04 <- ev2_4 %>% inner_join(ev2_8, by =c("seqcharge_file"="seqcharge_file"))
  ev2_04 <- ev2_04[which((ev2_04$Precursor.Translated & ev2_04$Precursor.Translated_d8 > 0) | (ev2_04$Ms1.Area & ev2_04$Ms1.Area_d8 > 0)),] #keep only non-zero values
  ev2_04 <- ev2_04[which((ev2_04$Precursor.Translated & ev2_04$Precursor.Translated_d8 > 0)),] #keep only non-zero values
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm_MS2" = Precursor.Translated/median(Precursor.Translated)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d8_norm_MS2" = Precursor.Translated_d8/median(Precursor.Translated_d8)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_d8_avg_int_MS2" = mean((Precursor.Translated +Precursor.Translated_d8)/2)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_norm_avg_MS2" = mean(d4_norm_MS2)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d8_norm_avg_MS2" = mean(d8_norm_MS2)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d4_d8_rat_MS2" = d4_norm_avg_MS2/d8_norm_avg_MS2)
  #ev2_04_1 <- ev2_04a[which((ev2_04a$Ms1.Area & ev2_04a$Ms1.Area_d8 > 0)),] #keep only non-zero values
  ev2_04_1 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm_MS1" = Ms1.Area/median(Ms1.Area)) %>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% group_by(Run) %>% mutate("d8_norm_MS1" = Ms1.Area_d8/median(Ms1.Area_d8)) %>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% group_by(seqcharge) %>% mutate("d4_d8_avg_int_MS1" = mean((Ms1.Area +Ms1.Area_d8)/2)) %>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% group_by(seqcharge) %>% mutate("d4_norm_avg_MS1" = mean(d4_norm_MS1)) %>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% group_by(seqcharge) %>% mutate("d8_norm_avg_MS1" = mean(d8_norm_MS1)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  ev2_04_both <- ev2_04_1 %>% mutate("d4_d8_rat_MS1" = d4_norm_avg_MS1/d8_norm_avg_MS1)
  #ev2_04_1 <- ev2_04_1 %>% dplyr::select(Protein.Names,seqcharge, d4_d8_avg_int_MS1, d4_d8_rat_MS1)
  
  # ev2_04_both <- ev2_04 %>% full_join(ev2_04_1, by =c("seqcharge" = "seqcharge"))
  ev2_04_both$d4_d8_rat <- rowMeans(subset(ev2_04_both, select = c(d4_d8_rat_MS2, d4_d8_rat_MS1)), na.rm = TRUE)
  ev2_04_both$d4_d8_avg_int <- rowMeans(subset(ev2_04_both, select = c(d4_d8_avg_int_MS2, d4_d8_avg_int_MS1)), na.rm = TRUE)
  
  ev2_04_both <- ev2_04_both%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04_both$Protein.Names), "HUMAN", 
                                                          ifelse(grepl("YEAST", ev2_04_both$Protein.Names), "YEAST", 
                                                                 ifelse(grepl("ECOLI", ev2_04_both$Protein.Names), "ECOLI", "unknown"))))
  
  # 
  # ev2_04_both <- ev2_04_both%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04_both$Protein.Names.x) | grepl("HUMAN", ev2_04_both$Protein.Names.y), "HUMAN", 
  #                                                         ifelse(grepl("YEAST", ev2_04_both$Protein.Names.x) | grepl("YEAST",ev2_04_both$Protein.Names.y), "YEAST", 
  #                                                                ifelse(grepl("ECOLI", ev2_04_both$Protein.Names.x) | grepl("ECOLI", ev2_04_both$Protein.Names.y), "ECOLI", "unknown"))))
  # 
  ev2_04_both <- ev2_04_both[!grepl("unknown", ev2_04_both$species),]
  ev2_04_both$Quant <- "MS2_MS1"
  ev2_04_both$Label <- "Unlabeled"
  
  return(ev2_04_both)
  
}


####################



################################################ MS1 quant.
MS1_t_mTRAQ_d0d4 <- function(ev2_0,ev2_4){
  
  ev2_04 <- ev2_0 %>% inner_join(ev2_4, by =c("seqcharge_file"="seqcharge_file"))
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm" = Ms1.Translated/median(Ms1.Translated)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm" = Ms1.Translated_d4/median(Ms1.Translated_d4)) %>% ungroup()
  ev2_04 <- ev2_04[which(ev2_04$Ms1.Translated & ev2_04$Ms1.Translated_d4 > 0),] #keep only non-zero values
  
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d4_avg_int" = mean((Ms1.Translated +Ms1.Translated_d4)/2)) %>% ungroup()
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

#
MS1_t_mTRAQ_d0d8 <- function(ev2_0,ev2_4){
  
  ev2_04 <- ev2_0 %>% inner_join(ev2_4, by =c("seqcharge_file"="seqcharge_file"))
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm" = Ms1.Translated/median(Ms1.Translated)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm" = Ms1.Translated_d8/median(Ms1.Translated_d8)) %>% ungroup()
  ev2_04 <- ev2_04[which(ev2_04$Ms1.Translated & ev2_04$Ms1.Translated_d8 > 0),] #keep only non-zero values
  
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d8_avg_int" = mean((Ms1.Translated +Ms1.Translated_d8)/2)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_norm_avg" = mean(d0_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_norm_avg" = mean(d4_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d0_d8_rat_individual" = d0_norm/d4_norm)
  ev2_04 <- ev2_04 %>% mutate("d0_d8_rat" = d0_norm_avg/d4_norm_avg)
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("pep_rat_CV_d0d8" = sd(d0_d8_rat_individual, na.rm=T)/mean(d0_d8_rat_individual)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  
  ev2_04 <- ev2_04%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04$Protein.Names), "HUMAN", 
                                                ifelse(grepl("YEAST", ev2_04$Protein.Names), "YEAST", 
                                                       ifelse(grepl("ECOLI", ev2_04$Protein.Names), "ECOLI", "unknown"))))
  
  ev2_04 <- ev2_04[!grepl("unknown", ev2_04$species),]
  ev2_04$Quant <- "MS1"
  ev2_04$Label <- "mTRAQ"
  
  return(ev2_04)
  
}

#
MS1_t_mTRAQ_d4d8 <- function(ev2_0,ev2_4){
  
  ev2_4 <- ev2_lim[grepl("mTRAQ4",ev2_lim$Precursor.Id),]
  ev2_04 <- ev2_4 %>% inner_join(ev2_8, by =c("seqcharge_file"="seqcharge_file"))
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm" = Ms1.Translated/median(Ms1.Translated)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm" = Ms1.Translated_d8/median(Ms1.Translated_d8)) %>% ungroup()
  ev2_04 <- ev2_04[which(ev2_04$Ms1.Translated & ev2_04$Ms1.Translated_d8 > 0),] #keep only non-zero values
  
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_d8_avg_int" = mean((Ms1.Translated +Ms1.Translated_d8)/2)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_norm_avg" = mean(d0_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_norm_avg" = mean(d4_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d4_d8_rat_individual" = d0_norm/d4_norm)
  ev2_04 <- ev2_04 %>% mutate("d4_d8_rat" = d0_norm_avg/d4_norm_avg)
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("pep_rat_CV_d4d8" = sd(d4_d8_rat_individual, na.rm=T)/mean(d4_d8_rat_individual)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  
  ev2_04 <- ev2_04%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04$Protein.Names), "HUMAN", 
                                                ifelse(grepl("YEAST", ev2_04$Protein.Names), "YEAST", 
                                                       ifelse(grepl("ECOLI", ev2_04$Protein.Names), "ECOLI", "unknown"))))
  
  ev2_04 <- ev2_04[!grepl("unknown", ev2_04$species),]
  ev2_04$Quant <- "MS1"
  ev2_04$Label <- "mTRAQ"
  
  return(ev2_04)
  
}



######################


################################################ MS1 and MS2 quant. averaged:

MS1pref_mTRAQ_d0d4 <- function(ev2_0,ev2_4){
  
  ev2_04a <- ev2_0 %>% inner_join(ev2_4, by =c("seqcharge_file"="seqcharge_file"))
  ev2_04 <- ev2_04a[which((ev2_04a$Precursor.Translated & ev2_04a$Precursor.Translated_d4 > 0)),] #keep only non-zero values
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm_MS2" = Precursor.Translated/median(Precursor.Translated)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm_MS2" = Precursor.Translated_d4/median(Precursor.Translated_d4)) %>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d0_d4_avg_int_MS2" = (Precursor.Translated +Precursor.Translated_d4)/2) %>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d0_d4_rat_MS2" = d0_norm_MS2/d4_norm_MS2)
  
  ev2_04_1 <- ev2_04a[which((ev2_04a$Ms1.Translated & ev2_04a$Ms1.Translated_d4 > 0)),] #keep only non-zero values
  ev2_04_1 <- ev2_04_1 %>% group_by(Run) %>% mutate("d0_norm_MS1" = Ms1.Translated/median(Ms1.Translated)) %>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% group_by(Run) %>% mutate("d4_norm_MS1" = Ms1.Translated_d4/median(Ms1.Translated_d4)) %>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% mutate("d0_d4_avg_int_MS1" = (Ms1.Translated +Ms1.Translated_d4)/2) %>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% mutate("d0_d4_rat_MS1" = d0_norm_MS1/d4_norm_MS1)
  ev2_04_1 <- ev2_04_1 %>% dplyr::select(seqcharge_file, Protein.Group, Protein.Names, seqcharge, d0_d4_avg_int_MS1, d0_d4_rat_MS1)
  
  ev2_04_both <- ev2_04 %>% full_join(ev2_04_1, by =c("seqcharge_file" = "seqcharge_file"))
  #if MS1 is real, say MS1... If MS1 is NA, then check MS2.. If MS2 is NA then leave it NA, otherwise say MS2.
  ev2_04_both <- ev2_04_both %>% mutate("Quant" = ifelse(is.na(d0_d4_rat_MS1) == FALSE, "MS1", 
                                                         ifelse(is.na(d0_d4_rat_MS2) == TRUE, NA, "MS2"))) 
  
  ev2_04_both <- ev2_04_both %>% mutate("d0_d4_avg_int_temp" = ifelse(grepl("MS1", Quant), d0_d4_avg_int_MS1, 
                                                             ifelse(grepl("MS2", Quant), d0_d4_avg_int_MS2, NA)))
  
  ev2_04_both <- ev2_04_both %>% mutate("d0_d4_rat_temp" = ifelse(grepl("MS1", Quant), d0_d4_rat_MS1, 
                                                                 ifelse(grepl("MS2", Quant), d0_d4_rat_MS2, NA)))
  
  ev2_04_both <- ev2_04_both %>% group_by(seqcharge.x) %>% mutate("d0_d4_rat" = median(d0_d4_rat_temp, na.rm=T)) %>% ungroup()
  ev2_04_both <- ev2_04_both %>% group_by(seqcharge.x) %>% mutate("d0_d4_avg_int" = median(d0_d4_avg_int_temp, na.rm=T)) %>% 
    ungroup() %>% distinct(seqcharge.x, .keep_all=T)
  
  ev2_04_both$Protein.Names <- ifelse(is.na(ev2_04_both$Protein.Names.x), ev2_04_both$Protein.Names.y, ev2_04_both$Protein.Names.x) 
  ev2_04_both$Protein.Group <- ifelse(is.na(ev2_04_both$Protein.Group.x), ev2_04_both$Protein.Group.y, ev2_04_both$Protein.Group.x) 

  ev2_04_both$Label <- "mTRAQ"
  
  ev2_04_both <- ev2_04_both%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04_both$Protein.Names), "HUMAN", 
                                                ifelse(grepl("YEAST", ev2_04_both$Protein.Names), "YEAST", 
                                                       ifelse(grepl("ECOLI", ev2_04_both$Protein.Names), "ECOLI", "unknown"))))
  
  ev2_04_both <- ev2_04_both[!grepl("unknown", ev2_04_both$species),]
  
  return(ev2_04_both)
  
}

#d0_d8
MS1pref_mTRAQ_d0d8 <- function(ev2_0,ev2_4){
  
  ev2_04a <- ev2_0 %>% inner_join(ev2_4, by =c("seqcharge_file"="seqcharge_file"))
  ev2_04 <- ev2_04a[which((ev2_04a$Precursor.Translated & ev2_04a$Precursor.Translated_d8 > 0)),] #keep only non-zero values
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm_MS2" = Precursor.Translated/median(Precursor.Translated)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d8_norm_MS2" = Precursor.Translated_d8/median(Precursor.Translated_d8)) %>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d0_d8_avg_int_MS2" = (Precursor.Translated +Precursor.Translated_d8)/2) %>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d0_d8_rat_MS2" = d0_norm_MS2/d8_norm_MS2)
  
  ev2_04_1 <- ev2_04a[which((ev2_04a$Ms1.Translated & ev2_04a$Ms1.Translated_d8 > 0)),] #keep only non-zero values
  ev2_04_1 <- ev2_04_1 %>% group_by(Run) %>% mutate("d0_norm_MS1" = Ms1.Translated/median(Ms1.Translated)) %>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% group_by(Run) %>% mutate("d8_norm_MS1" = Ms1.Translated_d8/median(Ms1.Translated_d8)) %>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% mutate("d0_d8_avg_int_MS1" = (Ms1.Translated +Ms1.Translated_d8)/2) %>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% mutate("d0_d8_rat_MS1" = d0_norm_MS1/d8_norm_MS1)
  ev2_04_1 <- ev2_04_1 %>% dplyr::select(seqcharge_file, Protein.Group, Protein.Names, seqcharge, d0_d8_avg_int_MS1, d0_d8_rat_MS1)
  
  ev2_04_both <- ev2_04 %>% full_join(ev2_04_1, by =c("seqcharge_file" = "seqcharge_file"))
  #if MS1 is real, say MS1... If MS1 is NA, then check MS2.. If MS2 is NA then leave it NA, otherwise say MS2.
  ev2_04_both <- ev2_04_both %>% mutate("Quant" = ifelse(is.na(d0_d8_rat_MS1) == FALSE, "MS1", 
                                                         ifelse(is.na(d0_d8_rat_MS2) == TRUE, NA, "MS2"))) 
  
  ev2_04_both <- ev2_04_both %>% mutate("d0_d8_avg_int_temp" = ifelse(grepl("MS1", Quant), d0_d8_avg_int_MS1, 
                                                                      ifelse(grepl("MS2", Quant), d0_d8_avg_int_MS2, NA)))
  
  ev2_04_both <- ev2_04_both %>% mutate("d0_d8_rat_temp" = ifelse(grepl("MS1", Quant), d0_d8_rat_MS1, 
                                                                  ifelse(grepl("MS2", Quant), d0_d8_rat_MS2, NA)))
  
  ev2_04_both <- ev2_04_both %>% group_by(seqcharge.x) %>% mutate("d0_d8_rat" = median(d0_d8_rat_temp, na.rm=T)) %>% ungroup()
  ev2_04_both <- ev2_04_both %>% group_by(seqcharge.x) %>% mutate("d0_d8_avg_int" = median(d0_d8_avg_int_temp, na.rm=T)) %>% 
    ungroup() %>% distinct(seqcharge.x, .keep_all=T)
  
  ev2_04_both$Protein.Names <- ifelse(is.na(ev2_04_both$Protein.Names.x), ev2_04_both$Protein.Names.y, ev2_04_both$Protein.Names.x) 
  ev2_04_both$Protein.Group <- ifelse(is.na(ev2_04_both$Protein.Group.x), ev2_04_both$Protein.Group.y, ev2_04_both$Protein.Group.x) 
  
  ev2_04_both$Label <- "mTRAQ"
  
  ev2_04_both <- ev2_04_both%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04_both$Protein.Names), "HUMAN", 
                                                          ifelse(grepl("YEAST", ev2_04_both$Protein.Names), "YEAST", 
                                                                 ifelse(grepl("ECOLI", ev2_04_both$Protein.Names), "ECOLI", "unknown"))))
  
  ev2_04_both <- ev2_04_both[!grepl("unknown", ev2_04_both$species),]
  
  return(ev2_04_both)
  
}
#d4_d8

MS1pref_mTRAQ_d4d8 <- function(ev2_0,ev2_4){
  
  ev2_0 <- ev2_lim[grepl("mTRAQ4",ev2_lim$Precursor.Id),] 
  
  ev2_04a <- ev2_0 %>% inner_join(ev2_4, by =c("seqcharge_file"="seqcharge_file"))
  ev2_04 <- ev2_04a[which((ev2_04a$Precursor.Translated & ev2_04a$Precursor.Translated_d8 > 0)),] #keep only non-zero values
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm_MS2" = Precursor.Translated/median(Precursor.Translated)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d8_norm_MS2" = Precursor.Translated_d8/median(Precursor.Translated_d8)) %>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d4_d8_avg_int_MS2" = (Precursor.Translated +Precursor.Translated_d8)/2) %>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d4_d8_rat_MS2" = d4_norm_MS2/d8_norm_MS2)
  
  ev2_04_1 <- ev2_04a[which((ev2_04a$Ms1.Translated & ev2_04a$Ms1.Translated_d8 > 0)),] #keep only non-zero values
  ev2_04_1 <- ev2_04_1 %>% group_by(Run) %>% mutate("d4_norm_MS1" = Ms1.Translated/median(Ms1.Translated)) %>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% group_by(Run) %>% mutate("d8_norm_MS1" = Ms1.Translated_d8/median(Ms1.Translated_d8)) %>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% mutate("d4_d8_avg_int_MS1" = (Ms1.Translated +Ms1.Translated_d8)/2) %>% ungroup()
  ev2_04_1 <- ev2_04_1 %>% mutate("d4_d8_rat_MS1" = d4_norm_MS1/d8_norm_MS1)
  ev2_04_1 <- ev2_04_1 %>% dplyr::select(seqcharge_file, Protein.Group, Protein.Names, seqcharge, d4_d8_avg_int_MS1, d4_d8_rat_MS1)
  
  ev2_04_both <- ev2_04 %>% full_join(ev2_04_1, by =c("seqcharge_file" = "seqcharge_file"))
  #if MS1 is real, say MS1... If MS1 is NA, then check MS2.. If MS2 is NA then leave it NA, otherwise say MS2.
  ev2_04_both <- ev2_04_both %>% mutate("Quant" = ifelse(is.na(d4_d8_rat_MS1) == FALSE, "MS1", 
                                                         ifelse(is.na(d4_d8_rat_MS2) == TRUE, NA, "MS2"))) 
  
  ev2_04_both <- ev2_04_both %>% mutate("d4_d8_avg_int_temp" = ifelse(grepl("MS1", Quant), d4_d8_avg_int_MS1, 
                                                                      ifelse(grepl("MS2", Quant), d4_d8_avg_int_MS2, NA)))
  
  ev2_04_both <- ev2_04_both %>% mutate("d4_d8_rat_temp" = ifelse(grepl("MS1", Quant), d4_d8_rat_MS1, 
                                                                  ifelse(grepl("MS2", Quant), d4_d8_rat_MS2, NA)))
  
  ev2_04_both <- ev2_04_both %>% group_by(seqcharge.x) %>% mutate("d4_d8_rat" = median(d4_d8_rat_temp, na.rm=T)) %>% ungroup()
  ev2_04_both <- ev2_04_both %>% group_by(seqcharge.x) %>% mutate("d4_d8_avg_int" = median(d4_d8_avg_int_temp, na.rm=T)) %>% 
    ungroup() %>% distinct(seqcharge.x, .keep_all=T)
  
  ev2_04_both$Protein.Names <- ifelse(is.na(ev2_04_both$Protein.Names.x), ev2_04_both$Protein.Names.y, ev2_04_both$Protein.Names.x) 
  ev2_04_both$Protein.Group <- ifelse(is.na(ev2_04_both$Protein.Group.x), ev2_04_both$Protein.Group.y, ev2_04_both$Protein.Group.x) 
  
  ev2_04_both$Label <- "mTRAQ"
  
  ev2_04_both <- ev2_04_both%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04_both$Protein.Names), "HUMAN", 
                                                          ifelse(grepl("YEAST", ev2_04_both$Protein.Names), "YEAST", 
                                                                 ifelse(grepl("ECOLI", ev2_04_both$Protein.Names), "ECOLI", "unknown"))))
  
  ev2_04_both <- ev2_04_both[!grepl("unknown", ev2_04_both$species),]
  
  return(ev2_04_both)
  
}


#############

################################################ MS1 quant.
MS1_t_mTRAQ_d0d4_alt <- function(ev2_0,ev2_4){
  
  ev2_04 <- ev2_0 %>% inner_join(ev2_4, by =c("seqcharge_file"="seqcharge_file"))
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm" = Ms1.Translated/median(Ms1.Translated)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm" = Ms1.Translated_d4/median(Ms1.Translated_d4)) %>% ungroup()
  ev2_04 <- ev2_04[which(ev2_04$Ms1.Translated & ev2_04$Ms1.Translated_d4 > 0),] #keep only non-zero values
  
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d4_avg_int" = mean((Ms1.Translated +Ms1.Translated_d4)/2)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_norm_avg" = mean(d0_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_norm_avg" = mean(d4_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d0_d4_rat_individual" = d0_norm/d4_norm)
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d4_rat" = median(d0_d4_rat_individual)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("pep_rat_CV_d0d4" = sd(d0_d4_rat_individual, na.rm=T)/mean(d0_d4_rat_individual)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  
  ev2_04 <- ev2_04%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04$Protein.Names), "HUMAN", 
                                                ifelse(grepl("YEAST", ev2_04$Protein.Names), "YEAST", 
                                                       ifelse(grepl("ECOLI", ev2_04$Protein.Names), "ECOLI", "unknown"))))
  
  ev2_04 <- ev2_04[!grepl("unknown", ev2_04$species),]
  ev2_04$Quant <- "MS1"
  ev2_04$Label <- "mTRAQ"
  
  return(ev2_04)
  
}

#
MS1_t_mTRAQ_d0d8_alt <- function(ev2_0,ev2_4){
  
  ev2_04 <- ev2_0 %>% inner_join(ev2_4, by =c("seqcharge_file"="seqcharge_file"))
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm" = Ms1.Translated/median(Ms1.Translated)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm" = Ms1.Translated_d8/median(Ms1.Translated_d8)) %>% ungroup()
  ev2_04 <- ev2_04[which(ev2_04$Ms1.Translated & ev2_04$Ms1.Translated_d8 > 0),] #keep only non-zero values
  
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d8_avg_int" = mean((Ms1.Translated +Ms1.Translated_d8)/2)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_norm_avg" = mean(d0_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_norm_avg" = mean(d4_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d0_d8_rat_individual" = d0_norm/d4_norm)
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d8_rat" = median(d0_d8_rat_individual)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("pep_rat_CV_d0d8" = sd(d0_d8_rat_individual, na.rm=T)/mean(d0_d8_rat_individual)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  
  ev2_04 <- ev2_04%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04$Protein.Names), "HUMAN", 
                                                ifelse(grepl("YEAST", ev2_04$Protein.Names), "YEAST", 
                                                       ifelse(grepl("ECOLI", ev2_04$Protein.Names), "ECOLI", "unknown"))))
  
  ev2_04 <- ev2_04[!grepl("unknown", ev2_04$species),]
  ev2_04$Quant <- "MS1"
  ev2_04$Label <- "mTRAQ"
  
  return(ev2_04)
  
}

#
MS1_t_mTRAQ_d4d8_alt <- function(ev2_0,ev2_4){
  
  ev2_4 <- ev2_lim[grepl("mTRAQ4",ev2_lim$Precursor.Id),]
  ev2_04 <- ev2_4 %>% inner_join(ev2_8, by =c("seqcharge_file"="seqcharge_file"))
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm" = Ms1.Translated/median(Ms1.Translated)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm" = Ms1.Translated_d8/median(Ms1.Translated_d8)) %>% ungroup()
  ev2_04 <- ev2_04[which(ev2_04$Ms1.Translated & ev2_04$Ms1.Translated_d8 > 0),] #keep only non-zero values
  
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_d8_avg_int" = mean((Ms1.Translated +Ms1.Translated_d8)/2)) %>% ungroup()
  #ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_norm_avg" = mean(d0_norm)) %>% ungroup()
  #ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_norm_avg" = mean(d4_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d4_d8_rat_individual" = d0_norm/d4_norm)
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_d8_rat" = median(d4_d8_rat_individual)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("pep_rat_CV_d4d8" = sd(d4_d8_rat_individual, na.rm=T)/mean(d4_d8_rat_individual)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  
  ev2_04 <- ev2_04%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04$Protein.Names), "HUMAN", 
                                                ifelse(grepl("YEAST", ev2_04$Protein.Names), "YEAST", 
                                                       ifelse(grepl("ECOLI", ev2_04$Protein.Names), "ECOLI", "unknown"))))
  
  ev2_04 <- ev2_04[!grepl("unknown", ev2_04$species),]
  ev2_04$Quant <- "MS1"
  ev2_04$Label <- "mTRAQ"
  
  return(ev2_04)
  
}


############################################# output minimum peptide raw int:


################################################ MS1 quant.
MS1_t_mTRAQ_d0d4_alt_min <- function(ev2_0,ev2_4){
  
  ev2_04 <- ev2_0 %>% inner_join(ev2_4, by =c("seqcharge_file"="seqcharge_file"))
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm" = Ms1.Translated/median(Ms1.Translated)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm" = Ms1.Translated_d4/median(Ms1.Translated_d4)) %>% ungroup()
  ev2_04 <- ev2_04[which(ev2_04$Ms1.Translated & ev2_04$Ms1.Translated_d4 > 0),] #keep only non-zero values
  
  ev2_04 <- ev2_04 %>% rowwise() %>% mutate("d0_d4_avg_int" = min(Ms1.Translated, Ms1.Translated_d4)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d4_avg_int" = median(d0_d4_avg_int)) %>% ungroup()
  #ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_norm_avg" = mean(d0_norm)) %>% ungroup()
  #ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_norm_avg" = mean(d4_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d0_d4_rat_individual" = d0_norm/d4_norm)
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d4_rat" = median(d0_d4_rat_individual)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("pep_rat_CV_d0d4" = sd(d0_d4_rat_individual, na.rm=T)/mean(d0_d4_rat_individual)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  
  ev2_04 <- ev2_04%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04$Protein.Names), "HUMAN", 
                                                ifelse(grepl("YEAST", ev2_04$Protein.Names), "YEAST", 
                                                       ifelse(grepl("ECOLI", ev2_04$Protein.Names), "ECOLI", "unknown"))))
  
  ev2_04 <- ev2_04[!grepl("unknown", ev2_04$species),]
  ev2_04$Quant <- "MS1"
  ev2_04$Label <- "mTRAQ"
  
  return(ev2_04)
  
}

#
MS1_t_mTRAQ_d0d8_alt_min <- function(ev2_0,ev2_4){
  
  ev2_04 <- ev2_0 %>% inner_join(ev2_4, by =c("seqcharge_file"="seqcharge_file"))
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm" = Ms1.Translated/median(Ms1.Translated)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm" = Ms1.Translated_d8/median(Ms1.Translated_d8)) %>% ungroup()
  ev2_04 <- ev2_04[which(ev2_04$Ms1.Translated & ev2_04$Ms1.Translated_d8 > 0),] #keep only non-zero values
  
  ev2_04 <- ev2_04 %>% rowwise() %>% mutate("d0_d8_avg_int" = min(Ms1.Translated, Ms1.Translated_d8)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d8_avg_int" = median(d0_d8_avg_int)) %>% ungroup()
  #ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_norm_avg" = mean(d0_norm)) %>% ungroup()
  #ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_norm_avg" = mean(d4_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d0_d8_rat_individual" = d0_norm/d4_norm)
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d8_rat" = median(d0_d8_rat_individual)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("pep_rat_CV_d0d8" = sd(d0_d8_rat_individual, na.rm=T)/mean(d0_d8_rat_individual)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  
  ev2_04 <- ev2_04%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04$Protein.Names), "HUMAN", 
                                                ifelse(grepl("YEAST", ev2_04$Protein.Names), "YEAST", 
                                                       ifelse(grepl("ECOLI", ev2_04$Protein.Names), "ECOLI", "unknown"))))
  
  ev2_04 <- ev2_04[!grepl("unknown", ev2_04$species),]
  ev2_04$Quant <- "MS1"
  ev2_04$Label <- "mTRAQ"
  
  return(ev2_04)
  
}

#
MS1_t_mTRAQ_d4d8_alt_min <- function(ev2_0,ev2_4){
  
  ev2_4 <- ev2_lim[grepl("mTRAQ4",ev2_lim$Precursor.Id),]
  ev2_04 <- ev2_4 %>% inner_join(ev2_8, by =c("seqcharge_file"="seqcharge_file"))
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm" = Ms1.Translated/median(Ms1.Translated)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm" = Ms1.Translated_d8/median(Ms1.Translated_d8)) %>% ungroup()
  ev2_04 <- ev2_04[which(ev2_04$Ms1.Translated & ev2_04$Ms1.Translated_d8 > 0),] #keep only non-zero values
  
  ev2_04 <- ev2_04 %>% rowwise() %>% mutate("d4_d8_avg_int" = min(Ms1.Translated, Ms1.Translated_d8)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_d8_avg_int" = median(d4_d8_avg_int)) %>% ungroup()
  #ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_norm_avg" = mean(d0_norm)) %>% ungroup()
  #ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_norm_avg" = mean(d4_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d4_d8_rat_individual" = d0_norm/d4_norm)
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_d8_rat" = median(d4_d8_rat_individual)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("pep_rat_CV_d4d8" = sd(d4_d8_rat_individual, na.rm=T)/mean(d4_d8_rat_individual)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  
  ev2_04 <- ev2_04%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04$Protein.Names), "HUMAN", 
                                                ifelse(grepl("YEAST", ev2_04$Protein.Names), "YEAST", 
                                                       ifelse(grepl("ECOLI", ev2_04$Protein.Names), "ECOLI", "unknown"))))
  
  ev2_04 <- ev2_04[!grepl("unknown", ev2_04$species),]
  ev2_04$Quant <- "MS1"
  ev2_04$Label <- "mTRAQ"
  
  return(ev2_04)
  
}

###########


############################################# output minimum peptide raw int:


################################################ MS1 quant.
MS1_mTRAQ_d0d4_alt_min <- function(ev2_0,ev2_4){
  
  ev2_04 <- ev2_0 %>% inner_join(ev2_4, by =c("seqcharge_file"="seqcharge_file"))
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm" = Ms1.Area/median(Ms1.Area)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm" = Ms1.Area_d4/median(Ms1.Area_d4)) %>% ungroup()
  ev2_04 <- ev2_04[which(ev2_04$Ms1.Area & ev2_04$Ms1.Area_d4 > 0),] #keep only non-zero values
  
  ev2_04 <- ev2_04 %>% rowwise() %>% mutate("d0_d4_avg_int" = min(Ms1.Area, Ms1.Area_d4)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d4_avg_int" = median(d0_d4_avg_int)) %>% ungroup()
  #ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_norm_avg" = mean(d0_norm)) %>% ungroup()
  #ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_norm_avg" = mean(d4_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d0_d4_rat_individual" = d0_norm/d4_norm)
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d4_rat" = median(d0_d4_rat_individual)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("pep_rat_CV_d0d4" = sd(d0_d4_rat_individual, na.rm=T)/mean(d0_d4_rat_individual)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  
  ev2_04 <- ev2_04%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04$Protein.Names), "HUMAN", 
                                                ifelse(grepl("YEAST", ev2_04$Protein.Names), "YEAST", 
                                                       ifelse(grepl("ECOLI", ev2_04$Protein.Names), "ECOLI", "unknown"))))
  
  ev2_04 <- ev2_04[!grepl("unknown", ev2_04$species),]
  ev2_04$Quant <- "MS1"
  ev2_04$Label <- "mTRAQ"
  
  return(ev2_04)
  
}

#
MS1_mTRAQ_d0d8_alt_min <- function(ev2_0,ev2_4){
  
  ev2_04 <- ev2_0 %>% inner_join(ev2_4, by =c("seqcharge_file"="seqcharge_file"))
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm" = Ms1.Area/median(Ms1.Area)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm" = Ms1.Area_d8/median(Ms1.Area_d8)) %>% ungroup()
  ev2_04 <- ev2_04[which(ev2_04$Ms1.Area & ev2_04$Ms1.Area_d8 > 0),] #keep only non-zero values
  
  ev2_04 <- ev2_04 %>% rowwise() %>% mutate("d0_d8_avg_int" = min(Ms1.Area, Ms1.Area_d8)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d8_avg_int" = median(d0_d8_avg_int)) %>% ungroup()
  #ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_norm_avg" = mean(d0_norm)) %>% ungroup()
  #ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_norm_avg" = mean(d4_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d0_d8_rat_individual" = d0_norm/d4_norm)
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d8_rat" = median(d0_d8_rat_individual)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("pep_rat_CV_d0d8" = sd(d0_d8_rat_individual, na.rm=T)/mean(d0_d8_rat_individual)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  
  ev2_04 <- ev2_04%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04$Protein.Names), "HUMAN", 
                                                ifelse(grepl("YEAST", ev2_04$Protein.Names), "YEAST", 
                                                       ifelse(grepl("ECOLI", ev2_04$Protein.Names), "ECOLI", "unknown"))))
  
  ev2_04 <- ev2_04[!grepl("unknown", ev2_04$species),]
  ev2_04$Quant <- "MS1"
  ev2_04$Label <- "mTRAQ"
  
  return(ev2_04)
  
}

#
MS1_mTRAQ_d4d8_alt_min <- function(ev2_0,ev2_4){
  
  ev2_4 <- ev2_lim[grepl("mTRAQ4",ev2_lim$Precursor.Id),]
  ev2_04 <- ev2_4 %>% inner_join(ev2_8, by =c("seqcharge_file"="seqcharge_file"))
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm" = Ms1.Area/median(Ms1.Area)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm" = Ms1.Area_d8/median(Ms1.Area_d8)) %>% ungroup()
  ev2_04 <- ev2_04[which(ev2_04$Ms1.Area & ev2_04$Ms1.Area_d8 > 0),] #keep only non-zero values
  
  ev2_04 <- ev2_04 %>% rowwise() %>% mutate("d4_d8_avg_int" = min(Ms1.Area, Ms1.Area_d8)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_d8_avg_int" = median(d4_d8_avg_int)) %>% ungroup()
  #ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_norm_avg" = mean(d0_norm)) %>% ungroup()
  #ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_norm_avg" = mean(d4_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d4_d8_rat_individual" = d0_norm/d4_norm)
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_d8_rat" = median(d4_d8_rat_individual)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("pep_rat_CV_d4d8" = sd(d4_d8_rat_individual, na.rm=T)/mean(d4_d8_rat_individual)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  
  ev2_04 <- ev2_04%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04$Protein.Names), "HUMAN", 
                                                ifelse(grepl("YEAST", ev2_04$Protein.Names), "YEAST", 
                                                       ifelse(grepl("ECOLI", ev2_04$Protein.Names), "ECOLI", "unknown"))))
  
  ev2_04 <- ev2_04[!grepl("unknown", ev2_04$species),]
  ev2_04$Quant <- "MS1"
  ev2_04$Label <- "mTRAQ"
  
  return(ev2_04)
  
}

#################################### MS2_t_alt_min

MS2_mTRAQ_d0d4_min <- function(ev2_0,ev2_4){
  
  ev2_04 <- ev2_0 %>% inner_join(ev2_4, by =c("seqcharge_file"="seqcharge_file"))
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm" = Precursor.Translated/median(Precursor.Translated)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm" = Precursor.Translated_d4/median(Precursor.Translated_d4)) %>% ungroup()
  ev2_04 <- ev2_04[which(ev2_04$Precursor.Translated & ev2_04$Precursor.Translated_d4 > 0),] #keep only non-zero values
  
  ev2_04 <- ev2_04 %>% rowwise() %>% mutate("d0_d4_avg_int" = min(Precursor.Translated, Precursor.Translated_d4)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d4_avg_int" = median(d0_d4_avg_int)) %>% ungroup()
  #ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_norm_avg" = mean(d0_norm)) %>% ungroup()
  #ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_norm_avg" = mean(d4_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d0_d4_rat_individual" = d0_norm/d4_norm)
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d4_rat" = median(d0_d4_rat_individual)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("pep_rat_CV_d0d4" = sd(d0_d4_rat_individual, na.rm=T)/mean(d0_d4_rat_individual)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  
  ev2_04 <- ev2_04%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04$Protein.Names), "HUMAN", 
                                                ifelse(grepl("YEAST", ev2_04$Protein.Names), "YEAST", 
                                                       ifelse(grepl("ECOLI", ev2_04$Protein.Names), "ECOLI", "unknown"))))
  
  ev2_04 <- ev2_04[!grepl("unknown", ev2_04$species),]
  ev2_04$Quant <- "MS2"
  ev2_04$Label <- "mTRAQ"
  
  return(ev2_04)
  
}

#
MS2_mTRAQ_d0d8_min <- function(ev2_0,ev2_4){
  
  ev2_04 <- ev2_0 %>% inner_join(ev2_4, by =c("seqcharge_file"="seqcharge_file"))
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm" = Precursor.Translated/median(Precursor.Translated)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm" = Precursor.Translated_d8/median(Precursor.Translated_d8)) %>% ungroup()
  ev2_04 <- ev2_04[which(ev2_04$Precursor.Translated & ev2_04$Precursor.Translated_d8 > 0),] #keep only non-zero values
  
  ev2_04 <- ev2_04 %>% rowwise() %>% mutate("d0_d8_avg_int" = min(Precursor.Translated, Precursor.Translated_d8)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d8_avg_int" = median(d0_d8_avg_int)) %>% ungroup()
  #ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_norm_avg" = mean(d0_norm)) %>% ungroup()
  #ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_norm_avg" = mean(d4_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d0_d8_rat_individual" = d0_norm/d4_norm)
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d8_rat" = median(d0_d8_rat_individual)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("pep_rat_CV_d0d8" = sd(d0_d8_rat_individual, na.rm=T)/mean(d0_d8_rat_individual)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  
  ev2_04 <- ev2_04%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04$Protein.Names), "HUMAN", 
                                                ifelse(grepl("YEAST", ev2_04$Protein.Names), "YEAST", 
                                                       ifelse(grepl("ECOLI", ev2_04$Protein.Names), "ECOLI", "unknown"))))
  
  ev2_04 <- ev2_04[!grepl("unknown", ev2_04$species),]
  ev2_04$Quant <- "MS2"
  ev2_04$Label <- "mTRAQ"
  
  return(ev2_04)
  
}


MS2_t_mTRAQ_d4d8_min <- function(ev2_0,ev2_4){
  
  ev2_4 <- ev2_lim[grepl("mTRAQ4",ev2_lim$Precursor.Id),]
  ev2_04 <- ev2_4 %>% inner_join(ev2_8, by =c("seqcharge_file"="seqcharge_file"))
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm" = Precursor.Translated/median(Precursor.Translated)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm" = Precursor.Translated_d8/median(Precursor.Translated_d8)) %>% ungroup()
  ev2_04 <- ev2_04[which(ev2_04$Precursor.Translated & ev2_04$Precursor.Translated_d8 > 0),] #keep only non-zero values
  
  ev2_04 <- ev2_04 %>% rowwise() %>% mutate("d4_d8_avg_int" = min(Precursor.Translated, Precursor.Translated_d8)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_d8_avg_int" = median(d4_d8_avg_int)) %>% ungroup()
  #ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_norm_avg" = mean(d0_norm)) %>% ungroup()
  #ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_norm_avg" = mean(d4_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("d4_d8_rat_individual" = d0_norm/d4_norm)
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_d8_rat" = median(d4_d8_rat_individual)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("pep_rat_CV_d4d8" = sd(d4_d8_rat_individual, na.rm=T)/mean(d4_d8_rat_individual)) %>% distinct(seqcharge,.keep_all=T)%>% ungroup()
  
  ev2_04 <- ev2_04%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04$Protein.Names), "HUMAN", 
                                                ifelse(grepl("YEAST", ev2_04$Protein.Names), "YEAST", 
                                                       ifelse(grepl("ECOLI", ev2_04$Protein.Names), "ECOLI", "unknown"))))
  
  ev2_04 <- ev2_04[!grepl("unknown", ev2_04$species),]
  ev2_04$Quant <- "MS1"
  ev2_04$Label <- "mTRAQ"
  
  return(ev2_04)
  
}


########################## data filtering:

dat_PG_Prep <- function(ev) {
  ev$seqcharge <- paste0(ev$Modified.Sequence, ev$Precursor.Charge)
  ev$seqcharge <- str_remove_all(ev$seqcharge, "\\(mTRAQ0\\)")
  ev$seqcharge <- str_remove_all(ev$seqcharge, "\\(mTRAQ4\\)")
  ev$seqcharge <- str_remove_all(ev$seqcharge, "\\(mTRAQ8\\)")
  ev$seqcharge_file <- paste0(ev$seqcharge,ev$Run)
  
  ev <- ev[which(ev$PG.Q.Value < 0.01),]
  length(unique(ev$Protein.Group))
  
  ev <- ev %>% mutate("label" = ifelse(grepl("mTRAQ0", Modified.Sequence), "mTRAQ0",
                                       ifelse(grepl("mTRAQ4", Modified.Sequence), "mTRAQ4", "mTRAQ8")))
  
  ev2_lim <- ev
  
  ev2_0 <- ev2_lim[grepl("mTRAQ0",ev2_lim$Precursor.Id),]
  ev2_4 <- ev2_lim[grepl("mTRAQ4",ev2_lim$Precursor.Id),] %>% dplyr::select("seqcharge","seqcharge_file","Precursor.Quantity","Precursor.Normalised","Precursor.Translated","Ms1.Translated")
  colnames(ev2_4) <- c("seqcharge_d4", "seqcharge_file","Precursor.Quantity_d4","Precursor.Normalised_d4","Precursor.Translated_d4","Ms1.Translated_d4")
  ev2_8 <- ev2_lim[grepl("mTRAQ8",ev2_lim$Precursor.Id),] %>% dplyr::select("seqcharge","seqcharge_file","Precursor.Quantity","Precursor.Normalised","Precursor.Translated","Ms1.Translated")
  colnames(ev2_8) <- c("seqcharge_d8", "seqcharge_file","Precursor.Quantity_d8","Precursor.Normalised_d8","Precursor.Translated_d8","Ms1.Translated_d8")
  
  dfs <- list(ev2_0,ev2_4,ev2_8, ev2_lim)
  return(dfs)
}


###########  MS2 alt

MS2_mTRAQ_d0d4_min_alt <- function(ev2_0,ev2_4){
  
  ev2_04 <- ev2_0 %>% inner_join(ev2_4, by =c("seqcharge_file"="seqcharge_file"))
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d0_norm" = Precursor.Translated/median(Precursor.Translated)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Run) %>% mutate("d4_norm" = Precursor.Translated_d4/median(Precursor.Translated_d4)) %>% ungroup()
  ev2_04 <- ev2_04[which(ev2_04$Precursor.Translated & ev2_04$Precursor.Translated_d4 > 0),] #keep only non-zero values
  
  ev2_04 <- ev2_04 %>% rowwise() %>% mutate("d0_d4_avg_int" = min(Precursor.Translated, Precursor.Translated_d4)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_d4_avg_int" = median(d0_d4_avg_int)) %>% ungroup()
  ev2_04 <- ev2_04 %>% mutate("PG_run" = paste0(Protein.Group,Run))
  #ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d0_norm_avg" = mean(d0_norm)) %>% ungroup()
  #ev2_04 <- ev2_04 %>% group_by(seqcharge) %>% mutate("d4_norm_avg" = mean(d4_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(PG_run) %>% mutate("med_PG_d0" = median(d0_norm)) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(PG_run) %>% mutate("med_PG_d4" = median(d4_norm)) %>% ungroup() %>% distinct(PG_run, .keep_all = T)
  ev2_04 <- ev2_04 %>% group_by(PG_run) %>% mutate("PG_d0_d4_rat_run" = med_PG_d0/med_PG_d4) %>% ungroup()
  ev2_04 <- ev2_04 %>% group_by(Protein.Group) %>% mutate("PG_d0_d4_rat" = median(med_PG_d0)/median(med_PG_d4)) %>% ungroup()
  

  ev2_04 <- ev2_04%>% mutate("species" = ifelse(grepl("HUMAN", ev2_04$Protein.Names), "HUMAN", 
                                                ifelse(grepl("YEAST", ev2_04$Protein.Names), "YEAST", 
                                                       ifelse(grepl("ECOLI", ev2_04$Protein.Names), "ECOLI", "unknown"))))
  
  ev2_04 <- ev2_04[!grepl("unknown", ev2_04$species),]
  ev2_04$Quant <- "MS1"
  ev2_04$Label <- "mTRAQ"
  
  return(ev2_04)
  
}




