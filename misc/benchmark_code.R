#!/usr/bin/env Rscript

#---------------------------------------------
#Code used for benchmark
#---------------------------------------------

rm(list=ls())
library(dplyr)
library(tidyr)
library(stringdist)
library(ggplot2)
library(ggsci)
library(Cairo)
library(jsonlite)
library(RColorBrewer)

options(bitmapType="cairo")

#Generate data----
N_FOLD_NEG_VAL<-5
K_NEG<-50
MAX_EPOCH<-70
SEED<-2000
FILE_DIR<-paste("feb20_pseudoVgene_upHc2Mouse","_val",N_FOLD_NEG_VAL,"neg_train",K_NEG,"neg_maxEpoch",MAX_EPOCH,sep = "")

set.seed(seed = SEED)

pairing_data_old<-read.table("cleaned_data_combined_completed.csv",sep=",",header = T)
pairing_data_new<-read.table("new_cleaned_data_combined.csv",sep=",",header = T)
pairing_data<-rbind(pairing_data_old[,-1],pairing_data_new)
rm(list=c("pairing_data_new","pairing_data_old"))

pairing_data$CDR3_BETA[pairing_data$CDR3_BETA=="`CASSGDGMNTEAFF"]<-"CASSGDGMNTEAFF"
pairing_data<-pairing_data[,-11]
pairing_data<-pairing_data[!duplicated(pairing_data),]
pairing_data$TRA_V[pairing_data$TRA_V==""]<-"XXX"
pairing_data$CDR3_ALPHA[pairing_data$CDR3_ALPHA==""]<-"XXX"
pairing_data$TRB_V[pairing_data$TRB_V==""]<-"XXX"
pairing_data$CDR3_BETA[pairing_data$CDR3_BETA==""]<-"XXX"
pairing_data$TRA_V[pairing_data$TRA_V=="TRAV14D-1*01" & pairing_data$TCR_SPECIES=="human"]<-"TRAV14/DV4*01"
pairing_data$MHC_ALPHA[pairing_data$MHC_ALPHA=="" & pairing_data$MHC_BETA=="HLA-DPB1*02:01"]<-"HLA-DPA1*01:03"
pairing_data<-cbind(pairing_data,data.frame(mhca1=gsub("HLA-","",pairing_data$MHC_ALPHA)))
pairing_data<-cbind(pairing_data,data.frame(mhca2=gsub("H2","H-2",pairing_data$mhca1)))
pairing_data<-cbind(pairing_data,data.frame(mhca3=pairing_data$mhca2))
pairing_data$mhca3[grepl(":",pairing_data$mhca3,fixed = T)]<-sapply(strsplit(pairing_data$mhca3[grepl(":",pairing_data$mhca3,fixed = T)],split = ":"),function(x) paste(x[1],x[2],sep=":"))

pairing_data<-cbind(pairing_data,data.frame(mhcb1=gsub("human_glublin","human_microglobulin",pairing_data$MHC_BETA)))
pairing_data<-cbind(pairing_data,data.frame(mhcb2=gsub("mouse_glublin","mouse_microglobulin",pairing_data$mhcb1)))
pairing_data<-cbind(pairing_data,data.frame(mhcb3=gsub("HLA-","",pairing_data$mhcb2)))
pairing_data<-cbind(pairing_data,data.frame(mhcb4=gsub("H2","H-2",pairing_data$mhcb3)))
pairing_data<-cbind(pairing_data,data.frame(mhcb5=pairing_data$mhcb4))
pairing_data$mhcb5[grepl(":",pairing_data$mhcb5,fixed = T)]<-sapply(strsplit(pairing_data$mhcb5[grepl(":",pairing_data$mhcb5,fixed = T)],split = ":"),function(x) paste(x[1],x[2],sep=":"))
pairing_data$mhca3[pairing_data$mhca3=="DPA1*01:01"]<-"DPA1*01:03"
pairing_data$mhca3[pairing_data$mhca3=="A*24:01"]<-"A*24:02"
pairing_data$mhca3[pairing_data$mhca3=="B*07:01"]<-"B*07:02"
pairing_data$mhca3[pairing_data$mhca3=="B*44:01"]<-"B*44:02"
pairing_data$mhca3[pairing_data$mhca3=="B*12:01"]<-"B*44:02"

#substitute pseudo V gene
pairing_data$TRA_V[pairing_data$TRA_V=="TRAV4-1*01" & pairing_data$TCR_SPECIES=="mouse"]<-"XXX"
pairing_data$TRA_V[pairing_data$TRA_V=="TRAV22*01" & pairing_data$TCR_SPECIES=="mouse"]<-"XXX"
pairing_data$TRA_V[pairing_data$TRA_V=="TRAV20*01" & pairing_data$TCR_SPECIES=="mouse"]<-"XXX"
pairing_data$TRB_V[pairing_data$TRB_V=="TRBV8*01" & pairing_data$TCR_SPECIES=="mouse"]<-"XXX"
pairing_data$TRB_V[pairing_data$TRB_V=="TRBV27*01" & pairing_data$TCR_SPECIES=="mouse"]<-"XXX"
pairing_data$TRB_V[pairing_data$TRB_V=="TRBV28*01" & pairing_data$TCR_SPECIES=="mouse"]<-"XXX"
pairing_data$TRB_V[pairing_data$TRB_V=="TRBV9*01" & pairing_data$TCR_SPECIES=="mouse"]<-"XXX"
pairing_data$TRB_V[pairing_data$TRB_V=="TRBV18*01" & pairing_data$TCR_SPECIES=="mouse"]<-"XXX"
pairing_data$TRB_V[pairing_data$TRB_V=="TRBV12-3*01" & pairing_data$TCR_SPECIES=="mouse"]<-"XXX"
pairing_data$TRA_V[pairing_data$TRA_V=="TRAV28*01" & pairing_data$TCR_SPECIES=="human"]<-"XXX"
pairing_data$TRA_V[pairing_data$TRA_V=="TRAV15*01" & pairing_data$TCR_SPECIES=="human"]<-"XXX"
pairing_data$TRA_V[pairing_data$TRA_V=="TRAV14-1*01" & pairing_data$TCR_SPECIES=="human"]<-"XXX"
pairing_data$TRA_V[pairing_data$TRA_V=="TRAV8-5*01" & pairing_data$TCR_SPECIES=="human"]<-"XXX"
pairing_data$TRA_V[pairing_data$TRA_V=="TRAV33*01" & pairing_data$TCR_SPECIES=="human"]<-"XXX"
pairing_data$TRB_V[pairing_data$TRB_V=="TRBV12-2*01" & pairing_data$TCR_SPECIES=="human"]<-"XXX"
pairing_data$TRB_V[pairing_data$TRB_V=="TRBV8-1*01" & pairing_data$TCR_SPECIES=="human"]<-"XXX"
pairing_data$TRB_V[pairing_data$TRB_V=="TRBV3-2*01" & pairing_data$TCR_SPECIES=="human"]<-"XXX"
pairing_data$TRB_V[pairing_data$TRB_V=="TRBV3-2*02" & pairing_data$TCR_SPECIES=="human"]<-"XXX"
pairing_data$TRB_V[pairing_data$TRB_V=="TRBV8-2*01" & pairing_data$TCR_SPECIES=="human"]<-"XXX"
pairing_data$TRB_V[pairing_data$TRB_V=="TRBV5-2*01" & pairing_data$TCR_SPECIES=="human"]<-"XXX"
pairing_data$TRB_V[pairing_data$TRB_V=="TRBV22-1*01" & pairing_data$TCR_SPECIES=="human"]<-"XXX"
pairing_data$TRB_V[pairing_data$TRB_V=="TRBV12-1*01" & pairing_data$TCR_SPECIES=="human"]<-"XXX"

pairing_data<-pairing_data[,c(1,2,3,4,5,6,13,18,9,10)]
pairing_data<-pairing_data[!duplicated(pairing_data),]
pairing_data<-pairing_data[!(grepl("XXX",pairing_data$CDR3_ALPHA) & grepl("XXX",pairing_data$CDR3_BETA)),]

#add v gene sequence
mouse_a<-Biostrings::readAAStringSet(filepath = "mouse_trav_pro.txt.completed.txt")
mouse_name_a <- names(mouse_a)
mouse_seq_a <- paste(mouse_a)
mouse_vgene_a<-data.frame(a_allele=sapply(strsplit(mouse_name_a,split = "|",fixed = T),function(x){x[2]}), mouse_seq_a)
mouse_vgene_a<-rbind(mouse_vgene_a,data.frame(a_allele="XXX",mouse_seq_a="XXX"))

mouse_b<-Biostrings::readAAStringSet(filepath = "mouse_trbv_pro.txt.completed.txt")
mouse_name_b <- names(mouse_b)
mouse_seq_b <- paste(mouse_b)
mouse_vgene_b<-data.frame(b_allele=sapply(strsplit(mouse_name_b,split = "|",fixed = T),function(x){x[2]}), mouse_seq_b)
mouse_vgene_b<-rbind(mouse_vgene_b,data.frame(b_allele="XXX",mouse_seq_b="XXX"))

human_a<-Biostrings::readAAStringSet(filepath = "human_trav_pro.txt.completed.txt")
human_name_a <- names(human_a)
human_seq_a <- paste(human_a)
vgene_human_a<-data.frame(a_allele=sapply(strsplit(human_name_a,split = "|",fixed = T),function(x){x[2]}), human_seq_a)
vgene_human_a<-rbind(vgene_human_a,data.frame(a_allele="XXX",human_seq_a="XXX"))

human_b<-Biostrings::readAAStringSet(filepath = "human_trbv_pro.txt.completed.txt")
human_name_b <- names(human_b)
human_seq_b <- paste(human_b)
vgene_human_b<-data.frame(b_allele=sapply(strsplit(human_name_b,split = "|",fixed = T),function(x){x[2]}), human_seq_b)
vgene_human_b<-rbind(vgene_human_b,data.frame(b_allele="XXX",human_seq_b="XXX"))

pairing_data_h<-pairing_data[pairing_data$TCR_SPECIES=="human" & (pairing_data$TRA_V %in% vgene_human_a$a_allele) & (pairing_data$TRB_V %in% vgene_human_b$b_allele),]
pairing_data_m<-pairing_data[pairing_data$TCR_SPECIES=="mouse" & (pairing_data$TRA_V %in% mouse_vgene_a$a_allele) & (pairing_data$TRB_V %in% mouse_vgene_b$b_allele),]
pairing_data_h<-merge(pairing_data_h,vgene_human_a,by.x="TRA_V",by.y="a_allele",all.x=T)
pairing_data_h<-merge(pairing_data_h,vgene_human_b,by.x="TRB_V",by.y="b_allele",all.x=T)
pairing_data_m<-merge(pairing_data_m,mouse_vgene_a,by.x="TRA_V",by.y="a_allele",all.x=T)
pairing_data_m<-merge(pairing_data_m,mouse_vgene_b,by.x="TRB_V",by.y="b_allele",all.x=T)
stopifnot(sum(is.na(pairing_data_h))==0)
stopifnot(sum(is.na(pairing_data_m))==0)
pairing_data_h<-subset(pairing_data_h,select=c(data_file,CDR3_ALPHA,CDR3_BETA,EPITOPE,mhca3,mhcb5,pMHC_SPECIES,TCR_SPECIES,human_seq_a,human_seq_b))
pairing_data_m<-subset(pairing_data_m,select=c(data_file,CDR3_ALPHA,CDR3_BETA,EPITOPE,mhca3,mhcb5,pMHC_SPECIES,TCR_SPECIES,mouse_seq_a,mouse_seq_b))
colnames(pairing_data_h)<-c("data_file","CDR3_ALPHA","CDR3_BETA","EPITOPE","mhca3","mhcb5","pMHC_SPECIES","TCR_SPECIES","TRA_V","TRB_V")
colnames(pairing_data_m)<-c("data_file","CDR3_ALPHA","CDR3_BETA","EPITOPE","mhca3","mhcb5","pMHC_SPECIES","TCR_SPECIES","TRA_V","TRB_V")
pairing_data<-rbind(pairing_data_h,pairing_data_m)
stopifnot(sum(is.na(pairing_data))==0)

# deal with missing sequences
pairing_data_complete<-pairing_data[!grepl("XXX",pairing_data$TRA_V) & !grepl("XXX",pairing_data$CDR3_ALPHA) & !grepl("XXX", pairing_data$TRB_V) & !grepl("XXX",pairing_data$CDR3_BETA),]
pairing_data_mis_va<-pairing_data[grepl("XXX",pairing_data$TRA_V) & !grepl("XXX",pairing_data$CDR3_ALPHA) & !grepl("XXX", pairing_data$TRB_V) & !grepl("XXX",pairing_data$CDR3_BETA),]
pairing_data_mis_vb<-pairing_data[!grepl("XXX",pairing_data$TRA_V) & !grepl("XXX",pairing_data$CDR3_ALPHA) & grepl("XXX", pairing_data$TRB_V) & !grepl("XXX",pairing_data$CDR3_BETA),]
pairing_data_mis_cdr3a<-pairing_data[!grepl("XXX",pairing_data$TRA_V) & grepl("XXX",pairing_data$CDR3_ALPHA) & !grepl("XXX", pairing_data$TRB_V) & !grepl("XXX",pairing_data$CDR3_BETA),]
pairing_data_mis_cdr3b<-pairing_data[!grepl("XXX",pairing_data$TRA_V) & !grepl("XXX",pairing_data$CDR3_ALPHA) & !grepl("XXX", pairing_data$TRB_V) & grepl("XXX",pairing_data$CDR3_BETA),]
pairing_data_vavb<-pairing_data[grepl("XXX",pairing_data$TRA_V) & !grepl("XXX",pairing_data$CDR3_ALPHA) & grepl("XXX", pairing_data$TRB_V) & !grepl("XXX",pairing_data$CDR3_BETA),]
pairing_data_vacdr3a<-pairing_data[grepl("XXX",pairing_data$TRA_V) & grepl("XXX",pairing_data$CDR3_ALPHA) & !grepl("XXX", pairing_data$TRB_V) & !grepl("XXX",pairing_data$CDR3_BETA),]
pairing_data_vacdr3b<-pairing_data[grepl("XXX",pairing_data$TRA_V) & !grepl("XXX",pairing_data$CDR3_ALPHA) & !grepl("XXX", pairing_data$TRB_V) & grepl("XXX",pairing_data$CDR3_BETA),]
pairing_data_vbcdr3a<-pairing_data[!grepl("XXX",pairing_data$TRA_V) & grepl("XXX",pairing_data$CDR3_ALPHA) & grepl("XXX", pairing_data$TRB_V) & !grepl("XXX",pairing_data$CDR3_BETA),]
pairing_data_vbcdr3b<-pairing_data[!grepl("XXX",pairing_data$TRA_V) & !grepl("XXX",pairing_data$CDR3_ALPHA) & grepl("XXX", pairing_data$TRB_V) & grepl("XXX",pairing_data$CDR3_BETA),]
pairing_data_vavbcdr3a<-pairing_data[grepl("XXX",pairing_data$TRA_V) & grepl("XXX",pairing_data$CDR3_ALPHA) & grepl("XXX", pairing_data$TRB_V) & !grepl("XXX",pairing_data$CDR3_BETA),]
pairing_data_vavbcdr3b<-pairing_data[grepl("XXX",pairing_data$TRA_V) & !grepl("XXX",pairing_data$CDR3_ALPHA) & grepl("XXX", pairing_data$TRB_V) & grepl("XXX",pairing_data$CDR3_BETA),]

#validation positive
val_pos_dataset<-c("Axelrod","Balko","Braunlein","Brunk","Ciacchi","Dahal","Foy","Hanada","hT27","Huisman","Luo2","Moore","Orphan","p53R175H","Peri","Poncette","Sooda","Stadinski","Ting","Ueno","10X2",
                   "ATLAS","Burrows","Cobo","DeWitt","Lichterfeld","Luo3","Ma","Mallajosyula","Nguyen","Ogunshola","Rowntree","Schinkelshoek","Shimizu","Sim","TetTCR_SeqHD_validate",
                   "Uchida","Watson","Wu","Yu","Assmus","Lu","NeoScreen","Wagner","Chen","Grant","Joglekar","Motozono","Wahl","TCR3d")
stopifnot(sum(unique(val_pos_dataset) %in% pairing_data_complete$data_file)==50)
validation<-pairing_data_complete[pairing_data_complete$data_file %in% val_pos_dataset,]
validation<-cbind(validation,label=rep(1))

#validation negative
h_neg_back_a<-read.table("human_alpha.txt",header = T,sep="\t",quote = "")
h_neg_back_b<-read.table("human_beta.txt",header = T,sep="\t",quote = "")
m_neg_back_a<-read.table("mouse_alpha.txt",header = T,sep="\t",quote = "") 
m_neg_back_b<-read.table("mouse_beta.txt",header = T,sep="\t",quote = "")
validation<-cbind(validation,valrandn=sample.int(n = nrow(validation)*100, size = nrow(validation)))
val_tcrh<-validation[validation$TCR_SPECIES=="human",]
val_tcrm<-validation[validation$TCR_SPECIES=="mouse",]
val_tcrh_tmp1<-val_tcrh
val_tcrh_tmp1<-val_tcrh_tmp1[,c(1,4,5,6,7,8,12,9,2,10,3,11)]
colnames(val_tcrh_tmp1)<-c("data_file","EPITOPE","mhca3","mhcb5","pMHC_SPECIES","TCR_SPECIES","valrandn","vaseq","cdr3a","vbseq","cdr3b","label")
val_tcrm_tmp1<-val_tcrm
val_tcrm_tmp1<-val_tcrm_tmp1[,c(1,4,5,6,7,8,12,9,2,10,3,11)]
colnames(val_tcrm_tmp1)<-c("data_file","EPITOPE","mhca3","mhcb5","pMHC_SPECIES","TCR_SPECIES","valrandn","vaseq","cdr3a","vbseq","cdr3b","label")

for (i in 1:(N_FOLD_NEG_VAL*2)) {
    val_tcrh_tmp1<-rbind(val_tcrh_tmp1, cbind(val_tcrh[,c("data_file","EPITOPE","mhca3","mhcb5","pMHC_SPECIES","TCR_SPECIES","valrandn")],
                                            h_neg_back_a[sample(x = 1:nrow(h_neg_back_a),size = nrow(val_tcrh)),c(3,2)],
                                            h_neg_back_b[sample(x = 1:nrow(h_neg_back_b),size = nrow(val_tcrh)),c(3,2)],
                                            label=rep(0)))
  val_tcrm_tmp1<-rbind(val_tcrm_tmp1, cbind(val_tcrm[,c("data_file","EPITOPE","mhca3","mhcb5","pMHC_SPECIES","TCR_SPECIES","valrandn")],
                                            m_neg_back_a[sample(x = 1:nrow(m_neg_back_a),size = nrow(val_tcrm)),c(3,2)],
                                            m_neg_back_b[sample(x = 1:nrow(m_neg_back_b),size = nrow(val_tcrm)),c(3,2)],
                                            label=rep(0)))
}

final_validation<-rbind(val_tcrh_tmp1,val_tcrm_tmp1)
colnames(final_validation)<-c("data_file","peptide","mhca","mhcb","pMHC_SPECIES","TCR_SPECIES","valrandn","va","cdr3a","vb","cdr3b","label")
tyu<-data.frame(pmhc=paste(final_validation$peptide,final_validation$mhca,final_validation$mhcb,final_validation$pMHC_SPECIES,sep = "_"),
                tcr=paste(final_validation$cdr3a,final_validation$cdr3b,final_validation$TCR_SPECIES,sep = "_"),
                label=final_validation$label)
stopifnot(sum(duplicated(tyu[tyu$label==0,]))==0)
stopifnot(length(intersect(tyu[tyu$label==1,'tcr'],tyu[tyu$label==0,'tcr']))==0)
stopifnot(sum(is.na(final_validation))==0)
stopifnot(sum(final_validation=="")==0)

#remove duplicates
fin_val_tmp<-subset(final_validation,final_validation$label==0,-c(label))
paring_tmp<-subset(pairing_data,select = -c(data_file))
colnames(paring_tmp)<-c("cdr3a","cdr3b","peptide","mhca","mhcb","pMHC_SPECIES","TCR_SPECIES","va","vb")
fin_val_tmp2<-dplyr::anti_join(fin_val_tmp,paring_tmp)

colnames(fin_val_tmp2)<-c("data_file","EPITOPE","mhca3","mhcb5","pMHC_SPECIES","TCR_SPECIES","valrandn","TRA_V","CDR3_ALPHA","TRB_V","CDR3_BETA")
fin_val_tmp2<-dplyr::anti_join(fin_val_tmp2,pairing_data_complete,by=c("EPITOPE","mhca3","mhcb5","pMHC_SPECIES","TCR_SPECIES","TRA_V","CDR3_ALPHA","TRB_V","CDR3_BETA"))
fin_val_tmp2<-dplyr::anti_join(fin_val_tmp2,pairing_data_mis_va,by=c("EPITOPE","mhca3","mhcb5","pMHC_SPECIES","TCR_SPECIES","CDR3_ALPHA","TRB_V","CDR3_BETA"))
fin_val_tmp2<-dplyr::anti_join(fin_val_tmp2,pairing_data_mis_vb,by=c("EPITOPE","mhca3","mhcb5","pMHC_SPECIES","TCR_SPECIES","TRA_V","CDR3_ALPHA","CDR3_BETA"))
fin_val_tmp2<-dplyr::anti_join(fin_val_tmp2,pairing_data_mis_cdr3a,by=c("EPITOPE","mhca3","mhcb5","pMHC_SPECIES","TCR_SPECIES","TRA_V","TRB_V","CDR3_BETA"))
fin_val_tmp2<-dplyr::anti_join(fin_val_tmp2,pairing_data_mis_cdr3b,by=c("EPITOPE","mhca3","mhcb5","pMHC_SPECIES","TCR_SPECIES","TRA_V","CDR3_ALPHA","TRB_V"))
fin_val_tmp2<-dplyr::anti_join(fin_val_tmp2,pairing_data_vavb,by=c("EPITOPE","mhca3","mhcb5","pMHC_SPECIES","TCR_SPECIES","CDR3_ALPHA","CDR3_BETA"))
fin_val_tmp2<-dplyr::anti_join(fin_val_tmp2,pairing_data_vacdr3a,by=c("EPITOPE","mhca3","mhcb5","pMHC_SPECIES","TCR_SPECIES","TRB_V","CDR3_BETA"))
fin_val_tmp2<-dplyr::anti_join(fin_val_tmp2,pairing_data_vacdr3b,by=c("EPITOPE","mhca3","mhcb5","pMHC_SPECIES","TCR_SPECIES","CDR3_ALPHA","TRB_V"))
fin_val_tmp2<-dplyr::anti_join(fin_val_tmp2,pairing_data_vbcdr3a,by=c("EPITOPE","mhca3","mhcb5","pMHC_SPECIES","TCR_SPECIES","TRA_V","CDR3_BETA"))
fin_val_tmp2<-dplyr::anti_join(fin_val_tmp2,pairing_data_vbcdr3b,by=c("EPITOPE","mhca3","mhcb5","pMHC_SPECIES","TCR_SPECIES","TRA_V","CDR3_ALPHA"))
fin_val_tmp2<-dplyr::anti_join(fin_val_tmp2,pairing_data_vavbcdr3a,by=c("EPITOPE","mhca3","mhcb5","pMHC_SPECIES","TCR_SPECIES","CDR3_BETA"))
fin_val_tmp2<-dplyr::anti_join(fin_val_tmp2,pairing_data_vavbcdr3b,by=c("EPITOPE","mhca3","mhcb5","pMHC_SPECIES","TCR_SPECIES","CDR3_ALPHA"))
colnames(fin_val_tmp2)<-c('data_file','peptide','mhca','mhcb','pMHC_SPECIES','TCR_SPECIES','valrandn','va','cdr3a','vb','cdr3b')
final_val_seq<-rbind(final_validation[final_validation$label==1,],cbind(fin_val_tmp2,label=rep(0)))

#format
final_val_seq_pos<-final_val_seq[final_val_seq$label==1,]
final_val_seq_neg<-final_val_seq[final_val_seq$label==0,]
final_val_seq_neg_sampled<-final_val_seq_neg %>% group_by(valrandn) %>% sample_n(N_FOLD_NEG_VAL)
final_val_seq2<-rbind(final_val_seq_pos,final_val_seq_neg_sampled)
final_val_seq3<-final_val_seq2 %>% group_by(valrandn) %>% sample_n(1+N_FOLD_NEG_VAL)
stopifnot(all(dim(final_val_seq2)==dim(final_val_seq3)))

#write data
final_validation<-as.data.frame(final_val_seq3)
final_validation<-subset(final_validation,select = c("data_file","peptide","mhca","mhcb","pMHC_SPECIES","TCR_SPECIES","va","cdr3a","vb","cdr3b","label"))
colnames(final_val_seq3)<-c("data_file","peptide","mhca","mhcb","pMHC_SPECIES","TCR_SPECIES","valrandn","vaseq","cdr3a","vbseq","cdr3b","label")
write.table(x = final_val_seq3,file = paste(FILE_DIR,"/validation.txt",sep=""),sep="\t",quote = F,row.names = F,col.names = T)
rm(list=setdiff(ls(), c('FILE_DIR',"final_validation","h_neg_back_a","h_neg_back_b","m_neg_back_a","m_neg_back_b",
                        "K_NEG","MAX_EPOCH","pairing_data","pairing_data_complete","pairing_data_mis_va","pairing_data_mis_vb","pairing_data_mis_cdr3a","pairing_data_mis_cdr3b",
                        "pairing_data_vavb","pairing_data_vacdr3a","pairing_data_vacdr3b","pairing_data_vbcdr3a","pairing_data_vbcdr3b","pairing_data_vavbcdr3a","pairing_data_vavbcdr3b")))

#positive training
final_val_copy<-final_validation[,c(7,8,9,10,2,3,4,5,6)]
colnames(final_val_copy)<-c("TRA_V","CDR3_ALPHA","TRB_V","CDR3_BETA","EPITOPE","mhca3","mhcb5","pMHC_SPECIES","TCR_SPECIES")
train_complete<-pairing_data_complete[!(pairing_data_complete$data_file %in% final_validation$data_file),]
train_complete<-dplyr::anti_join(train_complete[,-1],final_val_copy)
train_complete<-train_complete[!duplicated(train_complete),]

pairing_data_mis_va<-pairing_data_mis_va[!(pairing_data_mis_va$data_file %in% final_validation$data_file),]
pairing_data_mis_vb<-pairing_data_mis_vb[!(pairing_data_mis_vb$data_file %in% final_validation$data_file),]
pairing_data_mis_cdr3a<-pairing_data_mis_cdr3a[!(pairing_data_mis_cdr3a$data_file %in% final_validation$data_file),]
pairing_data_mis_cdr3b<-pairing_data_mis_cdr3b[!(pairing_data_mis_cdr3b$data_file %in% final_validation$data_file),]
pairing_data_vavb<-pairing_data_vavb[!(pairing_data_vavb$data_file %in% final_validation$data_file),]
pairing_data_vacdr3a<-pairing_data_vacdr3a[!(pairing_data_vacdr3a$data_file %in% final_validation$data_file),]
pairing_data_vacdr3b<-pairing_data_vacdr3b[!(pairing_data_vacdr3b$data_file %in% final_validation$data_file),]
pairing_data_vbcdr3a<-pairing_data_vbcdr3a[!(pairing_data_vbcdr3a$data_file %in% final_validation$data_file),]
pairing_data_vbcdr3b<-pairing_data_vbcdr3b[!(pairing_data_vbcdr3b$data_file %in% final_validation$data_file),]
pairing_data_vavbcdr3a<-pairing_data_vavbcdr3a[!(pairing_data_vavbcdr3a$data_file %in% final_validation$data_file),]
pairing_data_vavbcdr3b<-pairing_data_vavbcdr3b[!(pairing_data_vavbcdr3b$data_file %in% final_validation$data_file),]

pairing_data_mis_va<-pairing_data_mis_va[,-1]
pairing_data_mis_vb<-pairing_data_mis_vb[,-1]
pairing_data_mis_cdr3a<-pairing_data_mis_cdr3a[,-1]
pairing_data_mis_cdr3b<-pairing_data_mis_cdr3b[,-1]
pairing_data_vavb<-pairing_data_vavb[,-1]
pairing_data_vacdr3a<-pairing_data_vacdr3a[,-1]
pairing_data_vacdr3b<-pairing_data_vacdr3b[,-1]
pairing_data_vbcdr3a<-pairing_data_vbcdr3a[,-1]
pairing_data_vbcdr3b<-pairing_data_vbcdr3b[,-1]
pairing_data_vavbcdr3a<-pairing_data_vavbcdr3a[,-1]
pairing_data_vavbcdr3b<-pairing_data_vavbcdr3b[,-1]

#deal with missing value
val_va<-final_val_copy
val_va$TRA_V<-"XXX"
pairing_data_mis_va<-dplyr::anti_join(pairing_data_mis_va,val_va)

val_vb<-final_val_copy
val_vb$TRB_V<-"XXX"
pairing_data_mis_vb<-dplyr::anti_join(pairing_data_mis_vb,val_vb)

val_cdr3a<-final_val_copy
val_cdr3a$CDR3_ALPHA<-"XXX"
pairing_data_mis_cdr3a<-dplyr::anti_join(pairing_data_mis_cdr3a,val_cdr3a)

val_cdr3b<-final_val_copy
val_cdr3b$CDR3_BETA<-"XXX"
pairing_data_mis_cdr3b<-dplyr::anti_join(pairing_data_mis_cdr3b,val_cdr3b)

val_vavb<-final_val_copy
val_vavb$TRA_V<-"XXX"
val_vavb$TRB_V<-"XXX"
pairing_data_vavb<-dplyr::anti_join(pairing_data_vavb,val_vavb)

val_vacdr3a<-final_val_copy
val_vacdr3a$TRA_V<-"XXX"
val_vacdr3a$CDR3_ALPHA<-"XXX"
pairing_data_vacdr3a<-dplyr::anti_join(pairing_data_vacdr3a,val_vacdr3a)

val_vacdr3b<-final_val_copy
val_vacdr3b$TRA_V<-"XXX"
val_vacdr3b$CDR3_BETA<-"XXX"
pairing_data_vacdr3b<-dplyr::anti_join(pairing_data_vacdr3b,val_vacdr3b)

val_vbcdr3a<-final_val_copy
val_vbcdr3a$TRB_V<-"XXX"
val_vbcdr3a$CDR3_ALPHA<-"XXX"
pairing_data_vbcdr3a<-dplyr::anti_join(pairing_data_vbcdr3a,val_vbcdr3a)

val_vbcdr3b<-final_val_copy
val_vbcdr3b$TRB_V<-"XXX"
val_vbcdr3b$CDR3_BETA<-"XXX"
pairing_data_vbcdr3b<-dplyr::anti_join(pairing_data_vbcdr3b,val_vbcdr3b)

val_vavbcdr3a<-final_val_copy
val_vavbcdr3a$TRA_V<-"XXX"
val_vavbcdr3a$TRB_V<-"XXX"
val_vavbcdr3a$CDR3_ALPHA<-"XXX"
pairing_data_vavbcdr3a<-dplyr::anti_join(pairing_data_vavbcdr3a,val_vavbcdr3a)

val_vavbcdr3b<-final_val_copy
val_vavbcdr3b$TRA_V<-"XXX"
val_vavbcdr3b$TRB_V<-"XXX"
val_vavbcdr3b$CDR3_BETA<-"XXX"
pairing_data_vavbcdr3b<-dplyr::anti_join(pairing_data_vavbcdr3b,val_vavbcdr3b)

train_complete<-cbind(train_complete,class=rep(1))
pairing_data_mis_va<-cbind(pairing_data_mis_va,class=rep(2))
pairing_data_mis_cdr3a<-cbind(pairing_data_mis_cdr3a,class=rep(3))
pairing_data_mis_vb<-cbind(pairing_data_mis_vb,class=rep(4))
pairing_data_mis_cdr3b<-cbind(pairing_data_mis_cdr3b,class=rep(5))
pairing_data_vacdr3a<-cbind(pairing_data_vacdr3a,class=rep(6))
pairing_data_vavb<-cbind(pairing_data_vavb,class=rep(7))
pairing_data_vacdr3b<-cbind(pairing_data_vacdr3b,class=rep(8))
pairing_data_vbcdr3a<-cbind(pairing_data_vbcdr3a,class=rep(9))
pairing_data_vbcdr3b<-cbind(pairing_data_vbcdr3b,class=rep(10))
pairing_data_vavbcdr3a<-cbind(pairing_data_vavbcdr3a,class=rep(11))
pairing_data_vavbcdr3b<-cbind(pairing_data_vavbcdr3b,class=rep(12))

training_pos<-Reduce(rbind,list(train_complete,
                                pairing_data_mis_va,pairing_data_mis_vb,pairing_data_mis_cdr3a,pairing_data_mis_cdr3b,
                                pairing_data_vavb,pairing_data_vacdr3a,pairing_data_vacdr3b,pairing_data_vbcdr3a,pairing_data_vbcdr3b,
                                pairing_data_vavbcdr3a,pairing_data_vavbcdr3b))

training_pos<-training_pos[!duplicated(training_pos),]

stopifnot(sum(duplicated(training_pos))==0)
stopifnot(sum(is.na(training_pos))==0)
stopifnot(sum(training_pos=="")==0)
training_pos<-cbind(training_pos,label=rep(1))

#upsampling
PMHC_MOUSE_UPSAMPLING_FOLD<-10
PMHC_HC2_UPSAMPLING_FOLD<-10
training_pos_mouse_tmp<-training_pos[training_pos$pMHC_SPECIES=="mouse",]
training_pos_hc2_tmp<-training_pos[training_pos$pMHC_SPECIES=="human" & training_pos$mhcb5!="human_microglobulin",]

training_pos<-rbind(training_pos,training_pos_mouse_tmp[rep(seq_len(nrow(training_pos_mouse_tmp)), PMHC_MOUSE_UPSAMPLING_FOLD),])
training_pos<-rbind(training_pos,training_pos_hc2_tmp[rep(seq_len(nrow(training_pos_hc2_tmp)), PMHC_HC2_UPSAMPLING_FOLD),])

rownames(training_pos)<-NULL
training_pos_ori<-training_pos

rm(list=c("training_pos_mouse_tmp","training_pos_hc2_tmp"))
rm(list=setdiff(ls(), c('FILE_DIR',"final_validation","h_neg_back_a","h_neg_back_b","m_neg_back_a","m_neg_back_b",
                        "K_NEG","MAX_EPOCH","training_pos_ori")))

for (SEED in seq(1, MAX_EPOCH, 1)) {
  set.seed(seed = SEED)
  training_pos<-training_pos_ori
  
  training_pos<-cbind(training_pos,randn=paste(SEED,sample.int(n = nrow(training_pos)*100, size = nrow(training_pos)),sep = "_"))
  training_pos<-training_pos[sample(x = 1:nrow(training_pos),size = nrow(training_pos)),]
  training_pos_n<-training_pos[rep(seq_len(nrow(training_pos)), (K_NEG+150)),]
  
  rownames(training_pos_n)<-NULL

  training_pos_n<-subset(training_pos_n,select = -c(label))
  training_pos_n_h<-training_pos_n[training_pos_n$TCR_SPECIES=="human",]
  training_pos_n_m<-training_pos_n[training_pos_n$TCR_SPECIES=="mouse",]
  
  training_neg_tmp_h<-cbind(training_pos_n_h,
                            h_neg_back_a[sample(x = 1:nrow(h_neg_back_a),size = nrow(training_pos_n_h),replace = T),c(3,2)],
                            h_neg_back_b[sample(x = 1:nrow(h_neg_back_b),size = nrow(training_pos_n_h),replace = T),c(3,2)],
                            label=rep(0))
  
  training_neg_tmp_m<-cbind(training_pos_n_m,
                            m_neg_back_a[sample(x = 1:nrow(m_neg_back_a),size = nrow(training_pos_n_m),replace = T),c(3,2)],
                            m_neg_back_b[sample(x = 1:nrow(m_neg_back_b),size = nrow(training_pos_n_m),replace = T),c(3,2)],
                            label=rep(0))
  
  training_neg_tmp<-rbind(training_neg_tmp_h,training_neg_tmp_m)
  
  colnames(training_neg_tmp)<-c("CDR3_ALPHA","CDR3_BETA","EPITOPE","mhca3","mhcb5","pMHC_SPECIES","TCR_SPECIES","TRA_V","TRB_V","class","randn","va","cdr3a","vb","cdr3b","label")
  
  training_neg_tmp[training_neg_tmp$class==2,"va"]<-"XXX"
  training_neg_tmp[training_neg_tmp$class==3,"cdr3a"]<-"XXX"
  training_neg_tmp[training_neg_tmp$class==4,"vb"]<-"XXX"
  training_neg_tmp[training_neg_tmp$class==5,"cdr3b"]<-"XXX"
  training_neg_tmp[training_neg_tmp$class==6,"va"]<-"XXX"
  training_neg_tmp[training_neg_tmp$class==6,"cdr3a"]<-"XXX"
  training_neg_tmp[training_neg_tmp$class==7,"va"]<-"XXX"
  training_neg_tmp[training_neg_tmp$class==7,"vb"]<-"XXX"
  training_neg_tmp[training_neg_tmp$class==8,"va"]<-"XXX"
  training_neg_tmp[training_neg_tmp$class==8,"cdr3b"]<-"XXX"
  training_neg_tmp[training_neg_tmp$class==9,"vb"]<-"XXX"
  training_neg_tmp[training_neg_tmp$class==9,"cdr3a"]<-"XXX"
  training_neg_tmp[training_neg_tmp$class==10,"vb"]<-"XXX"
  training_neg_tmp[training_neg_tmp$class==10,"cdr3b"]<-"XXX"
  training_neg_tmp[training_neg_tmp$class==11,"va"]<-"XXX"
  training_neg_tmp[training_neg_tmp$class==11,"vb"]<-"XXX"
  training_neg_tmp[training_neg_tmp$class==11,"cdr3a"]<-"XXX"
  training_neg_tmp[training_neg_tmp$class==12,"va"]<-"XXX"
  training_neg_tmp[training_neg_tmp$class==12,"vb"]<-"XXX"
  training_neg_tmp[training_neg_tmp$class==12,"cdr3b"]<-"XXX"
  
  training_neg_tmp<-training_neg_tmp[!duplicated(subset(training_neg_tmp,select=c(EPITOPE,mhca3,mhcb5,pMHC_SPECIES,TCR_SPECIES,va,cdr3a,vb,cdr3b))),]
  
  # remove duplicates
  final_val_copy<-final_validation[,c(2,3,4,5,6,7,8,9,10)]
  colnames(final_val_copy)<-c("EPITOPE","mhca3","mhcb5","pMHC_SPECIES","TCR_SPECIES","va","cdr3a","vb","cdr3b")
  training_neg_tmp<-dplyr::anti_join(training_neg_tmp,final_val_copy)

  val_va<-final_val_copy
  val_va$va<-"XXX"
  training_neg_tmp<-dplyr::anti_join(training_neg_tmp,val_va)
  
  val_vb<-final_val_copy
  val_vb$vb<-"XXX"
  training_neg_tmp<-dplyr::anti_join(training_neg_tmp,val_vb)
  
  val_cdr3a<-final_val_copy
  val_cdr3a$cdr3a<-"XXX"
  training_neg_tmp<-dplyr::anti_join(training_neg_tmp,val_cdr3a)
  
  val_cdr3b<-final_val_copy
  val_cdr3b$cdr3b<-"XXX"
  training_neg_tmp<-dplyr::anti_join(training_neg_tmp,val_cdr3b)
  
  val_vavb<-final_val_copy
  val_vavb$va<-"XXX"
  val_vavb$vb<-"XXX"
  training_neg_tmp<-dplyr::anti_join(training_neg_tmp,val_vavb)
  
  val_vacdr3a<-final_val_copy
  val_vacdr3a$va<-"XXX"
  val_vacdr3a$cdr3a<-"XXX"
  training_neg_tmp<-dplyr::anti_join(training_neg_tmp,val_vacdr3a)
  
  val_vacdr3b<-final_val_copy
  val_vacdr3b$va<-"XXX"
  val_vacdr3b$cdr3b<-"XXX"
  training_neg_tmp<-dplyr::anti_join(training_neg_tmp,val_vacdr3b)
  
  val_vbcdr3a<-final_val_copy
  val_vbcdr3a$vb<-"XXX"
  val_vbcdr3a$cdr3a<-"XXX"
  training_neg_tmp<-dplyr::anti_join(training_neg_tmp,val_vbcdr3a)
  
  val_vbcdr3b<-final_val_copy
  val_vbcdr3b$vb<-"XXX"
  val_vbcdr3b$cdr3b<-"XXX"
  training_neg_tmp<-dplyr::anti_join(training_neg_tmp,val_vbcdr3b)
  
  val_vavbcdr3a<-final_val_copy
  val_vavbcdr3a$va<-"XXX"
  val_vavbcdr3a$vb<-"XXX"
  val_vavbcdr3a$cdr3a<-"XXX"
  training_neg_tmp<-dplyr::anti_join(training_neg_tmp,val_vavbcdr3a)
  
  val_vavbcdr3b<-final_val_copy
  val_vavbcdr3b$va<-"XXX"
  val_vavbcdr3b$vb<-"XXX"
  val_vavbcdr3b$cdr3b<-"XXX"
  training_neg_tmp<-dplyr::anti_join(training_neg_tmp,val_vavbcdr3b)
  
  training_pos_copy<-training_pos[,1:9]
  colnames(training_pos_copy)<-c("cdr3a","cdr3b","EPITOPE","mhca3","mhcb5","pMHC_SPECIES","TCR_SPECIES","va","vb")
  
  training_neg_tmp2_1<-dplyr::anti_join(training_neg_tmp[training_neg_tmp$class==1,],training_pos_copy)
  
  training_neg_tmp2_2<-dplyr::anti_join(training_neg_tmp[training_neg_tmp$class==2,],subset(training_pos_copy,select= -c(va)))
  training_neg_tmp2_3<-dplyr::anti_join(training_neg_tmp[training_neg_tmp$class==3,],subset(training_pos_copy,select= -c(cdr3a)))
  training_neg_tmp2_4<-dplyr::anti_join(training_neg_tmp[training_neg_tmp$class==4,],subset(training_pos_copy,select= -c(vb)))
  training_neg_tmp2_5<-dplyr::anti_join(training_neg_tmp[training_neg_tmp$class==5,],subset(training_pos_copy,select= -c(cdr3b)))
  
  training_neg_tmp2_6<-dplyr::anti_join(training_neg_tmp[training_neg_tmp$class==6,],subset(training_pos_copy,select= -c(va,cdr3a)))
  training_neg_tmp2_7<-dplyr::anti_join(training_neg_tmp[training_neg_tmp$class==7,],subset(training_pos_copy,select= -c(va,vb)))
  training_neg_tmp2_8<-dplyr::anti_join(training_neg_tmp[training_neg_tmp$class==8,],subset(training_pos_copy,select= -c(va,cdr3b)))
  training_neg_tmp2_9<-dplyr::anti_join(training_neg_tmp[training_neg_tmp$class==9,],subset(training_pos_copy,select= -c(vb,cdr3a)))
  training_neg_tmp2_10<-dplyr::anti_join(training_neg_tmp[training_neg_tmp$class==10,],subset(training_pos_copy,select= -c(vb,cdr3b)))
  
  training_neg_tmp2_11<-dplyr::anti_join(training_neg_tmp[training_neg_tmp$class==11,],subset(training_pos_copy,select= -c(va,vb,cdr3a)))
  training_neg_tmp2_12<-dplyr::anti_join(training_neg_tmp[training_neg_tmp$class==12,],subset(training_pos_copy,select= -c(va,vb,cdr3b)))
  
  training_neg<-Reduce(rbind,list(training_neg_tmp2_1,training_neg_tmp2_2,training_neg_tmp2_3,training_neg_tmp2_4,training_neg_tmp2_5,training_neg_tmp2_6,
                                  training_neg_tmp2_7,training_neg_tmp2_8,training_neg_tmp2_9,training_neg_tmp2_10,training_neg_tmp2_11,training_neg_tmp2_12))

  training_tmp1<-rbind(cbind(training_pos[,-11],va=training_pos$TRA_V,cdr3a=training_pos$CDR3_ALPHA,vb=training_pos$TRB_V,cdr3b=training_pos$CDR3_BETA,label=rep(1)),training_neg)
  
  training_tmp2<-subset(training_tmp1,select = -c(TRA_V,CDR3_ALPHA,TRB_V,CDR3_BETA))
  
  colnames(training_tmp2)<-c("peptide","mhca","mhcb","pMHC_SPECIES","TCR_SPECIES","class","randn","vaseq","cdr3a","vbseq","cdr3b","label")

  stopifnot(sum(is.na(training_tmp2))==0)
  
  #format
  training_tmp2_pos<-training_tmp2[training_tmp2$label==1,]
  training_tmp2_neg<-training_tmp2[training_tmp2$label==0,]
  training_tmp2_neg_sampled<-training_tmp2_neg %>% group_by(randn) %>% sample_n(K_NEG)
  training_tmp3<-rbind(training_tmp2_pos,training_tmp2_neg_sampled)
  training_tmp4<-training_tmp3 %>% group_by(randn) %>% sample_n(1+K_NEG)
  
  write.table(x = training_tmp4,file = paste("/project/",FILE_DIR,"/data",ceiling(SEED/5),"/training_individual_",SEED,".txt",sep=""),sep="\t",quote = F,row.names = F,col.names = T)

  rm(list=setdiff(ls(), c('FILE_DIR',"final_validation","h_neg_back_a","h_neg_back_b","m_neg_back_a","m_neg_back_b","K_NEG","MAX_EPOCH","training_pos_ori")))
}

#Fig.1
library(ggplot2)
library(ggsci)
library(Cairo)
library(jsonlite)
options(bitmapType="cairo")
# prepare input
human_a<-Biostrings::readAAStringSet(filepath = "human_trav_pro.txt.completed.txt")
human_name_a <- names(human_a)
human_seq_a <- paste(human_a)
vgene_human_a<-data.frame(a_allele=sapply(strsplit(human_name_a,split = "|",fixed = T),function(x){x[2]}), human_seq_a)
vgene_human_a<-rbind(vgene_human_a,data.frame(a_allele="XXX",human_seq_a="XXX"))
rownames(vgene_human_a)<-vgene_human_a$a_allele

human_b<-Biostrings::readAAStringSet(filepath = "human_trbv_pro.txt.completed.txt")
human_name_b <- names(human_b)
human_seq_b <- paste(human_b)
vgene_human_b<-data.frame(b_allele=sapply(strsplit(human_name_b,split = "|",fixed = T),function(x){x[2]}), human_seq_b)
vgene_human_b<-rbind(vgene_human_b,data.frame(b_allele="XXX",human_seq_b="XXX"))
rownames(vgene_human_b)<-vgene_human_b$b_allele

input_all<-data.frame()
for (pdir in list.dirs('BEAM-T',recursive = F)) {
  folder_name<-paste(paste(unlist(strsplit(pdir,'[/_]'))[c(11,12,13)],collapse = '_'),tail(unlist(strsplit(pdir,'[/_]')),n=1),sep = '_')
  k5.t<-read.table(
    list.files(path = pdir,pattern = 'annotations.csv',full.names = T,recursive = T),
    header = T,sep=',')
  k5.t<-k5.t[(k5.t$full_length=='true' | k5.t$full_length=='TRUE') & (k5.t$productive=='true' | k5.t$productive=='TRUE'),]
  k5.t$v_gene<-paste(k5.t$v_gene,'*01',sep='')
  stopifnot(sum(!(k5.t$v_gene %in% c(vgene_human_a$a_allele,vgene_human_b$b_allele)))==0)
  k5.t.a<-k5.t[k5.t$chain=='TRA',]
  k5.t.b<-k5.t[k5.t$chain=='TRB',]
  
  #only keep cell with both a and b chains
  k5.t.a<-k5.t.a[k5.t.a$clonotype_id %in% k5.t.b$clonotype_id,]
  k5.t.b<-k5.t.b[k5.t.b$clonotype_id %in% k5.t.a$clonotype_id,]
  
  k5.t.ab<-merge(k5.t.a[,c("clonotype_id", "consensus_id", "length","chain","v_gene","cdr3","umis")],
                 k5.t.b[,c("clonotype_id", "consensus_id", "length","chain","v_gene","cdr3","umis")],
                 by = 'clonotype_id',all=T)
  stopifnot(sum(is.na(k5.t.ab))==0)
  stopifnot(sum(is.null(k5.t.ab))==0)
  
  colnames(k5.t.ab)<-gsub('v_gene.x','va',colnames(k5.t.ab))
  colnames(k5.t.ab)<-gsub('v_gene.y','vb',colnames(k5.t.ab))
  colnames(k5.t.ab)<-gsub('cdr3.x','cdr3a',colnames(k5.t.ab))
  colnames(k5.t.ab)<-gsub('cdr3.y','cdr3b',colnames(k5.t.ab))
  
  k5.t.ab$vaseq<-vgene_human_a[k5.t.ab$va,'human_seq_a']
  k5.t.ab$vbseq<-vgene_human_b[k5.t.ab$vb,'human_seq_b']
  stopifnot(sum(is.na(k5.t.ab))==0)
  stopifnot(sum(is.null(k5.t.ab))==0)
  k5.t.ab$clonotype_id2<-paste(k5.t.ab$clonotype_id,ifelse(duplicated(k5.t.ab$clonotype_id),"_2",""),sep='')
  
  if(max(table(k5.t.ab$clonotype_id))>2){
    print(folder_name)
    print(table(k5.t.ab$clonotype_id)[table(k5.t.ab$clonotype_id)>2])
    k5.t.ab<-k5.t.ab[!(k5.t.ab$clonotype_id %in% names(table(k5.t.ab$clonotype_id)[table(k5.t.ab$clonotype_id)>2])),]
  }
  
  # antigen
  pmhc<-data.frame(peptide=c('LYVDSLFFL','IYNGKLFDL','SYQKVIELF','SYIGLKGLYF','VYVPLKELL','VYALPLKML','QYDPVAALF','AYAQKIFKI'),
                   mhca=rep('A*24:02'),
                   mhcb=rep('human_microglobulin'),
                   pMHC_SPECIES=rep('human'),
                   TCR_SPECIES=rep('human'))
  # make input file
  k5.input<-merge(k5.t.ab,pmhc,all = T)
  k5.input$batch<-rep(folder_name)
  input_all<-rbind(input_all,k5.input)
  rm(list = setdiff(ls(),c('input_all','vgene_human_a','vgene_human_b')))
}
write.table(input_all,'input_pmtnet_15folders.csv',append = F,quote = F,sep = ',',row.names = F,col.names = T)

#antigen specificity data sorting,use 10X's clonotype
all_antigen_specificity<-data.frame()

for (pdir in list.dirs('BEAM-T',recursive = F)) {
  folder_name<-paste(paste(unlist(strsplit(pdir,'[/_]'))[c(11,12,13)],collapse = '_'),tail(unlist(strsplit(pdir,'[/_]')),n=1),sep = '_')
  k5.antigen<-read.table(
    list.files(path = pdir,pattern = 'antigen_specificity_scores.csv',full.names = T,recursive = T),
    header = T,sep=',')
  k5.antigen<-k5.antigen[k5.antigen$raw_clonotype_id!='None',]
  k5.antigen<-k5.antigen[k5.antigen$antigen %in% c("BEAM05", "BEAM06", "BEAM07", "BEAM08", "BEAM09", "BEAM10", "BEAM11", "BEAM12", "BEAM13"),]
  antigen_tmp<-reshape2::dcast(k5.antigen[,c('raw_clonotype_id','antigen',"antigen_specificity_score")],raw_clonotype_id~antigen,mean)
  antigen_tmp$batch<-rep(folder_name)
  all_antigen_specificity<-rbind(all_antigen_specificity,antigen_tmp) 
}
write.table(all_antigen_specificity,'summary_antigen_specificity.csv',quote = F,sep = ',',row.names = F,col.names = T)

antigen.spec<-read.table('summary_antigen_specificity.csv',header = T,sep = ',')

# read clonal frequency
clonal_freq_all<-data.frame()
for (pdir in list.dirs('BEAM-T',recursive = F)) {
  folder_name<-paste(paste(unlist(strsplit(pdir,'[/_]'))[c(11,12,13)],collapse = '_'),tail(unlist(strsplit(pdir,'[/_]')),n=1),sep = '_')
  # TCR
  k5.t<-read.table(
    list.files(path = pdir,pattern = 'clonotypes.csv',full.names = T,recursive = T),
    header = T,sep=',')
  clonal_freq_all<-rbind(clonal_freq_all,cbind(k5.t[,c("clonotype_id","frequency","proportion")],batch=rep(folder_name)))
}
clonal_freq_all$batch_clon<-paste(clonal_freq_all$batch,clonal_freq_all$clonotype_id,sep='_')
rownames(clonal_freq_all)<-clonal_freq_all$batch_clon
rm(list=c('k5.t','folder_name','pdir'))

freq_cutoff<-5

#read results of pMTnet v2 
outv2<-read.table('results_pmtnet_15folders_v2.csv',header = T,sep = ',')
outv2$beam<-NA
outv2$beam[outv2$peptide=='LYVDSLFFL']<-'BEAM05'
outv2$beam[outv2$peptide=='IYNGKLFDL']<-'BEAM06'
outv2$beam[outv2$peptide=='SYQKVIELF']<-'BEAM08'
outv2$beam[outv2$peptide=='SYIGLKGLYF']<-'BEAM09'
outv2$beam[outv2$peptide=='VYVPLKELL']<-'BEAM10'
outv2$beam[outv2$peptide=='VYALPLKML']<-'BEAM11'
outv2$beam[outv2$peptide=='QYDPVAALF']<-'BEAM12'
outv2$beam[outv2$peptide=='AYAQKIFKI']<-'BEAM13'
stopifnot(sum(is.na(outv2))==0)

outv2$bat_clon_beam<-paste(outv2$batch,outv2$clonotype_id,outv2$beam,sep='_')
outv2<-outv2[order(outv2$bat_clon_beam,-outv2$umis.x,-outv2$umis.y),]
outv2<-outv2[!duplicated(outv2$bat_clon_beam),]
outv2$tcr<-paste(outv2$va,outv2$cdr3a,outv2$vb,outv2$cdr3b,sep='_')

outv2$batch_clon<-paste(outv2$batch,outv2$clonotype_id,sep='_')
outv2$freq<-NA
outv2$freq<-clonal_freq_all[outv2$batch_clon,'frequency']
stopifnot(sum(is.na(outv2))==0)
outv2<-outv2[outv2$freq>freq_cutoff,]

for (i in sort(unique(outv2$batch))) {
  out.tmp<-outv2[outv2$batch==i,]
  rownames(out.tmp)<-NULL
  out.tmp.summ<-reshape(out.tmp[,c('beam','avg_externak_rank2','clonotype_id2')],idvar = 'beam',timevar = 'clonotype_id2',direction = 'wide')
  rownames(out.tmp.summ)<-out.tmp.summ$beam
  out.tmp.summ<-out.tmp.summ[,-1]
  colnames(out.tmp.summ)<-unlist(lapply(strsplit(colnames(out.tmp.summ),'.',fixed = T), '[[',2))
  pheatmap::pheatmap(t(as.matrix(log10(out.tmp.summ+1e-30))),main=i)
}

outv2.plot<-outv2[,c('beam','avg_externak_rank2','clonotype_id2')]
outv2.plot$clonotype_id2<-paste(outv2$batch,outv2$clonotype_id2,sep='_')
out.tmp.summ<-reshape(outv2.plot[,c('beam','avg_externak_rank2','clonotype_id2')],idvar = 'beam',timevar = 'clonotype_id2',direction = 'wide')
rownames(out.tmp.summ)<-out.tmp.summ$beam
out.tmp.summ<-t(out.tmp.summ[,-1])

#read results of pMTnet v1 
outv1<-read.table('input_pmtnet_15folders.csv',header = T,sep=',')
pv1<-read.table('pMTnet_20240803164601__prediction.csv',header = T,sep = ',')

outv1$v1pred<-pv1$Rank
outv1$beam<-NA
outv1$beam[outv1$peptide=='LYVDSLFFL']<-'BEAM05'
outv1$beam[outv1$peptide=='IYNGKLFDL']<-'BEAM06'
outv1$beam[outv1$peptide=='SYQKVIELF']<-'BEAM08'
outv1$beam[outv1$peptide=='SYIGLKGLYF']<-'BEAM09'
outv1$beam[outv1$peptide=='VYVPLKELL']<-'BEAM10'
outv1$beam[outv1$peptide=='VYALPLKML']<-'BEAM11'
outv1$beam[outv1$peptide=='QYDPVAALF']<-'BEAM12'
outv1$beam[outv1$peptide=='AYAQKIFKI']<-'BEAM13'
stopifnot(sum(is.na(outv1))==0)

outv1$bat_clon_beam<-paste(outv1$batch,outv1$clonotype_id,outv1$beam,sep='_')
outv1<-outv1[order(outv1$bat_clon_beam,-outv1$umis.x,-outv1$umis.y),]
outv1<-outv1[!duplicated(outv1$bat_clon_beam),]

outv1$batch_clon<-paste(outv1$batch,outv1$clonotype_id,sep='_')
outv1$freq<-NA
outv1$prop<-NA
outv1$freq<-clonal_freq_all[outv1$batch_clon,'frequency']
outv1$prop<-clonal_freq_all[outv1$batch_clon,'proportion']
stopifnot(sum(is.na(outv1))==0)
outv1<-outv1[outv1$freq>freq_cutoff,]

#plot heatmap for each sample
for (i in sort(unique(outv1$batch))) {
  out.tmp<-outv1[outv1$batch==i,]
  rownames(out.tmp)<-NULL
  out.tmp.summ<-reshape(out.tmp[,c('beam','v1pred','clonotype_id2')],idvar = 'beam',timevar = 'clonotype_id2',direction = 'wide')
  rownames(out.tmp.summ)<-out.tmp.summ$beam
  
  out.tmp.summ<-out.tmp.summ[,-1]
  colnames(out.tmp.summ)<-unlist(lapply(strsplit(colnames(out.tmp.summ),'.',fixed = T), '[[',2))
  pheatmap::pheatmap(t(as.matrix(log10(out.tmp.summ+1e-30))),main=i)
}

#plot all samples
outv1.plot<-outv1[,c('beam','v1pred','clonotype_id2')]
outv1.plot$clonotype_id2<-paste(outv1$batch,outv1$clonotype_id2,sep='_')

out.tmp.summ<-reshape(outv1.plot[,c('beam','v1pred','clonotype_id2')],idvar = 'beam',timevar = 'clonotype_id2',direction = 'wide')
rownames(out.tmp.summ)<-out.tmp.summ$beam
out.tmp.summ<-t(out.tmp.summ[,-1])

# extract truly positive pairs according to specificity score
cutoff<-50
pos_all<-data.frame()
for (i in unique(antigen.spec$batch)) {
  antigen.tmp<-antigen.spec[antigen.spec$batch==i,]
  antigen.tmp.reshape<-reshape(antigen.tmp[,1:10], 
                               direction = "long",
                               varying = list(names(antigen.tmp)[2:10]),
                               v.names = "Value",
                               idvar = "raw_clonotype_id",
                               timevar = "beam",
                               times = paste('BEAM',sprintf("%02d",5:13),sep = ''))
  
  antigen.tmp.reshape<-antigen.tmp.reshape[antigen.tmp.reshape$Value>cutoff,]
  pos_all<-rbind(pos_all,cbind(antigen.tmp.reshape,batch=i))
}
pos_all$bat_clon_beam<-paste(pos_all$batch,pos_all$raw_clonotype_id,pos_all$beam,sep = '_')
rm(list=c('antigen.tmp','antigen.tmp.reshape'))

#extreact specificity score for all clonotype in all samples
spec_score_all<-data.frame()
for (i in unique(antigen.spec$batch)) {
  antigen.tmp<-antigen.spec[antigen.spec$batch==i,]
  antigen.tmp.reshape<-reshape(antigen.tmp[,1:10], 
                               direction = "long",
                               varying = list(names(antigen.tmp)[2:10]),
                               v.names = "Value",
                               idvar = "raw_clonotype_id",
                               timevar = "beam",
                               times = paste('BEAM',sprintf("%02d",5:13),sep = ''))
  spec_score_all<-rbind(spec_score_all,cbind(antigen.tmp.reshape,batch=i))
}

rm(list=c('antigen.tmp','antigen.tmp.reshape'))

spec_score_all$bat_clon_beam<-paste(spec_score_all$batch,spec_score_all$raw_clonotype_id,spec_score_all$beam,sep = '_')
rownames(spec_score_all)<-spec_score_all$bat_clon_beam
outv2.beam09<-outv2[outv2$peptide=='SYIGLKGLYF',]
outv2.beam09$spec_score<-NA
outv2.beam09$spec_score<-spec_score_all[outv2.beam09$bat_clon_beam,'Value']
stopifnot(sum(is.na(outv2.beam09))==0)

summary(outv2.beam09$spec_score[outv2.beam09$avg_externak_rank2>=0 & outv2.beam09$avg_externak_rank2<0.33])
summary(outv2.beam09$spec_score[outv2.beam09$avg_externak_rank2>=0.33 & outv2.beam09$avg_externak_rank2<0.66])
summary(outv2.beam09$spec_score[outv2.beam09$avg_externak_rank2>=0.66])

plot.data<-outv2.beam09[,c('spec_score','avg_externak_rank2')]
plot.data$group<-NA
plot.data$group[plot.data$avg_externak_rank2>=0 & plot.data$avg_externak_rank2<0.33]<-'high pred. affinity'
plot.data$group[plot.data$avg_externak_rank2>=0.33 & plot.data$avg_externak_rank2<0.66]<-'medium pred. affinity'
plot.data$group[plot.data$avg_externak_rank2>=0.66]<-'low pred. affinity'
plot.data$group<-factor(plot.data$group,levels = c('high pred. affinity','medium pred. affinity','low pred. affinity'))

ggplot(plot.data, aes(x=group, y=spec_score,fill=group,color=group)) +
  geom_bar(stat="summary", fun=mean, alpha=0.8,color=NA) +
  stat_summary(fun = mean,
               geom = "errorbar",width=0.3,
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)))+
  labs(x = "",y = "BEAM-T specificity score") +
  theme_classic() + ylim(0,100)+
  theme(axis.text=element_text(size=12), 
        axis.text.x = element_text(angle=15,vjust = 0.6),
        axis.title=element_text(size=14,face="bold"),
        legend.position="none")+
  scale_fill_npg()+scale_color_npg()

# plot auROC AUPR curve
teim<-data.frame(cdr3=input_all$cdr3b[nchar(input_all$cdr3b)<=20],epitope=input_all$peptide[nchar(input_all$cdr3b)<=20])
atmtcr<-data.frame(input_all$peptide,input_all$cdr3b,rep(0))
imrex<-data.frame(cdr3=input_all$cdr3b,antigen.epitope=input_all$peptide)
ergo2<-data.frame(TRA=input_all$cdr3a,TRB=input_all$cdr3b,
                  TRAV=rep(''),TRAJ=rep(''),TRBV=rep(''),
                  TRBJ=rep(''),T_Cell_Type=rep(''),
                  Peptide=input_all$peptide,
                  MHC=rep('HLA-A*24'))

write.table(teim,'input_teim.csv',col.names = T, row.names = F,sep = ',',quote = F)
write.table(atmtcr,'input_atmtcr.csv',col.names = F, row.names = F,sep = ',',quote = F)
write.table(imrex,'input_imrex.csv',col.names = T, row.names = F,sep = ';',quote = F)
write.table(ergo2,'input_ergo2.csv',col.names = T, row.names = F,sep = ',',quote = F)

atmtcr<-read.table('input_atmtcr.csv',header = F,sep=',')
atmtcr<-atmtcr[nchar(atmtcr$V2)<=20,]
write.table(atmtcr,'input_atmtcr_cdr20.csv',col.names = F, row.names = F,sep = ',',quote = F)

outv1<-read.table('input_pmtnet_15folders.csv',header = T,sep=',')
pv1<-read.table('pMTnet_20240803164601__prediction.csv',header = T,sep = ',')
outv1$v1pred<-pv1$Rank
outv1$beam<-NA
outv1$beam[outv1$peptide=='LYVDSLFFL']<-'BEAM05'
outv1$beam[outv1$peptide=='IYNGKLFDL']<-'BEAM06'
outv1$beam[outv1$peptide=='SYQKVIELF']<-'BEAM08'
outv1$beam[outv1$peptide=='SYIGLKGLYF']<-'BEAM09'
outv1$beam[outv1$peptide=='VYVPLKELL']<-'BEAM10'
outv1$beam[outv1$peptide=='VYALPLKML']<-'BEAM11'
outv1$beam[outv1$peptide=='QYDPVAALF']<-'BEAM12'
outv1$beam[outv1$peptide=='AYAQKIFKI']<-'BEAM13'
stopifnot(sum(is.na(outv1))==0)
outv1$bat_clon_beam<-paste(outv1$batch,outv1$clonotype_id,outv1$beam,sep='_')
outv1_subset<-outv1[nchar(outv1$cdr3b)<=20,]

out.teim<-read.table('sequence_level_binding.csv',header = T,sep=',')
stopifnot(sum(out.teim$cdr3!=outv1_subset$cdr3b)==0)
stopifnot(sum(out.teim$epitope!=outv1_subset$peptide)==0)
outv1_subset$teim<-out.teim$binding

out.atm<-read.table('pred_original_input_atmtcr_cdr20.csv',header = F,sep='\t')
out.atm<-out.atm[order(out.atm$V1,out.atm$V2),]
outv1_subset<-outv1_subset[order(outv1_subset$peptide,outv1_subset$cdr3b),]
stopifnot(sum(out.atm$V2!=outv1_subset$cdr3b)==0)
stopifnot(sum(out.atm$V1!=outv1_subset$peptide)==0)
outv1_subset$atm<-out.atm$V5

out.imrex<-read.table('output_imrex.csv',header = T,sep=',')#requires the same inpuit as TEIM that CDR3<=20AA
out.imrex<-out.imrex[order(out.imrex$antigen.epitope,out.imrex$cdr3),]
stopifnot(sum(out.imrex$cdr3!=outv1_subset$cdr3b)==0)
stopifnot(sum(out.imrex$antigen.epitope!=outv1_subset$peptide)==0)
outv1_subset$imrex<-out.imrex$prediction_score

out.ergo2ae.mac<-read.table('ergo2-beamt-results-au-mc.csv',header = T,sep=',')
stopifnot(sum(out.ergo2ae.mac$TRB!=outv1$cdr3b)==0)
stopifnot(sum(out.ergo2ae.mac$TRA!=outv1$cdr3a)==0)
stopifnot(sum(out.ergo2ae.mac$Peptide!=outv1$peptide)==0)
outv1$ergo_ae_mac<-out.ergo2ae.mac$Score

out.ergo2ae.vdj<-read.table('ergo2-beamt-results-au-vdj.csv',header = T,sep=',')
stopifnot(sum(out.ergo2ae.vdj$TRB!=outv1$cdr3b)==0)
stopifnot(sum(out.ergo2ae.vdj$TRA!=outv1$cdr3a)==0)
stopifnot(sum(out.ergo2ae.vdj$Peptide!=outv1$peptide)==0)
outv1$ergo_ae_vdj<-out.ergo2ae.vdj$Score

out.ergo2ls.mac<-read.table('ergo2-beamt-results-lstm-mc.csv',header = T,sep=',')
stopifnot(sum(out.ergo2ls.mac$TRB!=outv1$cdr3b)==0)
stopifnot(sum(out.ergo2ls.mac$TRA!=outv1$cdr3a)==0)
stopifnot(sum(out.ergo2ls.mac$Peptide!=outv1$peptide)==0)
outv1$ergo_ls_mac<-out.ergo2ls.mac$Score

out.ergo2ls.vdj<-read.table('ergo2-beamt-results-lstm-vdj.csv',header = T,sep=',')
stopifnot(sum(out.ergo2ls.vdj$TRA!=outv1$cdr3a)==0)
stopifnot(sum(out.ergo2ls.vdj$TRB!=outv1$cdr3b)==0)
stopifnot(sum(out.ergo2ls.vdj$Peptide!=outv1$peptide)==0)
outv1$ergo_ls_vdj<-out.ergo2ls.vdj$Score

epact.beam<-read.table('preds_epact_origin.csv',header = T,sep=',')
stopifnot(sum(epact.beam$CDR3.alpha.aa!=outv1$cdr3a)==0)
stopifnot(sum(epact.beam$CDR3.beta.aa!=outv1$cdr3b)==0)
stopifnot(sum(epact.beam$Epitope.peptide!=outv1$peptide)==0)
outv1$epact<-epact.beam$Pred

nettcr.beam<-read.table('fig1g_result_nettcr.csv',header = T,sep=',')
nettcr.beam<-nettcr.beam[order(nettcr.beam$cdr3a,nettcr.beam$cdr3b,nettcr.beam$peptide),]
outv1<-outv1[order(outv1$cdr3a,outv1$cdr3b,outv1$peptide),]
stopifnot(sum(nettcr.beam$cdr3a!=outv1$cdr3a)==0)
stopifnot(sum(nettcr.beam$cdr3b!=outv1$cdr3b)==0)
stopifnot(sum(nettcr.beam$peptide!=outv1$peptide)==0)
outv1$nettcr<-nettcr.beam$prediction

lines <- readLines("ranking-tulip-origin.jsonl")
tulip.beam <- lapply(lines, fromJSON)
tulip.beam <- rbind(
  data.frame(peptide=rep(tulip.beam[[1]]$peptide),cdr3a=tulip.beam[[1]]$CDR3a,cdr3b=tulip.beam[[1]]$CDR3b,rank=tulip.beam[[1]]$rank),
  data.frame(peptide=rep(tulip.beam[[2]]$peptide),cdr3a=tulip.beam[[2]]$CDR3a,cdr3b=tulip.beam[[2]]$CDR3b,rank=tulip.beam[[2]]$rank),
  data.frame(peptide=rep(tulip.beam[[3]]$peptide),cdr3a=tulip.beam[[3]]$CDR3a,cdr3b=tulip.beam[[3]]$CDR3b,rank=tulip.beam[[3]]$rank),
  data.frame(peptide=rep(tulip.beam[[4]]$peptide),cdr3a=tulip.beam[[4]]$CDR3a,cdr3b=tulip.beam[[4]]$CDR3b,rank=tulip.beam[[4]]$rank),
  data.frame(peptide=rep(tulip.beam[[5]]$peptide),cdr3a=tulip.beam[[5]]$CDR3a,cdr3b=tulip.beam[[5]]$CDR3b,rank=tulip.beam[[5]]$rank),
  data.frame(peptide=rep(tulip.beam[[6]]$peptide),cdr3a=tulip.beam[[6]]$CDR3a,cdr3b=tulip.beam[[6]]$CDR3b,rank=tulip.beam[[6]]$rank),
  data.frame(peptide=rep(tulip.beam[[7]]$peptide),cdr3a=tulip.beam[[7]]$CDR3a,cdr3b=tulip.beam[[7]]$CDR3b,rank=tulip.beam[[7]]$rank),
  data.frame(peptide=rep(tulip.beam[[8]]$peptide),cdr3a=tulip.beam[[8]]$CDR3a,cdr3b=tulip.beam[[8]]$CDR3b,rank=tulip.beam[[8]]$rank))

tulip.beam<-tulip.beam[order(tulip.beam$cdr3a,tulip.beam$cdr3b,tulip.beam$peptide),]
stopifnot(sum(tulip.beam$cdr3a!=outv1$cdr3a)==0)
stopifnot(sum(tulip.beam$cdr3b!=outv1$cdr3b)==0)
stopifnot(sum(tulip.beam$peptide!=outv1$peptide)==0)
outv1$tulip<-tulip.beam$rank

#remove TCRs with double TRA/TRB
outv1<-outv1[order(outv1$bat_clon_beam,-outv1$umis.x,-outv1$umis.y),]
outv1<-outv1[!duplicated(outv1$bat_clon_beam),]
outv1_subset<-outv1_subset[order(outv1_subset$bat_clon_beam,-outv1_subset$umis.x,-outv1_subset$umis.y),]
outv1_subset<-outv1_subset[!duplicated(outv1_subset$bat_clon_beam),]

#filter according to clone size
outv1$batch_clon<-paste(outv1$batch,outv1$clonotype_id,sep='_')
outv1$freq<-NA
outv1$freq<-clonal_freq_all[outv1$batch_clon,'frequency']
stopifnot(sum(is.na(outv1))==0)
outv1<-outv1[outv1$freq>freq_cutoff,]

outv1_subset$batch_clon<-paste(outv1_subset$batch,outv1_subset$clonotype_id,sep='_')
outv1_subset$freq<-NA
outv1_subset$freq<-clonal_freq_all[outv1_subset$batch_clon,'frequency']
stopifnot(sum(is.na(outv1_subset))==0)
outv1_subset<-outv1_subset[outv1_subset$freq>freq_cutoff,]

#calculate auc
outv1.pos<-outv1[(outv1$bat_clon_beam %in% pos_all$bat_clon_beam),]
outv1.neg<-outv1[!(outv1$bat_clon_beam %in% pos_all$bat_clon_beam),]

print(length(unique(outv1.pos$batch_clon[outv1.pos$peptide=='SYIGLKGLYF'])))
print(length(unique(outv1.neg$batch_clon[outv1.neg$peptide=='SYIGLKGLYF'])))
print(PRROC::roc.curve(1-outv1.pos$v1pred[outv1.pos$peptide=='SYIGLKGLYF'],1-outv1.neg$v1pred[outv1.neg$peptide=='SYIGLKGLYF']))
print(PRROC::roc.curve(outv1.pos$ergo_ae_mac[outv1.pos$peptide=='SYIGLKGLYF'],outv1.neg$ergo_ae_mac[outv1.neg$peptide=='SYIGLKGLYF']))
print(PRROC::roc.curve(outv1.pos$ergo_ae_vdj[outv1.pos$peptide=='SYIGLKGLYF'],outv1.neg$ergo_ae_vdj[outv1.neg$peptide=='SYIGLKGLYF']))
print(PRROC::roc.curve(outv1.pos$ergo_ls_mac[outv1.pos$peptide=='SYIGLKGLYF'],outv1.neg$ergo_ls_mac[outv1.neg$peptide=='SYIGLKGLYF']))
print(PRROC::roc.curve(outv1.pos$ergo_ls_vdj[outv1.pos$peptide=='SYIGLKGLYF'],outv1.neg$ergo_ls_vdj[outv1.neg$peptide=='SYIGLKGLYF']))
print(PRROC::roc.curve(outv1.pos$nettcr[outv1.pos$peptide=='SYIGLKGLYF'],outv1.neg$nettcr[outv1.neg$peptide=='SYIGLKGLYF']))
print(PRROC::roc.curve(outv1.pos$epact[outv1.pos$peptide=='SYIGLKGLYF'],outv1.neg$epact[outv1.neg$peptide=='SYIGLKGLYF']))
print(PRROC::roc.curve(1000-outv1.pos$tulip[outv1.pos$peptide=='SYIGLKGLYF'],1000-outv1.neg$tulip[outv1.neg$peptide=='SYIGLKGLYF']))

outv1_subset.pos<-outv1_subset[(outv1_subset$bat_clon_beam %in% pos_all$bat_clon_beam),]
outv1_subset.neg<-outv1_subset[!(outv1_subset$bat_clon_beam %in% pos_all$bat_clon_beam),]

print(PRROC::roc.curve(outv1_subset.pos$atm[outv1_subset.pos$peptide=='SYIGLKGLYF'],outv1_subset.neg$atm[outv1_subset.neg$peptide=='SYIGLKGLYF']))
print(PRROC::roc.curve(outv1_subset.pos$teim[outv1_subset.pos$peptide=='SYIGLKGLYF'],outv1_subset.neg$teim[outv1_subset.neg$peptide=='SYIGLKGLYF']))
print(PRROC::roc.curve(outv1_subset.pos$imrex[outv1_subset.pos$peptide=='SYIGLKGLYF'],outv1_subset.neg$imrex[outv1_subset.neg$peptide=='SYIGLKGLYF']))

#auprc
freq_cutoff<-5
cutoff<-15

antigen.spec<-read.table('summary_antigen_specificity.csv',header = T,sep = ',')
clonal_freq_all<-data.frame()
for (pdir in list.dirs('BEAM-T-Cassian',recursive = F)) {
  folder_name<-paste(paste(unlist(strsplit(pdir,'[/_]'))[c(11,12,13)],collapse = '_'),tail(unlist(strsplit(pdir,'[/_]')),n=1),sep = '_')
  # TCR
  k5.t<-read.table(
    list.files(path = pdir,pattern = 'clonotypes.csv',full.names = T,recursive = T),
    header = T,sep=',')
  clonal_freq_all<-rbind(clonal_freq_all,cbind(k5.t[,c("clonotype_id","frequency","proportion")],batch=rep(folder_name)))
}
clonal_freq_all$batch_clon<-paste(clonal_freq_all$batch,clonal_freq_all$clonotype_id,sep='_')
rownames(clonal_freq_all)<-clonal_freq_all$batch_clon
rm(list=c('k5.t','folder_name','pdir'))

outv2<-read.table('results_pmtnet_15folders_v2.csv',header = T,sep = ',')
outv2$beam<-NA
outv2$beam[outv2$peptide=='LYVDSLFFL']<-'BEAM05'
outv2$beam[outv2$peptide=='IYNGKLFDL']<-'BEAM06'
outv2$beam[outv2$peptide=='SYQKVIELF']<-'BEAM08'
outv2$beam[outv2$peptide=='SYIGLKGLYF']<-'BEAM09'
outv2$beam[outv2$peptide=='VYVPLKELL']<-'BEAM10'
outv2$beam[outv2$peptide=='VYALPLKML']<-'BEAM11'
outv2$beam[outv2$peptide=='QYDPVAALF']<-'BEAM12'
outv2$beam[outv2$peptide=='AYAQKIFKI']<-'BEAM13'
stopifnot(sum(is.na(outv2))==0)
outv2$bat_clon_beam<-paste(outv2$batch,outv2$clonotype_id,outv2$beam,sep='_')
outv2<-outv2[order(outv2$bat_clon_beam,-outv2$umis.x,-outv2$umis.y),]
outv2<-outv2[!duplicated(outv2$bat_clon_beam),]
outv2$tcr<-paste(outv2$va,outv2$cdr3a,outv2$vb,outv2$cdr3b,sep='_')
outv2$batch_clon<-paste(outv2$batch,outv2$clonotype_id,sep='_')
outv2$freq<-NA
outv2$freq<-clonal_freq_all[outv2$batch_clon,'frequency']
stopifnot(sum(is.na(outv2))==0)
outv2<-outv2[outv2$freq>freq_cutoff,]

pos_all<-data.frame()
for (i in unique(antigen.spec$batch)) {
  antigen.tmp<-antigen.spec[antigen.spec$batch==i,]
  antigen.tmp.reshape<-reshape(antigen.tmp[,1:10], 
                               direction = "long",
                               varying = list(names(antigen.tmp)[2:10]),
                               v.names = "Value",
                               idvar = "raw_clonotype_id",
                               timevar = "beam",
                               times = paste('BEAM',sprintf("%02d",5:13),sep = ''))
  antigen.tmp.reshape<-antigen.tmp.reshape[antigen.tmp.reshape$Value>cutoff,]
  pos_all<-rbind(pos_all,cbind(antigen.tmp.reshape,batch=i))
}
pos_all$bat_clon_beam<-paste(pos_all$batch,pos_all$raw_clonotype_id,pos_all$beam,sep = '_')
rm(list=c('antigen.tmp','antigen.tmp.reshape'))

outv2.pos<-outv2[(outv2$bat_clon_beam %in% pos_all$bat_clon_beam),]
outv2.neg<-outv2[!(outv2$bat_clon_beam %in% pos_all$bat_clon_beam),]

outv1<-read.table('input_pmtnet_15folders.csv',header = T,sep=',')
pv1<-read.table('pMTnet_20240803164601__prediction.csv',header = T,sep = ',')
outv1$v1pred<-pv1$Rank
outv1$beam<-NA
outv1$beam[outv1$peptide=='LYVDSLFFL']<-'BEAM05'
outv1$beam[outv1$peptide=='IYNGKLFDL']<-'BEAM06'
outv1$beam[outv1$peptide=='SYQKVIELF']<-'BEAM08'
outv1$beam[outv1$peptide=='SYIGLKGLYF']<-'BEAM09'
outv1$beam[outv1$peptide=='VYVPLKELL']<-'BEAM10'
outv1$beam[outv1$peptide=='VYALPLKML']<-'BEAM11'
outv1$beam[outv1$peptide=='QYDPVAALF']<-'BEAM12'
outv1$beam[outv1$peptide=='AYAQKIFKI']<-'BEAM13'
stopifnot(sum(is.na(outv1))==0)
outv1$bat_clon_beam<-paste(outv1$batch,outv1$clonotype_id,outv1$beam,sep='_')
outv1_subset<-outv1[nchar(outv1$cdr3b)<=20,]

out.teim<-read.table('sequence_level_binding.csv',header = T,sep=',')
stopifnot(sum(out.teim$cdr3!=outv1_subset$cdr3b)==0)
stopifnot(sum(out.teim$epitope!=outv1_subset$peptide)==0)
outv1_subset$teim<-out.teim$binding

out.atm<-read.table('pred_original_input_atmtcr_cdr20.csv',header = F,sep='\t')
out.atm<-out.atm[order(out.atm$V1,out.atm$V2),]
outv1_subset<-outv1_subset[order(outv1_subset$peptide,outv1_subset$cdr3b),]
stopifnot(sum(out.atm$V2!=outv1_subset$cdr3b)==0)
stopifnot(sum(out.atm$V1!=outv1_subset$peptide)==0)
outv1_subset$atm<-out.atm$V5

out.imrex<-read.table('output_imrex.csv',header = T,sep=',')#requires the same inpuit as TEIM that CDR3<=20AA
out.imrex<-out.imrex[order(out.imrex$antigen.epitope,out.imrex$cdr3),]
stopifnot(sum(out.imrex$cdr3!=outv1_subset$cdr3b)==0)
stopifnot(sum(out.imrex$antigen.epitope!=outv1_subset$peptide)==0)
outv1_subset$imrex<-out.imrex$prediction_score

out.ergo2ae.mac<-read.table('ergo2-beamt-results-au-mc.csv',header = T,sep=',')
stopifnot(sum(out.ergo2ae.mac$TRB!=outv1$cdr3b)==0)
stopifnot(sum(out.ergo2ae.mac$TRA!=outv1$cdr3a)==0)
stopifnot(sum(out.ergo2ae.mac$Peptide!=outv1$peptide)==0)
outv1$ergo_ae_mac<-out.ergo2ae.mac$Score

out.ergo2ae.vdj<-read.table('ergo2-beamt-results-au-vdj.csv',header = T,sep=',')
stopifnot(sum(out.ergo2ae.vdj$TRB!=outv1$cdr3b)==0)
stopifnot(sum(out.ergo2ae.vdj$TRA!=outv1$cdr3a)==0)
stopifnot(sum(out.ergo2ae.vdj$Peptide!=outv1$peptide)==0)
outv1$ergo_ae_vdj<-out.ergo2ae.vdj$Score

out.ergo2ls.mac<-read.table('ergo2-beamt-results-lstm-mc.csv',header = T,sep=',')
stopifnot(sum(out.ergo2ls.mac$TRB!=outv1$cdr3b)==0)
stopifnot(sum(out.ergo2ls.mac$TRA!=outv1$cdr3a)==0)
stopifnot(sum(out.ergo2ls.mac$Peptide!=outv1$peptide)==0)
outv1$ergo_ls_mac<-out.ergo2ls.mac$Score

out.ergo2ls.vdj<-read.table('ergo2-beamt-results-lstm-vdj.csv',header = T,sep=',')
stopifnot(sum(out.ergo2ls.vdj$TRA!=outv1$cdr3a)==0)
stopifnot(sum(out.ergo2ls.vdj$TRB!=outv1$cdr3b)==0)
stopifnot(sum(out.ergo2ls.vdj$Peptide!=outv1$peptide)==0)
outv1$ergo_ls_vdj<-out.ergo2ls.vdj$Score

epact.beam<-read.table('preds_epact_origin.csv',header = T,sep=',')
stopifnot(sum(epact.beam$CDR3.alpha.aa!=outv1$cdr3a)==0)
stopifnot(sum(epact.beam$CDR3.beta.aa!=outv1$cdr3b)==0)
stopifnot(sum(epact.beam$Epitope.peptide!=outv1$peptide)==0)
outv1$epact<-epact.beam$Pred

nettcr.beam<-read.table('fig1g_result_nettcr.csv',header = T,sep=',')
nettcr.beam<-nettcr.beam[order(nettcr.beam$cdr3a,nettcr.beam$cdr3b,nettcr.beam$peptide),]
outv1<-outv1[order(outv1$cdr3a,outv1$cdr3b,outv1$peptide),]
stopifnot(sum(nettcr.beam$cdr3a!=outv1$cdr3a)==0)
stopifnot(sum(nettcr.beam$cdr3b!=outv1$cdr3b)==0)
stopifnot(sum(nettcr.beam$peptide!=outv1$peptide)==0)
outv1$nettcr<-nettcr.beam$prediction

lines <- readLines("ranking-tulip-origin.jsonl")
tulip.beam <- lapply(lines, fromJSON)
tulip.beam <- rbind(
  data.frame(peptide=rep(tulip.beam[[1]]$peptide),cdr3a=tulip.beam[[1]]$CDR3a,cdr3b=tulip.beam[[1]]$CDR3b,rank=tulip.beam[[1]]$rank),
  data.frame(peptide=rep(tulip.beam[[2]]$peptide),cdr3a=tulip.beam[[2]]$CDR3a,cdr3b=tulip.beam[[2]]$CDR3b,rank=tulip.beam[[2]]$rank),
  data.frame(peptide=rep(tulip.beam[[3]]$peptide),cdr3a=tulip.beam[[3]]$CDR3a,cdr3b=tulip.beam[[3]]$CDR3b,rank=tulip.beam[[3]]$rank),
  data.frame(peptide=rep(tulip.beam[[4]]$peptide),cdr3a=tulip.beam[[4]]$CDR3a,cdr3b=tulip.beam[[4]]$CDR3b,rank=tulip.beam[[4]]$rank),
  data.frame(peptide=rep(tulip.beam[[5]]$peptide),cdr3a=tulip.beam[[5]]$CDR3a,cdr3b=tulip.beam[[5]]$CDR3b,rank=tulip.beam[[5]]$rank),
  data.frame(peptide=rep(tulip.beam[[6]]$peptide),cdr3a=tulip.beam[[6]]$CDR3a,cdr3b=tulip.beam[[6]]$CDR3b,rank=tulip.beam[[6]]$rank),
  data.frame(peptide=rep(tulip.beam[[7]]$peptide),cdr3a=tulip.beam[[7]]$CDR3a,cdr3b=tulip.beam[[7]]$CDR3b,rank=tulip.beam[[7]]$rank),
  data.frame(peptide=rep(tulip.beam[[8]]$peptide),cdr3a=tulip.beam[[8]]$CDR3a,cdr3b=tulip.beam[[8]]$CDR3b,rank=tulip.beam[[8]]$rank))

tulip.beam<-tulip.beam[order(tulip.beam$cdr3a,tulip.beam$cdr3b,tulip.beam$peptide),]
stopifnot(sum(tulip.beam$cdr3a!=outv1$cdr3a)==0)
stopifnot(sum(tulip.beam$cdr3b!=outv1$cdr3b)==0)
stopifnot(sum(tulip.beam$peptide!=outv1$peptide)==0)
outv1$tulip<-tulip.beam$rank

#remove TCRs with double TRA/TRB
outv1<-outv1[order(outv1$bat_clon_beam,-outv1$umis.x,-outv1$umis.y),]
outv1<-outv1[!duplicated(outv1$bat_clon_beam),]

outv1_subset<-outv1_subset[order(outv1_subset$bat_clon_beam,-outv1_subset$umis.x,-outv1_subset$umis.y),]
outv1_subset<-outv1_subset[!duplicated(outv1_subset$bat_clon_beam),]

#filter according to clone size
outv1$batch_clon<-paste(outv1$batch,outv1$clonotype_id,sep='_')
outv1$freq<-NA
outv1$freq<-clonal_freq_all[outv1$batch_clon,'frequency']
stopifnot(sum(is.na(outv1))==0)
outv1<-outv1[outv1$freq>freq_cutoff,]

outv1_subset$batch_clon<-paste(outv1_subset$batch,outv1_subset$clonotype_id,sep='_')
outv1_subset$freq<-NA
outv1_subset$freq<-clonal_freq_all[outv1_subset$batch_clon,'frequency']
stopifnot(sum(is.na(outv1_subset))==0)
outv1_subset<-outv1_subset[outv1_subset$freq>freq_cutoff,]

#calculate auc
outv1.pos<-outv1[(outv1$bat_clon_beam %in% pos_all$bat_clon_beam),]
outv1.neg<-outv1[!(outv1$bat_clon_beam %in% pos_all$bat_clon_beam),]

outv1_subset.pos<-outv1_subset[(outv1_subset$bat_clon_beam %in% pos_all$bat_clon_beam),]
outv1_subset.neg<-outv1_subset[!(outv1_subset$bat_clon_beam %in% pos_all$bat_clon_beam),]

pheatmap::pheatmap(as.matrix(plot.data),cluster_rows = F,cluster_cols = F,)



#Fig. 2
#pMtnet-Omni
pmtnet_all=read.csv("ensembl_results_10x_drop_fig2_data_stage13_epoch2.csv",stringsAsFactors = F) # complete and unduplicated 10X pairs
pmtnet_all=pmtnet_all[pmtnet_all$group=="training",]

human_a<-Biostrings::readAAStringSet(filepath = "human_trav_pro.txt.completed.txt")
human_name_a <- names(human_a)
human_seq_a <- paste(human_a)
va<-data.frame(a_allele=sapply(strsplit(human_name_a,split = "|",fixed = T),function(x){x[2]}), human_seq_a)
va<-rbind(va,data.frame(a_allele="XXX",human_seq_a="XXX"))

human_b<-Biostrings::readAAStringSet(filepath = "human_trbv_pro.txt.completed.txt")
human_name_b <- names(human_b)
human_seq_b <- paste(human_b)
vb<-data.frame(b_allele=sapply(strsplit(human_name_b,split = "|",fixed = T),function(x){x[2]}), human_seq_b)
vb<-rbind(vb,data.frame(b_allele="XXX",human_seq_b="XXX"))


for (i in 1:dim(pmtnet_all)[1])
{
  v=unique(va[va$human_seq_a==pmtnet_all$vaseq[i],"a_allele"])
  stopifnot(length(unique(unlist(lapply(strsplit(v,split = "*",fixed = T),`[[`,1))))==1)
  pmtnet_all$va[i]=unique(unlist(lapply(strsplit(v,split = "*",fixed = T),`[[`,1)))
}

for (i in 1:dim(pmtnet_all)[1])
{
  v=unique(vb[vb$human_seq_b==pmtnet_all$vbseq[i],"b_allele"])
  if (length(unique(unlist(lapply(strsplit(v,split = "*",fixed = T),`[[`,1))))>1) {
    if(all(v==c("TRBV6-2*01", "TRBV6-3*01"))){
      pmtnet_all$vb[i]="TRBV6-3"
    }else{
      print(v)
      stop("Error")
    }
  } else if (length(unique(unlist(lapply(strsplit(v,split = "*",fixed = T),`[[`,1))))==1){
    pmtnet_all$vb[i]=unique(unlist(lapply(strsplit(v,split = "*",fixed = T),`[[`,1)))
  } else{
    print(v)
    stop("Error")
  }
  
}

pmtnet_all$va[pmtnet_all$va=="TRAV14/DV4"]<-"TRAV14DV4"
pmtnet_all$va[pmtnet_all$va=="TRAV38-2/DV8"]<-"TRAV38-2DV8"
pmtnet_all$va[pmtnet_all$va=="TRAV29/DV5"]<-"TRAV29DV5"
pmtnet_all$va[pmtnet_all$va=="TRAV36/DV7"]<-"TRAV36DV7"
pmtnet_all$va[pmtnet_all$va=="TRAV23/DV6"]<-"TRAV23DV6"

contig<-read.table("vdj_v1_hs_aggregated_donor1_all_contig_annotations.csv",header = T,sep=",")
stopifnot(sum(!(pmtnet_all$va %in% contig$v_gene))==0)
stopifnot(sum(!(pmtnet_all$vb %in% contig$v_gene))==0)

pmtnet_all=pmtnet_all[(!is.na(pmtnet_all$va)) & (!is.na(pmtnet_all$vb)),]
pmtnet_all$clone=apply(pmtnet_all[,c("va","vb","cdr3a","cdr3b")],
                       1,function(x) paste(x,collapse="_"))
pmtnet_all$pMHC=sapply(strsplit(pmtnet_all$X_pmhc,"_"),
                       function(x) paste(sub("\\:","",sub("\\*","",x[2],perl=T),perl=T),x[1],sep="_"))
rownames(pmtnet_all)=
  apply(pmtnet_all[,c("clone","pMHC")],1,function(x) paste(x,collapse="_"))
rm(list=setdiff(ls(),c("pmtnet_all")))

all_comparison_all_info<-read.table("all_comparison_mean2.txt",header = T,sep = "\t",stringsAsFactors = F)
all_comparison_all_info$better_pred<-NA
all_comparison_all_info$worse_pred<-NA
all_comparison_all_info$better_pred=pmtnet_all[apply(all_comparison_all_info[,c("better_clone","pMHC")],1,function(x) paste(x,collapse="_")),"avg_rank"]
all_comparison_all_info$worse_pred=pmtnet_all[apply(all_comparison_all_info[,c("worse_clone","pMHC")],1,function(x) paste(x,collapse="_")),"avg_rank"]
stopifnot(sum(is.na(all_comparison_all_info[,c("worseCDR3a","worseCDR3b","worse_umi","betterCDR3a","betterCDR3b","better_umi","Va","Vb","pMHC","better_clone","worse_clone")]))==0)

get_ci_est_small<-function(fold,data,type){
  if(sum((data$better_umi/data$worse_umi)>fold)>1){
    tmp<-confintr::ci_proportion(x=as.numeric(data[(data$better_umi/data$worse_umi)>fold,"better_pred"]<data[(data$better_umi/data$worse_umi)>fold,"worse_pred"]),type=CI_TYPE)
    if(type=="lower"){
      return(tmp$interval[1])
    } else if(type=="upper"){
      return(tmp$interval[2])
    }else{
      stop("Type is wrong!")
    }
  }else{
    return(NA)
  }
}

umi_cutoff<-30
min_umi_better<-umi_cutoff
min_umi_worse<-umi_cutoff

data_tmp1<-all_comparison_all_info[(all_comparison_all_info$dista + all_comparison_all_info$distb)>=1 & 
                                     (all_comparison_all_info$dista + all_comparison_all_info$distb)<=3 & 
                                     all_comparison_all_info$better_umi>min_umi_better & all_comparison_all_info$worse_umi>min_umi_worse,]
data_tmp2<-all_comparison_all_info[(all_comparison_all_info$dista + all_comparison_all_info$distb)>=4 & 
                                     (all_comparison_all_info$dista + all_comparison_all_info$distb)<=5 & 
                                     all_comparison_all_info$better_umi>min_umi_better & all_comparison_all_info$worse_umi>min_umi_worse,]
data_tmp3<-all_comparison_all_info[(all_comparison_all_info$dista + all_comparison_all_info$distb)>=1 & 
                                     (all_comparison_all_info$dista + all_comparison_all_info$distb)<=5 & 
                                     all_comparison_all_info$better_umi>min_umi_better & all_comparison_all_info$worse_umi>min_umi_worse,]

acc12<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(mean(data_tmp1[(data_tmp1$better_umi/data_tmp1$worse_umi)>x,"better_pred"]<data_tmp1[(data_tmp1$better_umi/data_tmp1$worse_umi)>x,"worse_pred"]))}))
acc35<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(mean(data_tmp2[(data_tmp2$better_umi/data_tmp2$worse_umi)>x,"better_pred"]<data_tmp2[(data_tmp2$better_umi/data_tmp2$worse_umi)>x,"worse_pred"]))}))
acc15<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(mean(data_tmp3[(data_tmp3$better_umi/data_tmp3$worse_umi)>x,"better_pred"]<data_tmp3[(data_tmp3$better_umi/data_tmp3$worse_umi)>x,"worse_pred"]))}))
lower12<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_small,data=data_tmp1,type="lower"))
lower35<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_small,data=data_tmp2,type="lower"))
lower15<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_small,data=data_tmp3,type="lower"))
upper12<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_small,data=data_tmp1,type="upper"))
upper35<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_small,data=data_tmp2,type="upper"))
upper15<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_small,data=data_tmp3,type="upper"))
num12<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(length(data_tmp1[(data_tmp1$better_umi/data_tmp1$worse_umi)>x,"better_pred"]))}))
num35<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(length(data_tmp2[(data_tmp2$better_umi/data_tmp2$worse_umi)>x,"better_pred"]))}))
num15<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(length(data_tmp3[(data_tmp3$better_umi/data_tmp3$worse_umi)>x,"better_pred"]))}))

plot_data2<-data.frame(fold=rep(c(seq(2,6.5,0.125)),3), acc=c(acc12,acc35,acc15),lower=c(lower12,lower35,lower15),upper=c(upper12,upper35,upper15),num=c(num12,num35,num15),group=c(rep("1~3_AA",37),rep("4~5_AA",37),rep("1~5_AA",37)))

# TEIM
validation<-read.table("validation.txt",header = T,sep = "\t")
validation<-validation[order(validation$cdr3a,validation$cdr3b,validation$peptide,validation$vaseq,validation$vbseq,validation$mhca,validation$mhcb),]
data_cdr3b_epitope<-validation[validation$pMHC_SPECIES=="human" & validation$TCR_SPECIES=="human",c("cdr3b","peptide","label")]
colnames(data_cdr3b_epitope)<-c("cdr3","epitope","label")
data_cdr3b_epitope<-data_cdr3b_epitope[(nchar(data_cdr3b_epitope$epitope)<=12)&(nchar(data_cdr3b_epitope$cdr3)<=20),]
data_cdr3b_epitope<-data_cdr3b_epitope[!duplicated(data_cdr3b_epitope),]
write.table(data_cdr3b_epitope[,1:2],"input_for_teim_our_val.txt",quote = F,sep = ",",row.names = F)

teim_out<-read.table("sequence_level_binding.csv",header = T,sep = ",")
stopifnot(sum(teim_out$cdr3!=data_cdr3b_epitope$cdr3)==0)
stopifnot(sum(teim_out$epitope!=data_cdr3b_epitope$epitope)==0)

pROC::auc(data_cdr3b_epitope$label,teim_out$binding)
PRROC::roc.curve(teim_out$binding[data_cdr3b_epitope$label==1],teim_out$binding[data_cdr3b_epitope$label==0])
PRROC::pr.curve(teim_out$binding[data_cdr3b_epitope$label==1],teim_out$binding[data_cdr3b_epitope$label==0])
pred_obj <- ROCR::prediction(teim_out$binding, data_cdr3b_epitope$label)
perf_obj <- ROCR::performance(pred_obj, measure = "prec", x.measure = "rec")
ROCR::plot(perf_obj, ylim = c(0,1))
qwe<-ROCR::performance(pred_obj, measure = "aucpr")
qwe@y.values

#10X data
input<-read.table("input_data_continuous_better_worse_for_benchmark_software_have_vname.txt",header = T,sep="\t")
input_hc1<-input[input$pMHC_SPECIES=="human" & input$TCR_SPECIES=="human" & input$MHC_BETA=="human_glublin",c('TRA_V','CDR3_ALPHA','TRB_V','CDR3_BETA','EPITOPE','MHC_ALPHA')]
colnames(input_hc1)<-c("TRA_V","CDR3_ALPHA","TRB_V","cdr3","epitope","MHC_ALPHA")
input_hc1<-input_hc1[(nchar(input_hc1$epitope)<=12)&(nchar(input_hc1$cdr3)<=20),]
input_hc1<-input_hc1[!duplicated(input_hc1),]

all_comparison_all_info<-read.table("all_comparison_mean2.txt",header = T,sep = "\t",stringsAsFactors = F)
all_comparison_all_info$better_pred<-NA
all_comparison_all_info$worse_pred<-NA

teim<-read.table("sequence_level_binding.csv",header = T,sep = ",")
stopifnot(sum(input_hc1$cdr3!=teim$cdr3)==0)
stopifnot(sum(input_hc1$epitope!=teim$epitope)==0)
teim<-cbind(Score=teim$binding,input_hc1)

teim$TRAV<-unlist(lapply(strsplit(teim$TRA_V,"*",fixed = T),`[`,1))
teim$TRBV<-unlist(lapply(strsplit(teim$TRB_V,"*",fixed = T),`[`,1))
teim$pmhc<-paste(unlist(lapply(lapply(strsplit(teim$MHC_ALPHA,"[-:\\*]"),`[`,c(2,3,4)), paste,collapse="")),teim$epitope,sep = "_")
teim$tcr<-paste(teim$TRAV,teim$TRBV,teim$CDR3_ALPHA,teim$cdr3,sep="_")
teim$pmhctcr<-paste(teim$pmhc,teim$tcr,sep="_")
rownames(teim)<-teim$pmhctcr

all_comparison_all_info$better_pred<-teim[paste(all_comparison_all_info$pMHC,all_comparison_all_info$Va,all_comparison_all_info$Vb,all_comparison_all_info$betterCDR3a,all_comparison_all_info$betterCDR3b,sep = "_"),"Score"]
all_comparison_all_info$worse_pred<-teim[paste(all_comparison_all_info$pMHC,all_comparison_all_info$Va,all_comparison_all_info$Vb,all_comparison_all_info$worseCDR3a,all_comparison_all_info$worseCDR3b,sep = "_"),"Score"]
all_comparison_all_info<-na.omit(all_comparison_all_info)

data_tmp1<-all_comparison_all_info[(all_comparison_all_info$dista + all_comparison_all_info$distb)>=1 & 
                                     (all_comparison_all_info$dista + all_comparison_all_info$distb)<=3 & 
                                     all_comparison_all_info$better_umi>30 & all_comparison_all_info$worse_umi>30,]
data_tmp2<-all_comparison_all_info[(all_comparison_all_info$dista + all_comparison_all_info$distb)>=4 & 
                                     (all_comparison_all_info$dista + all_comparison_all_info$distb)<=5 & 
                                     all_comparison_all_info$better_umi>30 & all_comparison_all_info$worse_umi>30,]
data_tmp3<-all_comparison_all_info[(all_comparison_all_info$dista + all_comparison_all_info$distb)>=1 & 
                                     (all_comparison_all_info$dista + all_comparison_all_info$distb)<=5 & 
                                     all_comparison_all_info$better_umi>30 & all_comparison_all_info$worse_umi>30,]

get_ci_est_large<-function(fold,data,type){
  if(sum((data$better_umi/data$worse_umi)>fold)>1){
    tmp<-confintr::ci_proportion(x=as.numeric(data[(data$better_umi/data$worse_umi)>fold,"better_pred"]>data[(data$better_umi/data$worse_umi)>fold,"worse_pred"]),type="Wilson")
    if(type=="lower"){
      return(tmp$interval[1])
    } else if(type=="upper"){
      return(tmp$interval[2])
    }else{
      stop("Type is wrong!")
    }
  }else{
    return(NA)
  }
}
acc12<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(mean(data_tmp1[(data_tmp1$better_umi/data_tmp1$worse_umi)>x,"better_pred"]>data_tmp1[(data_tmp1$better_umi/data_tmp1$worse_umi)>x,"worse_pred"]))}))
acc35<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(mean(data_tmp2[(data_tmp2$better_umi/data_tmp2$worse_umi)>x,"better_pred"]>data_tmp2[(data_tmp2$better_umi/data_tmp2$worse_umi)>x,"worse_pred"]))}))
acc15<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(mean(data_tmp3[(data_tmp3$better_umi/data_tmp3$worse_umi)>x,"better_pred"]>data_tmp3[(data_tmp3$better_umi/data_tmp3$worse_umi)>x,"worse_pred"]))}))
lower12<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp1,type="lower"))
lower35<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp2,type="lower"))
lower15<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp3,type="lower"))
upper12<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp1,type="upper"))
upper35<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp2,type="upper"))
upper15<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp3,type="upper"))
num12<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(length(data_tmp1[(data_tmp1$better_umi/data_tmp1$worse_umi)>x,"better_pred"]))}))
num35<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(length(data_tmp2[(data_tmp2$better_umi/data_tmp2$worse_umi)>x,"better_pred"]))}))
num15<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(length(data_tmp3[(data_tmp3$better_umi/data_tmp3$worse_umi)>x,"better_pred"]))}))

plot_data2_teim1<-data.frame(fold=rep(c(seq(2,6.5,0.125)),3), acc=c(acc12,acc35,acc15),lower=c(lower12,lower35,lower15),upper=c(upper12,upper35,upper15),num=c(num12,num35,num15),group=c(rep("1~3_AA",37),rep("4~5_AA",37),rep("1~5_AA",37)))

#pmtnet v1
pmtnetv1_data<-read.table("pmtnet_v1_data.csv",header = T,sep=",")
pmtnetv1_data$cdr3pephla<-paste(pmtnetv1_data$CDR3,pmtnetv1_data$Antigen,pmtnetv1_data$HLA,sep="_")
pmtnetv1_data<-pmtnetv1_data[!duplicated(pmtnetv1_data$cdr3pephla),]
pmtnetv1_data<-pmtnetv1_data[pmtnetv1_data$cdr3pephla %in% result1_hc1$cdr3bpephla,]

pmtnetv1<-read.table("pMTnet_prediction.csv",header = T,sep=",")
pmtnetv1$cdr3pephla<-paste(pmtnetv1$CDR3,pmtnetv1$Antigen,pmtnetv1$HLA,sep="_")
pmtnetv1<-pmtnetv1[!duplicated(pmtnetv1$cdr3pephla),]
pmtnetv1<-pmtnetv1[pmtnetv1$cdr3pephla %in% result1_hc1$cdr3bpephla,]

stopifnot(sum(pmtnetv1$CDR3!=pmtnetv1_data$CDR3)==0)
stopifnot(sum(pmtnetv1$Antigen!=pmtnetv1_data$Antigen)==0)
stopifnot(sum(pmtnetv1$HLA!=pmtnetv1_data$HLA)==0)

pROC::auc(pmtnetv1_data$label,pmtnetv1$Rank)
PRROC::roc.curve(1-pmtnetv1$Rank[pmtnetv1_data$label==1],1-pmtnetv1$Rank[pmtnetv1_data$label==0])
PRROC::pr.curve(1-pmtnetv1$Rank[pmtnetv1_data$label==1],1-pmtnetv1$Rank[pmtnetv1_data$label==0])
plot(PRROC::pr.curve(1-pmtnetv1$Rank[pmtnetv1_data$label==1],1-pmtnetv1$Rank[pmtnetv1_data$label==0],curve = T))
pred_obj <- ROCR::prediction(1-pmtnetv1$Rank, pmtnetv1_data$label)
perf_obj <- ROCR::performance(pred_obj, measure = "prec", x.measure = "rec")
ROCR::plot(perf_obj, ylim = c(0,1))
qwe<-ROCR::performance(pred_obj, measure = "aucpr")
qwe@y.values

#better worse
get_ci_est_small<-function(fold,data,type){
  if(sum((data$better_umi/data$worse_umi)>fold)>1){
    tmp<-confintr::ci_proportion(x=as.numeric(data[(data$better_umi/data$worse_umi)>fold,"better_pred"]<data[(data$better_umi/data$worse_umi)>fold,"worse_pred"]),type="Wilson")
    if(type=="lower"){
      return(tmp$interval[1])
    } else if(type=="upper"){
      return(tmp$interval[2])
    }else{
      stop("Type is wrong!")
    }
  }else{
    return(NA)
  }
}

input<-read.table("input_data_continuous_better_worse_for_benchmark_software_have_vname.txt",header = T,sep="\t")
input_hc1<-input[input$pMHC_SPECIES=="human" & input$TCR_SPECIES=="human" & input$MHC_BETA=="human_glublin",c('TRA_V','CDR3_ALPHA','TRB_V','CDR3_BETA','EPITOPE','MHC_ALPHA')]
colnames(input_hc1)<-c("TRA_V","CDR3_ALPHA","TRB_V","CDR3","Antigen","HLA")

input_hc1<-input_hc1[!duplicated(input_hc1),]
input_hc1$HLA<-unlist(lapply(strsplit(input_hc1$HLA,"-"),`[`,2))

all_comparison_all_info<-read.table("all_comparison_mean2.txt",header = T,sep = "\t",stringsAsFactors = F)
all_comparison_all_info$better_pred<-NA
all_comparison_all_info$worse_pred<-NA

pmtnetv1<-read.table("prediction_pmtnetv1_better_worse.csv",header = T,sep = ",")
stopifnot(sum(pmtnetv1$CDR3!=input_hc1$CDR3)==0)
stopifnot(sum(pmtnetv1$Antigen!=input_hc1$Antigen)==0)
stopifnot(sum(pmtnetv1$HLA!=input_hc1$HLA)==0)

pmtnetv1<-cbind(input_hc1,rank=pmtnetv1$Rank)
pmtnetv1$TRAV<-unlist(lapply(strsplit(pmtnetv1$TRA_V,"*",fixed = T),`[`,1))
pmtnetv1$TRBV<-unlist(lapply(strsplit(pmtnetv1$TRB_V,"*",fixed = T),`[`,1))
pmtnetv1$pmhc<-paste(unlist(lapply(lapply(strsplit(pmtnetv1$HLA,"[:\\*]"),`[`,c(1,2,3)), paste,collapse="")),pmtnetv1$Antigen,sep = "_")
pmtnetv1$tcr<-paste(pmtnetv1$TRAV,pmtnetv1$TRBV,pmtnetv1$CDR3_ALPHA,pmtnetv1$CDR3,sep="_")
pmtnetv1$pmhctcr<-paste(pmtnetv1$pmhc,pmtnetv1$tcr,sep="_")
rownames(pmtnetv1)<-pmtnetv1$pmhctcr

all_comparison_all_info$better_pred<-pmtnetv1[paste(all_comparison_all_info$pMHC,all_comparison_all_info$Va,all_comparison_all_info$Vb,all_comparison_all_info$betterCDR3a,all_comparison_all_info$betterCDR3b,sep = "_"),"rank"]
all_comparison_all_info$worse_pred<-pmtnetv1[paste(all_comparison_all_info$pMHC,all_comparison_all_info$Va,all_comparison_all_info$Vb,all_comparison_all_info$worseCDR3a,all_comparison_all_info$worseCDR3b,sep = "_"),"rank"]
all_comparison_all_info<-na.omit(all_comparison_all_info)

data_tmp1<-all_comparison_all_info[(all_comparison_all_info$dista + all_comparison_all_info$distb)>=1 & 
                                     (all_comparison_all_info$dista + all_comparison_all_info$distb)<=3 & 
                                     all_comparison_all_info$better_umi>30 & all_comparison_all_info$worse_umi>30,]
data_tmp2<-all_comparison_all_info[(all_comparison_all_info$dista + all_comparison_all_info$distb)>=4 & 
                                     (all_comparison_all_info$dista + all_comparison_all_info$distb)<=5 & 
                                     all_comparison_all_info$better_umi>30 & all_comparison_all_info$worse_umi>30,]
data_tmp3<-all_comparison_all_info[(all_comparison_all_info$dista + all_comparison_all_info$distb)>=1 & 
                                     (all_comparison_all_info$dista + all_comparison_all_info$distb)<=5 & 
                                     all_comparison_all_info$better_umi>30 & all_comparison_all_info$worse_umi>30,]

acc12<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(mean(data_tmp1[(data_tmp1$better_umi/data_tmp1$worse_umi)>x,"better_pred"]<data_tmp1[(data_tmp1$better_umi/data_tmp1$worse_umi)>x,"worse_pred"]))}))
acc35<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(mean(data_tmp2[(data_tmp2$better_umi/data_tmp2$worse_umi)>x,"better_pred"]<data_tmp2[(data_tmp2$better_umi/data_tmp2$worse_umi)>x,"worse_pred"]))}))
acc15<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(mean(data_tmp3[(data_tmp3$better_umi/data_tmp3$worse_umi)>x,"better_pred"]<data_tmp3[(data_tmp3$better_umi/data_tmp3$worse_umi)>x,"worse_pred"]))}))
lower12<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_small,data=data_tmp1,type="lower"))
lower35<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_small,data=data_tmp2,type="lower"))
lower15<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_small,data=data_tmp3,type="lower"))
upper12<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_small,data=data_tmp1,type="upper"))
upper35<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_small,data=data_tmp2,type="upper"))
upper15<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_small,data=data_tmp3,type="upper"))
num12<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(length(data_tmp1[(data_tmp1$better_umi/data_tmp1$worse_umi)>x,"better_pred"]))}))
num35<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(length(data_tmp2[(data_tmp2$better_umi/data_tmp2$worse_umi)>x,"better_pred"]))}))
num15<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(length(data_tmp3[(data_tmp3$better_umi/data_tmp3$worse_umi)>x,"better_pred"]))}))

plot_data2_pmtnetv11<-data.frame(fold=rep(c(seq(2,6.5,0.125)),3), acc=c(acc12,acc35,acc15),lower=c(lower12,lower35,lower15),upper=c(upper12,upper35,upper15),num=c(num12,num35,num15),group=c(rep("1~3_AA",37),rep("4~5_AA",37),rep("1~5_AA",37)))

#ATM-TCR
all_comparison_all_info<-read.table("all_comparison_mean2.txt",header = T,sep = "\t",stringsAsFactors = F)
input<-read.table("input_data_continuous_better_worse_for_benchmark_software_have_vname.txt",header = T,sep="\t")
input<-input[input$TRA_V!="" & input$CDR3_ALPHA!="" & input$TRB_V!="" & input$CDR3_BETA!="",]
input<-input[!duplicated(input),]
input$TRA_V<-unlist(lapply(strsplit(input$TRA_V,"*",fixed = T),`[`,1))
input$TRB_V<-unlist(lapply(strsplit(input$TRB_V,"*",fixed = T),`[`,1))
input$TRA_V[input$TRA_V=="TRAV14/DV4"]<-"TRAV14DV4"
input$TRA_V[input$TRA_V=="TRAV38-2/DV8"]<-"TRAV38-2DV8"
input$TRA_V[input$TRA_V=="TRAV29/DV5"]<-"TRAV29DV5"
input$TRA_V[input$TRA_V=="TRAV36/DV7"]<-"TRAV36DV7"
input$TRA_V[input$TRA_V=="TRAV23/DV6"]<-"TRAV23DV6"

input$cdr3pep<-paste(input$CDR3_BETA,input$EPITOPE,sep = "_")
input$pmhc<-paste(unlist(lapply(lapply(strsplit(input$MHC_ALPHA,"[-:\\*]"),`[`,c(2,3,4)), paste,collapse="")),input$EPITOPE,sep = "_")
input$tcr<-paste(input$TRA_V,input$TRB_V,input$CDR3_ALPHA,input$CDR3_BETA,sep="_")
input$pmhctcr<-paste(input$pmhc,input$tcr,sep="_")

atm<-read.table("atm_tcr_results.csv",header = F,sep="\t")
atm<-atm[atm$V2 %in% input$CDR3_BETA,]
stopifnot(sum(!(atm$V1 %in% input$EPITOPE))==0)
atm<-atm[!duplicated(atm),]
atm$cdr3pep<-paste(atm$V2,atm$V1,sep = "_")
input_combine<-merge(input,atm,by="cdr3pep")
rownames(input_combine)<-input_combine$pmhctcr

all_comparison_all_info$better_pred<-NA
all_comparison_all_info$worse_pred<-NA
all_comparison_all_info$better_pred<-input_combine[paste(all_comparison_all_info$pMHC,all_comparison_all_info$Va,all_comparison_all_info$Vb,all_comparison_all_info$betterCDR3a,all_comparison_all_info$betterCDR3b,sep = "_"),"V5"]
all_comparison_all_info$worse_pred<-input_combine[paste(all_comparison_all_info$pMHC,all_comparison_all_info$Va,all_comparison_all_info$Vb,all_comparison_all_info$worseCDR3a,all_comparison_all_info$worseCDR3b,sep = "_"),"V5"]
all_comparison_all_info<-na.omit(all_comparison_all_info)

data_tmp1<-all_comparison_all_info[(all_comparison_all_info$dista + all_comparison_all_info$distb)>=1 & 
                                     (all_comparison_all_info$dista + all_comparison_all_info$distb)<=3 & 
                                     all_comparison_all_info$better_umi>min_umi_better & all_comparison_all_info$worse_umi>min_umi_worse,]
data_tmp2<-all_comparison_all_info[(all_comparison_all_info$dista + all_comparison_all_info$distb)>=4 & 
                                     (all_comparison_all_info$dista + all_comparison_all_info$distb)<=5 & 
                                     all_comparison_all_info$better_umi>min_umi_better & all_comparison_all_info$worse_umi>min_umi_worse,]
data_tmp3<-all_comparison_all_info[(all_comparison_all_info$dista + all_comparison_all_info$distb)>=1 & 
                                     (all_comparison_all_info$dista + all_comparison_all_info$distb)<=5 & 
                                     all_comparison_all_info$better_umi>min_umi_better & all_comparison_all_info$worse_umi>min_umi_worse,]

acc12<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(mean(data_tmp1[(data_tmp1$better_umi/data_tmp1$worse_umi)>x,"better_pred"]>data_tmp1[(data_tmp1$better_umi/data_tmp1$worse_umi)>x,"worse_pred"]))}))
acc35<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(mean(data_tmp2[(data_tmp2$better_umi/data_tmp2$worse_umi)>x,"better_pred"]>data_tmp2[(data_tmp2$better_umi/data_tmp2$worse_umi)>x,"worse_pred"]))}))
acc15<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(mean(data_tmp3[(data_tmp3$better_umi/data_tmp3$worse_umi)>x,"better_pred"]>data_tmp3[(data_tmp3$better_umi/data_tmp3$worse_umi)>x,"worse_pred"]))}))
lower12<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp1,type="lower"))
lower35<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp2,type="lower"))
lower15<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp3,type="lower"))
upper12<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp1,type="upper"))
upper35<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp2,type="upper"))
upper15<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp3,type="upper"))
num12<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(length(data_tmp1[(data_tmp1$better_umi/data_tmp1$worse_umi)>x,"better_pred"]))}))
num35<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(length(data_tmp2[(data_tmp2$better_umi/data_tmp2$worse_umi)>x,"better_pred"]))}))
num15<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(length(data_tmp3[(data_tmp3$better_umi/data_tmp3$worse_umi)>x,"better_pred"]))}))

plot_data2_atm<-data.frame(fold=rep(c(seq(2,6.5,0.125)),3), acc=c(acc12,acc35,acc15),lower=c(lower12,lower35,lower15),upper=c(upper12,upper35,upper15),num=c(num12,num35,num15),group=c(rep("1~3_AA",37),rep("4~5_AA",37),rep("1~5_AA",37)))
if(nrow(plot_data2_atm[!is.na(plot_data2_atm$acc),])<3){next}
rm(list=setdiff(ls(),c("plot_data2","umi_cutoff","plot_data2_imrex","plot_data2_atm","get_ci_est_large")))

# ergo 1
all_comparison_all_info<-read.table("all_comparison_mean2.txt",header = T,sep = "\t",stringsAsFactors = F)
all_comparison_all_info$better_pred<-NA
all_comparison_all_info$worse_pred<-NA
ergo<-read.table("ergo_ii_data_result_ae_mcpas.csv",header = T,sep = ",")
ergo<-ergo[!duplicated(ergo[,1:7]),]

ergo$TRAV<-unlist(lapply(strsplit(ergo$TRAV,"*",fixed = T),`[`,1))
ergo$TRBV<-unlist(lapply(strsplit(ergo$TRBV,"*",fixed = T),`[`,1))
ergo$pmhc<-paste(unlist(lapply(lapply(strsplit(ergo$MHC,"[-:\\*]"),`[`,c(2,3,4)), paste,collapse="")),ergo$Peptide,sep = "_")
ergo$tcr<-paste(ergo$TRAV,ergo$TRBV,ergo$TRA,ergo$TRB,sep="_")
ergo$pmhctcr<-paste(ergo$pmhc,ergo$tcr,sep="_")
rownames(ergo)<-ergo$pmhctcr

all_comparison_all_info$better_pred<-ergo[paste(all_comparison_all_info$pMHC,all_comparison_all_info$Va,all_comparison_all_info$Vb,all_comparison_all_info$betterCDR3a,all_comparison_all_info$betterCDR3b,sep = "_"),"Score"]
all_comparison_all_info$worse_pred<-ergo[paste(all_comparison_all_info$pMHC,all_comparison_all_info$Va,all_comparison_all_info$Vb,all_comparison_all_info$worseCDR3a,all_comparison_all_info$worseCDR3b,sep = "_"),"Score"]
all_comparison_all_info<-na.omit(all_comparison_all_info)

data_tmp1<-all_comparison_all_info[(all_comparison_all_info$dista + all_comparison_all_info$distb)>=1 & 
                                     (all_comparison_all_info$dista + all_comparison_all_info$distb)<=3 & 
                                     all_comparison_all_info$better_umi>min_umi_better & all_comparison_all_info$worse_umi>min_umi_worse,]
data_tmp2<-all_comparison_all_info[(all_comparison_all_info$dista + all_comparison_all_info$distb)>=4 & 
                                     (all_comparison_all_info$dista + all_comparison_all_info$distb)<=5 & 
                                     all_comparison_all_info$better_umi>min_umi_better & all_comparison_all_info$worse_umi>min_umi_worse,]
data_tmp3<-all_comparison_all_info[(all_comparison_all_info$dista + all_comparison_all_info$distb)>=1 & 
                                     (all_comparison_all_info$dista + all_comparison_all_info$distb)<=5 & 
                                     all_comparison_all_info$better_umi>min_umi_better & all_comparison_all_info$worse_umi>min_umi_worse,]

acc12<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(mean(data_tmp1[(data_tmp1$better_umi/data_tmp1$worse_umi)>x,"better_pred"]>data_tmp1[(data_tmp1$better_umi/data_tmp1$worse_umi)>x,"worse_pred"]))}))
acc35<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(mean(data_tmp2[(data_tmp2$better_umi/data_tmp2$worse_umi)>x,"better_pred"]>data_tmp2[(data_tmp2$better_umi/data_tmp2$worse_umi)>x,"worse_pred"]))}))
acc15<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(mean(data_tmp3[(data_tmp3$better_umi/data_tmp3$worse_umi)>x,"better_pred"]>data_tmp3[(data_tmp3$better_umi/data_tmp3$worse_umi)>x,"worse_pred"]))}))
lower12<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp1,type="lower"))
lower35<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp2,type="lower"))
lower15<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp3,type="lower"))
upper12<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp1,type="upper"))
upper35<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp2,type="upper"))
upper15<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp3,type="upper"))
num12<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(length(data_tmp1[(data_tmp1$better_umi/data_tmp1$worse_umi)>x,"better_pred"]))}))
num35<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(length(data_tmp2[(data_tmp2$better_umi/data_tmp2$worse_umi)>x,"better_pred"]))}))
num15<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(length(data_tmp3[(data_tmp3$better_umi/data_tmp3$worse_umi)>x,"better_pred"]))}))

plot_data2_ergo1<-data.frame(fold=rep(c(seq(2,6.5,0.125)),3), acc=c(acc12,acc35,acc15),lower=c(lower12,lower35,lower15),upper=c(upper12,upper35,upper15),num=c(num12,num35,num15),group=c(rep("1~3_AA",37),rep("4~5_AA",37),rep("1~5_AA",37)))
rm(list=setdiff(ls(),c("plot_data2","umi_cutoff","plot_data2_imrex","plot_data2_atm","get_ci_est_large","plot_data2_ergo1")))

# ergo 2
all_comparison_all_info<-read.table("all_comparison_mean2.txt",header = T,sep = "\t",stringsAsFactors = F)
all_comparison_all_info$better_pred<-NA
all_comparison_all_info$worse_pred<-NA

ergo<-read.table("ergo_ii_data_result_ae_vdjdb.csv",header = T,sep = ",")
ergo<-ergo[!duplicated(ergo[,1:7]),]

ergo$TRAV<-unlist(lapply(strsplit(ergo$TRAV,"*",fixed = T),`[`,1))
ergo$TRBV<-unlist(lapply(strsplit(ergo$TRBV,"*",fixed = T),`[`,1))
ergo$pmhc<-paste(unlist(lapply(lapply(strsplit(ergo$MHC,"[-:\\*]"),`[`,c(2,3,4)), paste,collapse="")),ergo$Peptide,sep = "_")
ergo$tcr<-paste(ergo$TRAV,ergo$TRBV,ergo$TRA,ergo$TRB,sep="_")
ergo$pmhctcr<-paste(ergo$pmhc,ergo$tcr,sep="_")
rownames(ergo)<-ergo$pmhctcr

all_comparison_all_info$better_pred<-ergo[paste(all_comparison_all_info$pMHC,all_comparison_all_info$Va,all_comparison_all_info$Vb,all_comparison_all_info$betterCDR3a,all_comparison_all_info$betterCDR3b,sep = "_"),"Score"]
all_comparison_all_info$worse_pred<-ergo[paste(all_comparison_all_info$pMHC,all_comparison_all_info$Va,all_comparison_all_info$Vb,all_comparison_all_info$worseCDR3a,all_comparison_all_info$worseCDR3b,sep = "_"),"Score"]
all_comparison_all_info<-na.omit(all_comparison_all_info)

data_tmp1<-all_comparison_all_info[(all_comparison_all_info$dista + all_comparison_all_info$distb)>=1 & 
                                     (all_comparison_all_info$dista + all_comparison_all_info$distb)<=3 & 
                                     all_comparison_all_info$better_umi>min_umi_better & all_comparison_all_info$worse_umi>min_umi_worse,]
data_tmp2<-all_comparison_all_info[(all_comparison_all_info$dista + all_comparison_all_info$distb)>=4 & 
                                     (all_comparison_all_info$dista + all_comparison_all_info$distb)<=5 & 
                                     all_comparison_all_info$better_umi>min_umi_better & all_comparison_all_info$worse_umi>min_umi_worse,]
data_tmp3<-all_comparison_all_info[(all_comparison_all_info$dista + all_comparison_all_info$distb)>=1 & 
                                     (all_comparison_all_info$dista + all_comparison_all_info$distb)<=5 & 
                                     all_comparison_all_info$better_umi>min_umi_better & all_comparison_all_info$worse_umi>min_umi_worse,]

acc12<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(mean(data_tmp1[(data_tmp1$better_umi/data_tmp1$worse_umi)>x,"better_pred"]>data_tmp1[(data_tmp1$better_umi/data_tmp1$worse_umi)>x,"worse_pred"]))}))
acc35<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(mean(data_tmp2[(data_tmp2$better_umi/data_tmp2$worse_umi)>x,"better_pred"]>data_tmp2[(data_tmp2$better_umi/data_tmp2$worse_umi)>x,"worse_pred"]))}))
acc15<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(mean(data_tmp3[(data_tmp3$better_umi/data_tmp3$worse_umi)>x,"better_pred"]>data_tmp3[(data_tmp3$better_umi/data_tmp3$worse_umi)>x,"worse_pred"]))}))
lower12<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp1,type="lower"))
lower35<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp2,type="lower"))
lower15<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp3,type="lower"))
upper12<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp1,type="upper"))
upper35<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp2,type="upper"))
upper15<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp3,type="upper"))
num12<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(length(data_tmp1[(data_tmp1$better_umi/data_tmp1$worse_umi)>x,"better_pred"]))}))
num35<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(length(data_tmp2[(data_tmp2$better_umi/data_tmp2$worse_umi)>x,"better_pred"]))}))
num15<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(length(data_tmp3[(data_tmp3$better_umi/data_tmp3$worse_umi)>x,"better_pred"]))}))

plot_data2_ergo2<-data.frame(fold=rep(c(seq(2,6.5,0.125)),3), acc=c(acc12,acc35,acc15),lower=c(lower12,lower35,lower15),upper=c(upper12,upper35,upper15),num=c(num12,num35,num15),group=c(rep("1~3_AA",37),rep("4~5_AA",37),rep("1~5_AA",37)))
rm(list=setdiff(ls(),c("plot_data2","umi_cutoff","plot_data2_imrex","plot_data2_atm","get_ci_est_large","plot_data2_ergo1","plot_data2_ergo2")))

# ergo 3
all_comparison_all_info<-read.table("all_comparison_mean2.txt",header = T,sep = "\t",stringsAsFactors = F)
all_comparison_all_info$better_pred<-NA
all_comparison_all_info$worse_pred<-NA

ergo<-read.table("ergo_ii_data_result_lstm_mcpas.csv",header = T,sep = ",")
ergo<-ergo[!duplicated(ergo[,1:7]),]

ergo$TRAV<-unlist(lapply(strsplit(ergo$TRAV,"*",fixed = T),`[`,1))
ergo$TRBV<-unlist(lapply(strsplit(ergo$TRBV,"*",fixed = T),`[`,1))
ergo$pmhc<-paste(unlist(lapply(lapply(strsplit(ergo$MHC,"[-:\\*]"),`[`,c(2,3,4)), paste,collapse="")),ergo$Peptide,sep = "_")
ergo$tcr<-paste(ergo$TRAV,ergo$TRBV,ergo$TRA,ergo$TRB,sep="_")
ergo$pmhctcr<-paste(ergo$pmhc,ergo$tcr,sep="_")
rownames(ergo)<-ergo$pmhctcr

all_comparison_all_info$better_pred<-ergo[paste(all_comparison_all_info$pMHC,all_comparison_all_info$Va,all_comparison_all_info$Vb,all_comparison_all_info$betterCDR3a,all_comparison_all_info$betterCDR3b,sep = "_"),"Score"]
all_comparison_all_info$worse_pred<-ergo[paste(all_comparison_all_info$pMHC,all_comparison_all_info$Va,all_comparison_all_info$Vb,all_comparison_all_info$worseCDR3a,all_comparison_all_info$worseCDR3b,sep = "_"),"Score"]
all_comparison_all_info<-na.omit(all_comparison_all_info)

data_tmp1<-all_comparison_all_info[(all_comparison_all_info$dista + all_comparison_all_info$distb)>=1 & 
                                     (all_comparison_all_info$dista + all_comparison_all_info$distb)<=3 & 
                                     all_comparison_all_info$better_umi>min_umi_better & all_comparison_all_info$worse_umi>min_umi_worse,]
data_tmp2<-all_comparison_all_info[(all_comparison_all_info$dista + all_comparison_all_info$distb)>=4 & 
                                     (all_comparison_all_info$dista + all_comparison_all_info$distb)<=5 & 
                                     all_comparison_all_info$better_umi>min_umi_better & all_comparison_all_info$worse_umi>min_umi_worse,]
data_tmp3<-all_comparison_all_info[(all_comparison_all_info$dista + all_comparison_all_info$distb)>=1 & 
                                     (all_comparison_all_info$dista + all_comparison_all_info$distb)<=5 & 
                                     all_comparison_all_info$better_umi>min_umi_better & all_comparison_all_info$worse_umi>min_umi_worse,]

acc12<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(mean(data_tmp1[(data_tmp1$better_umi/data_tmp1$worse_umi)>x,"better_pred"]>data_tmp1[(data_tmp1$better_umi/data_tmp1$worse_umi)>x,"worse_pred"]))}))
acc35<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(mean(data_tmp2[(data_tmp2$better_umi/data_tmp2$worse_umi)>x,"better_pred"]>data_tmp2[(data_tmp2$better_umi/data_tmp2$worse_umi)>x,"worse_pred"]))}))
acc15<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(mean(data_tmp3[(data_tmp3$better_umi/data_tmp3$worse_umi)>x,"better_pred"]>data_tmp3[(data_tmp3$better_umi/data_tmp3$worse_umi)>x,"worse_pred"]))}))
lower12<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp1,type="lower"))
lower35<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp2,type="lower"))
lower15<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp3,type="lower"))
upper12<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp1,type="upper"))
upper35<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp2,type="upper"))
upper15<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp3,type="upper"))
num12<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(length(data_tmp1[(data_tmp1$better_umi/data_tmp1$worse_umi)>x,"better_pred"]))}))
num35<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(length(data_tmp2[(data_tmp2$better_umi/data_tmp2$worse_umi)>x,"better_pred"]))}))
num15<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(length(data_tmp3[(data_tmp3$better_umi/data_tmp3$worse_umi)>x,"better_pred"]))}))

plot_data2_ergo3<-data.frame(fold=rep(c(seq(2,6.5,0.125)),3), acc=c(acc12,acc35,acc15),lower=c(lower12,lower35,lower15),upper=c(upper12,upper35,upper15),num=c(num12,num35,num15),group=c(rep("1~3_AA",37),rep("4~5_AA",37),rep("1~5_AA",37)))

rm(list=setdiff(ls(),c("plot_data2","umi_cutoff","plot_data2_imrex","plot_data2_atm","get_ci_est_large","plot_data2_ergo1","plot_data2_ergo2","plot_data2_ergo3")))

# ergo4
all_comparison_all_info<-read.table("all_comparison_mean2.txt",header = T,sep = "\t",stringsAsFactors = F)
all_comparison_all_info$better_pred<-NA
all_comparison_all_info$worse_pred<-NA

ergo<-read.table("ergo_ii_data_result_lstm_vdjdb.csv",header = T,sep = ",")
ergo<-ergo[!duplicated(ergo[,1:7]),]

ergo$TRAV<-unlist(lapply(strsplit(ergo$TRAV,"*",fixed = T),`[`,1))
ergo$TRBV<-unlist(lapply(strsplit(ergo$TRBV,"*",fixed = T),`[`,1))
ergo$pmhc<-paste(unlist(lapply(lapply(strsplit(ergo$MHC,"[-:\\*]"),`[`,c(2,3,4)), paste,collapse="")),ergo$Peptide,sep = "_")
ergo$tcr<-paste(ergo$TRAV,ergo$TRBV,ergo$TRA,ergo$TRB,sep="_")
ergo$pmhctcr<-paste(ergo$pmhc,ergo$tcr,sep="_")
rownames(ergo)<-ergo$pmhctcr

all_comparison_all_info$better_pred<-ergo[paste(all_comparison_all_info$pMHC,all_comparison_all_info$Va,all_comparison_all_info$Vb,all_comparison_all_info$betterCDR3a,all_comparison_all_info$betterCDR3b,sep = "_"),"Score"]
all_comparison_all_info$worse_pred<-ergo[paste(all_comparison_all_info$pMHC,all_comparison_all_info$Va,all_comparison_all_info$Vb,all_comparison_all_info$worseCDR3a,all_comparison_all_info$worseCDR3b,sep = "_"),"Score"]
all_comparison_all_info<-na.omit(all_comparison_all_info)

data_tmp1<-all_comparison_all_info[(all_comparison_all_info$dista + all_comparison_all_info$distb)>=1 & 
                                     (all_comparison_all_info$dista + all_comparison_all_info$distb)<=3 & 
                                     all_comparison_all_info$better_umi>min_umi_better & all_comparison_all_info$worse_umi>min_umi_worse,]
data_tmp2<-all_comparison_all_info[(all_comparison_all_info$dista + all_comparison_all_info$distb)>=4 & 
                                     (all_comparison_all_info$dista + all_comparison_all_info$distb)<=5 & 
                                     all_comparison_all_info$better_umi>min_umi_better & all_comparison_all_info$worse_umi>min_umi_worse,]
data_tmp3<-all_comparison_all_info[(all_comparison_all_info$dista + all_comparison_all_info$distb)>=1 & 
                                     (all_comparison_all_info$dista + all_comparison_all_info$distb)<=5 & 
                                     all_comparison_all_info$better_umi>min_umi_better & all_comparison_all_info$worse_umi>min_umi_worse,]

acc12<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(mean(data_tmp1[(data_tmp1$better_umi/data_tmp1$worse_umi)>x,"better_pred"]>data_tmp1[(data_tmp1$better_umi/data_tmp1$worse_umi)>x,"worse_pred"]))}))
acc35<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(mean(data_tmp2[(data_tmp2$better_umi/data_tmp2$worse_umi)>x,"better_pred"]>data_tmp2[(data_tmp2$better_umi/data_tmp2$worse_umi)>x,"worse_pred"]))}))
acc15<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(mean(data_tmp3[(data_tmp3$better_umi/data_tmp3$worse_umi)>x,"better_pred"]>data_tmp3[(data_tmp3$better_umi/data_tmp3$worse_umi)>x,"worse_pred"]))}))
lower12<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp1,type="lower"))
lower35<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp2,type="lower"))
lower15<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp3,type="lower"))
upper12<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp1,type="upper"))
upper35<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp2,type="upper"))
upper15<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp3,type="upper"))
num12<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(length(data_tmp1[(data_tmp1$better_umi/data_tmp1$worse_umi)>x,"better_pred"]))}))
num35<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(length(data_tmp2[(data_tmp2$better_umi/data_tmp2$worse_umi)>x,"better_pred"]))}))
num15<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(length(data_tmp3[(data_tmp3$better_umi/data_tmp3$worse_umi)>x,"better_pred"]))}))

plot_data2_ergo4<-data.frame(fold=rep(c(seq(2,6.5,0.125)),3), acc=c(acc12,acc35,acc15),lower=c(lower12,lower35,lower15),upper=c(upper12,upper35,upper15),num=c(num12,num35,num15),group=c(rep("1~3_AA",37),rep("4~5_AA",37),rep("1~5_AA",37)))

rm(list=setdiff(ls(),c("plot_data2","umi_cutoff","plot_data2_imrex","plot_data2_atm","get_ci_est_large","plot_data2_ergo1","plot_data2_ergo2","plot_data2_ergo3","plot_data2_ergo4")))

#Imrex
get_ci_est_large<-function(fold,data,type){
  if(sum((data$better_umi/data$worse_umi)>fold)>1){
    tmp<-confintr::ci_proportion(x=as.numeric(data[(data$better_umi/data$worse_umi)>fold,"better_pred"]>data[(data$better_umi/data$worse_umi)>fold,"worse_pred"]),type=CI_TYPE)
    if(type=="lower"){
      return(tmp$interval[1])
    } else if(type=="upper"){
      return(tmp$interval[2])
    }else{
      stop("Type is wrong!")
    }
  }else{
    return(NA)
  }
}
all_comparison_all_info<-read.table("all_comparison_mean2.txt",header = T,sep = "\t",stringsAsFactors = F)

input<-read.table("input_data_continuous_better_worse_for_benchmark_software_have_vname.txt",header = T,sep="\t")
input<-input[input$TRA_V!="" & input$CDR3_ALPHA!="" & input$TRB_V!="" & input$CDR3_BETA!="",]
input<-input[!duplicated(input),]
input$TRA_V<-unlist(lapply(strsplit(input$TRA_V,"*",fixed = T),`[`,1))
input$TRB_V<-unlist(lapply(strsplit(input$TRB_V,"*",fixed = T),`[`,1))
input$TRA_V[input$TRA_V=="TRAV14/DV4"]<-"TRAV14DV4"
input$TRA_V[input$TRA_V=="TRAV38-2/DV8"]<-"TRAV38-2DV8"
input$TRA_V[input$TRA_V=="TRAV29/DV5"]<-"TRAV29DV5"
input$TRA_V[input$TRA_V=="TRAV36/DV7"]<-"TRAV36DV7"
input$TRA_V[input$TRA_V=="TRAV23/DV6"]<-"TRAV23DV6"

input$cdr3pep<-paste(input$CDR3_BETA,input$EPITOPE,sep = "_")
input$pmhc<-paste(unlist(lapply(lapply(strsplit(input$MHC_ALPHA,"[-:\\*]"),`[`,c(2,3,4)), paste,collapse="")),input$EPITOPE,sep = "_")
input$tcr<-paste(input$TRA_V,input$TRB_V,input$CDR3_ALPHA,input$CDR3_BETA,sep="_")
input$pmhctcr<-paste(input$pmhc,input$tcr,sep="_")

imrex<-read.table("imrex_data_results.csv",header = T,sep = ",")
imrex<-imrex[!duplicated(imrex),]
stopifnot(sum(!(imrex$cdr3 %in% input$CDR3_BETA))==0)
stopifnot(sum(!(imrex$antigen.epitope %in% input$EPITOPE))==0)

imrex$cdr3pep<-paste(imrex$cdr3,imrex$antigen.epitope,sep = "_")
input_combine<-merge(input,imrex,by="cdr3pep")
rownames(input_combine)<-input_combine$pmhctcr

all_comparison_all_info$better_pred<-NA
all_comparison_all_info$worse_pred<-NA
all_comparison_all_info$better_pred<-input_combine[paste(all_comparison_all_info$pMHC,all_comparison_all_info$Va,all_comparison_all_info$Vb,all_comparison_all_info$betterCDR3a,all_comparison_all_info$betterCDR3b,sep = "_"),"prediction_score"]
all_comparison_all_info$worse_pred<-input_combine[paste(all_comparison_all_info$pMHC,all_comparison_all_info$Va,all_comparison_all_info$Vb,all_comparison_all_info$worseCDR3a,all_comparison_all_info$worseCDR3b,sep = "_"),"prediction_score"]
all_comparison_all_info<-na.omit(all_comparison_all_info)

data_tmp1<-all_comparison_all_info[(all_comparison_all_info$dista + all_comparison_all_info$distb)>=1 & 
                                     (all_comparison_all_info$dista + all_comparison_all_info$distb)<=3 & 
                                     all_comparison_all_info$better_umi>min_umi_better & all_comparison_all_info$worse_umi>min_umi_worse,]
data_tmp2<-all_comparison_all_info[(all_comparison_all_info$dista + all_comparison_all_info$distb)>=4 & 
                                     (all_comparison_all_info$dista + all_comparison_all_info$distb)<=5 & 
                                     all_comparison_all_info$better_umi>min_umi_better & all_comparison_all_info$worse_umi>min_umi_worse,]
data_tmp3<-all_comparison_all_info[(all_comparison_all_info$dista + all_comparison_all_info$distb)>=1 & 
                                     (all_comparison_all_info$dista + all_comparison_all_info$distb)<=5 & 
                                     all_comparison_all_info$better_umi>min_umi_better & all_comparison_all_info$worse_umi>min_umi_worse,]

acc12<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(mean(data_tmp1[(data_tmp1$better_umi/data_tmp1$worse_umi)>x,"better_pred"]>data_tmp1[(data_tmp1$better_umi/data_tmp1$worse_umi)>x,"worse_pred"]))}))
acc35<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(mean(data_tmp2[(data_tmp2$better_umi/data_tmp2$worse_umi)>x,"better_pred"]>data_tmp2[(data_tmp2$better_umi/data_tmp2$worse_umi)>x,"worse_pred"]))}))
acc15<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(mean(data_tmp3[(data_tmp3$better_umi/data_tmp3$worse_umi)>x,"better_pred"]>data_tmp3[(data_tmp3$better_umi/data_tmp3$worse_umi)>x,"worse_pred"]))}))
lower12<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp1,type="lower"))
lower35<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp2,type="lower"))
lower15<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp3,type="lower"))
upper12<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp1,type="upper"))
upper35<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp2,type="upper"))
upper15<-unlist(lapply(as.list(seq(2,6.5,0.125)),get_ci_est_large,data=data_tmp3,type="upper"))
num12<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(length(data_tmp1[(data_tmp1$better_umi/data_tmp1$worse_umi)>x,"better_pred"]))}))
num35<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(length(data_tmp2[(data_tmp2$better_umi/data_tmp2$worse_umi)>x,"better_pred"]))}))
num15<-unlist(lapply(as.list(seq(2,6.5,0.125)),FUN = function(x){return(length(data_tmp3[(data_tmp3$better_umi/data_tmp3$worse_umi)>x,"better_pred"]))}))

plot_data2_imrex<-data.frame(fold=rep(c(seq(2,6.5,0.125)),3), acc=c(acc12,acc35,acc15),lower=c(lower12,lower35,lower15),upper=c(upper12,upper35,upper15),num=c(num12,num35,num15),group=c(rep("1~3_AA",37),rep("4~5_AA",37),rep("1~5_AA",37)))
if(nrow(plot_data2[!is.na(plot_data2_imrex$acc),])<3){next}
rm(list=setdiff(ls(),c("plot_data2","umi_cutoff","plot_data2_imrex","get_ci_est_large")))

plot_data2$method<-rep("pMTnet Omni")
plot_data2_atm$method<-rep("ATM-TCR")
plot_data2_ergo1$method<-rep("ERGO-II:AE+MCPAS")
plot_data2_ergo2$method<-rep("ERGO-II:AE+VDJDB")
plot_data2_ergo3$method<-rep("ERGO-II:LSTM+MCPAS")
plot_data2_ergo4$method<-rep("ERGO-II:LSTM+VDJDB")
plot_data2_imrex$method<-rep("ImRex")
plot_data2_pmtnetv11$method<-rep("pMTnet")
plot_data2_teim1$method<-rep("TEIM")

text_size<-20
pd <- position_dodge(0.1)

plot_all_methods<-data.table::rbindlist(list(plot_data2[plot_data2$group=="1~5_AA",],
                                             plot_data2_atm[plot_data2_atm$group=="1~5_AA",],
                                             plot_data2_ergo1[plot_data2_ergo1$group=="1~5_AA",],
                                             plot_data2_ergo2[plot_data2_ergo2$group=="1~5_AA",],
                                             plot_data2_ergo3[plot_data2_ergo3$group=="1~5_AA",],
                                             plot_data2_ergo4[plot_data2_ergo4$group=="1~5_AA",],
                                             plot_data2_imrex[plot_data2_imrex$group=="1~5_AA",],
                                             plot_data2_pmtnetv11[plot_data2_pmtnetv11$group=="1~5_AA",],
                                             plot_data2_teim1[plot_data2_teim1$group=="1~5_AA",]))
plot_all_methods$method<-factor(plot_all_methods$method,levels = c('pMTnet Omni','pMTnet',"TEIM",'ATM-TCR','ImRex','ERGO-II:AE+VDJDB','ERGO-II:AE+MCPAS','ERGO-II:LSTM+VDJDB','ERGO-II:LSTM+MCPAS'))

ggplot(plot_all_methods[plot_all_methods$fold %in% c(4,5,6),],aes(x=fold,y=acc,color=method))+ geom_line(linewidth=2,position = pd)+ geom_point(size=2.5,position = pd)+
  theme_bw()+ylim(c(0,1))+
  scale_color_manual(values = c("tomato2","hotpink1","mediumpurple1","khaki","tan2","lightyellow3","lemonchiffon2","lightsteelblue2","lightskyblue2"))+
  xlab("umi fold")+ylab("Percent of correct prediction")+ggtitle("1~5AA")+guides(color=guide_legend(ncol=2))+
  theme(axis.title.y = element_text(size=text_size),
        axis.title.x = element_text(size=text_size),
        axis.text=element_text(size=text_size),
        legend.position=c(0.5,0.16),legend.text = element_text(size=text_size),legend.background = element_blank(),legend.title = element_blank())

plot_all_methods<-data.table::rbindlist(list(plot_data2[plot_data2$group=="1~3_AA",],
                                             plot_data2_atm[plot_data2_atm$group=="1~3_AA",],
                                             plot_data2_ergo1[plot_data2_ergo1$group=="1~3_AA",],
                                             plot_data2_ergo2[plot_data2_ergo2$group=="1~3_AA",],
                                             plot_data2_ergo3[plot_data2_ergo3$group=="1~3_AA",],
                                             plot_data2_ergo4[plot_data2_ergo4$group=="1~3_AA",],
                                             plot_data2_imrex[plot_data2_imrex$group=="1~3_AA",],
                                             plot_data2_pmtnetv11[plot_data2_pmtnetv11$group=="1~3_AA",],
                                             plot_data2_teim1[plot_data2_teim1$group=="1~3_AA",]))
plot_all_methods$method<-factor(plot_all_methods$method,levels = c('pMTnet Omni','pMTnet',"TEIM",'ATM-TCR','ImRex','ERGO-II:AE+VDJDB','ERGO-II:AE+MCPAS','ERGO-II:LSTM+VDJDB','ERGO-II:LSTM+MCPAS'))

ggplot(plot_all_methods[plot_all_methods$fold>=4 & plot_all_methods$fold<=6,],aes(x=fold,y=acc,color=method))+ geom_line(linewidth=2,position = pd)+ geom_point(size=2.5,position = pd)+
  theme_bw()+ylim(c(0,1))+
  scale_color_manual(values = c("tomato2","hotpink1","mediumpurple1","khaki","tan2","lightyellow3","lemonchiffon2","lightsteelblue2","lightskyblue2"))+
  xlab("umi fold")+ylab("Percent of correct prediction")+ggtitle("1~3AA")+guides(color=guide_legend(ncol=2))+
  theme(axis.title.y = element_text(size=text_size),
        axis.title.x = element_text(size=text_size),
        axis.text=element_text(size=text_size),
        legend.position=c(0.5,0.16),legend.text = element_text(size=text_size),legend.background = element_blank(),legend.title = element_blank())

plot_all_methods<-data.table::rbindlist(list(plot_data2[plot_data2$group=="4~5_AA",],
                                             plot_data2_atm[plot_data2_atm$group=="4~5_AA",],
                                             plot_data2_ergo1[plot_data2_ergo1$group=="4~5_AA",],
                                             plot_data2_ergo2[plot_data2_ergo2$group=="4~5_AA",],
                                             plot_data2_ergo3[plot_data2_ergo3$group=="4~5_AA",],
                                             plot_data2_ergo4[plot_data2_ergo4$group=="4~5_AA",],
                                             plot_data2_imrex[plot_data2_imrex$group=="4~5_AA",],
                                             plot_data2_pmtnetv11[plot_data2_pmtnetv11$group=="4~5_AA",],
                                             plot_data2_teim1[plot_data2_teim1$group=="4~5_AA",]))
plot_all_methods$method<-factor(plot_all_methods$method,levels = c('pMTnet Omni','pMTnet',"TEIM",'ATM-TCR','ImRex','ERGO-II:AE+VDJDB','ERGO-II:AE+MCPAS','ERGO-II:LSTM+VDJDB','ERGO-II:LSTM+MCPAS'))
ggplot(plot_all_methods[plot_all_methods$fold>=4 & plot_all_methods$fold<=6,],aes(x=fold,y=acc,color=method))+ geom_line(linewidth=2,position = pd)+ geom_point(size=2.5,position = pd)+
  theme_bw()+ylim(c(0,1))+
  scale_color_manual(values = c("tomato2","hotpink1","mediumpurple1","khaki","tan2","lightyellow3","lemonchiffon2","lightsteelblue2","lightskyblue2"))+
  xlab("umi fold")+ylab("Percent of correct prediction")+ggtitle("4~5AA")+guides(color=guide_legend(ncol=2))+
  theme(axis.title.y = element_text(size=text_size),
        axis.title.x = element_text(size=text_size),
        axis.text=element_text(size=text_size),
        legend.position=c(0.5,0.16),legend.text = element_text(size=text_size),legend.background = element_blank(),legend.title = element_blank())

#Fig. 3 PRAME
prame<-read.table('valdata_prame.csv',header = T,sep = ',')
kras1<-read.table('input_top1_one.tsv',header = T,sep='\t')
kras2<-read.table('input_top2.tsv',header = T,sep='\t')
kras3<-read.table('input_uniqueTRATRB_part.tsv',header = T,sep='\t')
kras<-rbind(kras1[,c("peptide","mhca","mhcb","va","cdr3a","vb","cdr3b","TCR_SPECIES","pMHC_SPECIES","vaseq","vbseq")],
            kras2[,c("peptide","mhca","mhcb","va","cdr3a","vb","cdr3b","TCR_SPECIES","pMHC_SPECIES","vaseq","vbseq")],
            kras3[,c("peptide","mhca","mhcb","va","cdr3a","vb","cdr3b","TCR_SPECIES","pMHC_SPECIES","vaseq","vbseq")])
rm(list=c('kras1','kras2','kras3'))

pmtnetv1<-data.frame(CDR3=c(prame$cdr3b,kras$cdr3b),Antigen=c(prame$peptide,kras$peptide),HLA=c(prame$mhca,kras$mhca))
teim<-data.frame(cdr3=c(prame$cdr3b,kras$cdr3b),epitope=c(prame$peptide,kras$peptide))
atmtcr<-data.frame(c(prame$peptide,kras$peptide),c(prame$cdr3b,kras$cdr3b),rep(0)) #no header
imrex<-data.frame(cdr3=c(prame$cdr3b,kras$cdr3b),antigen.epitope=c(prame$peptide,kras$peptide)) #sep=；
ergo2<-data.frame(TRA=c(prame$cdr3a,kras$cdr3a),TRB=c(prame$cdr3b,kras$cdr3b),Peptide=c(prame$peptide,kras$peptide),
                  MHC=paste('HLA-',unlist(lapply(strsplit(c(prame$mhca,kras$mhca),split = ":"),'[[',1)),sep = "")) # write  'TRA', 'TRB', 'TRAV', 'TRAJ', 'TRBV', 'TRBJ', 'T_Cell_Type', 'Peptide' and 'MHC', then change T_Cell_Type to T-Cell-Type
write.table(pmtnetv1,'input_pmtnetv1.csv',col.names = T, row.names = F,sep = ',',quote = F)
write.table(teim,'input_teim.csv',col.names = T, row.names = F,sep = ',',quote = F)
write.table(atmtcr,'input_atmtcr.csv',col.names = F, row.names = F,sep = ',',quote = F)
write.table(imrex,'input_imrex.csv',col.names = T, row.names = F,sep = ';',quote = F)
write.table(ergo2,'input_ergo2.csv',col.names = T, row.names = F,sep = ',',quote = F)

pmtnetv1<-read.table('pmtnetv1_prediction_prame_kras.csv',header = T,sep=',')
teim<-read.table('out_teim.csv',header = T,sep=',')
atm<-read.table('out_atmtcr.csv',header = F,sep='\t')
imrex<-read.table('output_imrex.csv',header = T,sep=',')
ergo2ae.mac<-read.table('ergo2-ae-mcpas.csv',header = T,sep=',')
ergo2ae.vdj<-read.table('ergo2-ae-vdjdb.csv',header = T,sep=',')
ergo2ls.mac<-read.table('ergo2-lstm-mcpas.csv',header = T,sep=',')
ergo2ls.vdj<-read.table('ergo2-lstm-vdjdb.csv',header = T,sep=',')
pmtnetv2prame<-read.table("results_numTrain0_numSeedSplit10_pre.csv",header = T,sep = ',')
pmtnetv2kras<-read.table('results_pre_alex_subset.csv',header = T,sep = ',')
pmtnetv2<-rbind(pmtnetv2prame[,c("peptide","cdr3a","cdr3b","avg_percentile_rank","label")],pmtnetv2kras[,c("peptide","cdr3a","cdr3b","avg_percentile_rank","label")])

epact.prame<-read.table('preds_epact_origin.csv',header = T,sep=',')
nettcr.prame<-read.table('prame_nettcr_result.csv',header = T,sep=',')
lines <- readLines("ranking-tulip-origin.jsonl")
tulip.prame <- lapply(lines, fromJSON)
tulip.prame <- data.frame(peptide=rep(tulip.prame[[1]]$peptide),cdr3a=tulip.prame[[1]]$CDR3a,cdr3b=tulip.prame[[1]]$CDR3b,rank=tulip.prame[[1]]$rank)

prame<-read.table('valdata_prame.csv',header = T,sep = ',')
kras1<-read.table('input_top1_one.tsv',header = T,sep='\t')
kras2<-read.table('input_top2.tsv',header = T,sep='\t')
kras3<-read.table('input_uniqueTRATRB_part.tsv',header = T,sep='\t')
kras<-rbind(kras1[,c("peptide","mhca","mhcb","va","cdr3a","vb","cdr3b","TCR_SPECIES","pMHC_SPECIES","vaseq","vbseq")],
            kras2[,c("peptide","mhca","mhcb","va","cdr3a","vb","cdr3b","TCR_SPECIES","pMHC_SPECIES","vaseq","vbseq")],
            kras3[,c("peptide","mhca","mhcb","va","cdr3a","vb","cdr3b","TCR_SPECIES","pMHC_SPECIES","vaseq","vbseq")])
kras$label<-c(rep(1,6),rep(0,60))

rawdata<-data.frame(cdr3a=c(prame$cdr3a,kras$cdr3a),cdr3b=c(prame$cdr3b,kras$cdr3b),peptide=c(prame$peptide,kras$peptide),mhc=c(prame$mhca,kras$mhca),label=c(prame$label,kras$label),ID=c(prame$ID,rep('non',66)))
rm(list=c('prame','kras3','kras','pmtnetv2prame','pmtnetv2kras'))

assign_label<-function(df,num){
  df$label<-0
  df$label[df[,num] %in% c(kras1$cdr3b,kras2$cdr3b)]<-1
  df$label[df[,num]=="CASSWWTGGASPIRF"]<-rawdata$label[rawdata$cdr3b=="CASSWWTGGASPIRF"]
  df$label[df[,num]=="CASSWWTGGSAPISF"]<-rawdata$label[rawdata$cdr3b=="CASSWWTGGSAPISF"]
  df$label[df[,num]=="CASSWWTGGSSPISF"]<-rawdata$label[rawdata$cdr3b=="CASSWWTGGSSPISF"]
  df$label[df[,num]=="CASSWWTGGSAPIRF"]<-rawdata$label[rawdata$cdr3b=="CASSWWTGGSAPIRF"]
  df$label[df[,num]=="CASSPWTGGSAPIRF"]<-rawdata$label[rawdata$cdr3b=="CASSPWTGGSAPIRF"]
  df$label[df[,num]=="CASSWWTSGSAPIRF"]<-rawdata$label[rawdata$cdr3b=="CASSWWTSGSAPIRF"]
  df$label[df[,num]=="CASSWWTGGSAEIRF"]<-rawdata$label[rawdata$cdr3b=="CASSWWTGGSAEIRF"]
  df$label[df[,num]=="CASSWWTGGSAPIYF"]<-rawdata$label[rawdata$cdr3b=="CASSWWTGGSAPIYF"]
  df$label[df[,num]=="CASSWWTGGSSPIRF"]<-rawdata$label[rawdata$cdr3b=="CASSWWTGGSSPIRF"]
  df$label[df[,num]=="CASSWWTGGAAPIRF"]<-rawdata$label[rawdata$cdr3b=="CASSWWTGGAAPIRF"]
  return(df)
}

stopifnot(sum(pmtnetv1$CDR3!=rawdata$cdr3b)==0)
stopifnot(sum(teim$cdr3!=rawdata$cdr3b)==0)
stopifnot(sum(imrex$cdr3!=rawdata$cdr3b)==0)
stopifnot(sum(ergo2ae.mac$TRB!=rawdata$cdr3b)==0)

teim$label<-rawdata$label
imrex$label<-rawdata$label
pmtnetv1$label<-rawdata$label
ergo2ae.mac$label<-rawdata$label
ergo2ae.vdj$label<-rawdata$label
ergo2ls.mac$label<-rawdata$label
ergo2ls.vdj$label<-rawdata$label

teim$ID<-rawdata$ID
imrex$ID<-rawdata$ID
pmtnetv1$ID<-rawdata$ID
ergo2ae.mac$ID<-rawdata$ID
ergo2ae.vdj$ID<-rawdata$ID
ergo2ls.mac$ID<-rawdata$ID
ergo2ls.vdj$ID<-rawdata$ID

stopifnot(sum(atm$V2!=rawdata$cdr3b)==0)

atm<-atm[order(atm$V2,atm$V1),]
rawdata<-rawdata[order(rawdata$cdr3b,rawdata$peptide),]
atm$label<-rawdata$label
atm$ID<-rawdata$ID

rawdata_p<-rawdata[rawdata$peptide=='SLLQHLIGL',]
rawdata_p<-rawdata_p[order(rawdata_p$cdr3a,rawdata_p$cdr3b),]
nettcr.prame<-nettcr.prame[order(nettcr.prame$cdr3a,nettcr.prame$cdr3b),]
epact.prame<-epact.prame[order(epact.prame$CDR3.alpha.aa,epact.prame$CDR3.beta.aa),]
tulip.prame<-tulip.prame[order(tulip.prame$cdr3a,tulip.prame$cdr3b),]

sum(rawdata_p$cdr3a!=nettcr.prame$cdr3a)
sum(rawdata_p$cdr3b!=nettcr.prame$cdr3b)
sum(rawdata_p$cdr3a!=epact.prame$CDR3.alpha.aa)
sum(rawdata_p$cdr3b!=epact.prame$CDR3.beta.aa)
sum(rawdata_p$cdr3a!=tulip.prame$cdr3a)
sum(rawdata_p$cdr3b!=tulip.prame$cdr3b)

epact.prame$label<-rawdata_p$label
epact.prame$ID<-rawdata_p$ID

tulip.prame$label<-rawdata_p$label
tulip.prame$ID<-rawdata_p$ID

set_kd_t2<-function(data){
  data['Kd']= NA
  data[data['ID']=='a28b50','Kd']=391
  data[data['ID']=='a28b60','Kd']=261
  data[data['ID']=='a28b74','Kd']=182
  data[data['ID']=='a28b75','Kd']=214
  data[data['ID']=='a28b57','Kd']=83
  data[data['ID']=='a28b58','Kd']=79
  data[data['ID']=='a79b46','Kd']=31.8
  data[data['ID']=='a109b46','Kd']=170
  data[data['ID']=='a79b63','Kd']=79
  data[data['ID']=='a79b64','Kd']=138
  data[data['ID']=='a79b66','Kd']=89
  data[data['ID']=='a79b67','Kd']=47
  data[data['ID']=='a79b69','Kd']=52
  data[data['ID']=='a79b71','Kd']=87
  data[data['ID']=='a79b58','Kd']=23.1
  data[data['ID']=='a79b73','Kd']=132
  data[data['ID']=='a79b74','Kd']=53.3
  data[data['ID']=='a79b75','Kd']=57.7
  data[data['ID']=='a79b76','Kd']=11.8
  data[data['ID']=='a79b77','Kd']=77.9
  
  data['t2']= NA
  data[data['ID']=='a28b50','t2']=1.8
  data[data['ID']=='a28b60','t2']=2.8
  data[data['ID']=='a28b74','t2']=3.7
  data[data['ID']=='a28b75','t2']=5.1
  data[data['ID']=='a28b57','t2']=8.3
  data[data['ID']=='a28b58','t2']=8.9
  data[data['ID']=='a79b46','t2']=29.2
  data[data['ID']=='a109b46','t2']=7.31
  data[data['ID']=='a79b63','t2']=10.8
  data[data['ID']=='a79b64','t2']=6.38
  data[data['ID']=='a79b66','t2']=9.16
  data[data['ID']=='a79b67','t2']=12.69
  data[data['ID']=='a79b69','t2']=20.41
  data[data['ID']=='a79b71','t2']=14.89
  data[data['ID']=='a79b58','t2']=28.7
  data[data['ID']=='a79b73','t2']=4.6
  data[data['ID']=='a79b74','t2']=12.5
  data[data['ID']=='a79b75','t2']=16.9
  data[data['ID']=='a79b76','t2']=58.3
  data[data['ID']=='a79b77','t2']=8.6
  
  return(data)
}

teim<-set_kd_t2(teim)
imrex<-set_kd_t2(imrex)
pmtnetv1<-set_kd_t2(pmtnetv1)
ergo2ae.mac<-set_kd_t2(ergo2ae.mac)
ergo2ae.vdj<-set_kd_t2(ergo2ae.vdj)
ergo2ls.mac<-set_kd_t2(ergo2ls.mac)
ergo2ls.vdj<-set_kd_t2(ergo2ls.vdj)
atm<-set_kd_t2(atm)

nettcr.prame<-set_kd_t2(nettcr.prame)
epact.prame<-set_kd_t2(epact.prame)
tulip.prame<-set_kd_t2(tulip.prame)

# make sure that results of all tools have same order
pmtnetv2<-pmtnetv2[order(pmtnetv2$cdr3b,pmtnetv2$peptide),]
pmtnetv1<-pmtnetv1[order(pmtnetv1$CDR3,pmtnetv1$Antigen),]
teim<-teim[order(teim$cdr3,teim$epitope),]
atm<-atm[order(atm$V2,atm$V1),]
imrex<-imrex[order(imrex$cdr3,imrex$antigen.epitope),]
ergo2ae.mac<-ergo2ae.mac[order(ergo2ae.mac$TRB,ergo2ae.mac$Peptide),]
ergo2ae.vdj<-ergo2ae.vdj[order(ergo2ae.vdj$TRB,ergo2ae.vdj$Peptide),]
ergo2ls.mac<-ergo2ls.mac[order(ergo2ls.mac$TRB,ergo2ls.mac$Peptide),]
ergo2ls.vdj<-ergo2ls.vdj[order(ergo2ls.vdj$TRB,ergo2ls.vdj$Peptide),]

stopifnot(sum(pmtnetv2$cdr3b!=pmtnetv1$CDR3)==0)
stopifnot(sum(pmtnetv2$peptide!=pmtnetv1$Antigen)==0)
stopifnot(sum(pmtnetv2$cdr3b!=teim$cdr3)==0)
stopifnot(sum(pmtnetv2$peptide!=teim$epitope)==0)
stopifnot(sum(pmtnetv2$cdr3b!=atm$V2)==0)
stopifnot(sum(pmtnetv2$peptide!=atm$V1)==0)
stopifnot(sum(pmtnetv2$cdr3b!=imrex$cdr3)==0)
stopifnot(sum(pmtnetv2$peptide!=imrex$antigen.epitope)==0)
stopifnot(sum(pmtnetv2$cdr3b!=ergo2ae.mac$TRB)==0)
stopifnot(sum(pmtnetv2$peptide!=ergo2ae.mac$Peptide)==0)

plot.data<-data.frame(pmtnetv2=c(ifelse(pmtnetv2$avg_percentile_rank[pmtnetv2$peptide=='SLLQHLIGL' & pmtnetv2$label==1]<1.5e-2,1,0),sum(pmtnetv2$avg_percentile_rank[pmtnetv2$peptide=='SLLQHLIGL' & pmtnetv2$label==0]>1.5e-2)/750),
                      pmtnetv1=c(ifelse(pmtnetv1$Rank[pmtnetv1$Antigen=='SLLQHLIGL' & pmtnetv1$label==1]<1.5e-2,1,0),sum(pmtnetv1$Rank[pmtnetv1$Antigen=='SLLQHLIGL' & pmtnetv1$label==0]>1.5e-2)/750),
                      teim=c(ifelse(teim$binding[teim$epitope=='SLLQHLIGL' & teim$label==1]>quantile(teim$binding,1-1.5e-2),1,0),sum(teim$binding[teim$epitope=='SLLQHLIGL' & teim$label==0]<quantile(teim$binding[teim$epitope=='SLLQHLIGL'],1-1.5e-2))/750),
                      atm=c(ifelse(atm$V5[atm$V1=='SLLQHLIGL' & atm$label==1]>quantile(atm$V5,1-1.5e-2),1,0),sum(atm$V5[atm$V1=='SLLQHLIGL' & atm$label==0]<quantile(atm$V5[atm$V1=='SLLQHLIGL'],1-1.5e-2))/750),
                      imrex=c(ifelse(imrex$prediction_score[imrex$antigen.epitope=='SLLQHLIGL' & imrex$label==1]>quantile(imrex$prediction_score,1-1.5e-2),1,0),sum(imrex$prediction_score[imrex$antigen.epitope=='SLLQHLIGL' & imrex$label==0]<quantile(imrex$prediction_score[imrex$antigen.epitope=='SLLQHLIGL'],1-1.5e-2))/750),
                      ergo2aemac=c(ifelse(ergo2ae.mac$Score[ergo2ae.mac$Peptide=='SLLQHLIGL' & ergo2ae.mac$label==1]>quantile(ergo2ae.mac$Score,1-1.5e-2),1,0),sum(ergo2ae.mac$Score[ergo2ae.mac$Peptide=='SLLQHLIGL' & ergo2ae.mac$label==0]<quantile(ergo2ae.mac$Score[ergo2ae.mac$Peptide=='SLLQHLIGL'],1-1.5e-2))/750),
                      ergo2aevdj=c(ifelse(ergo2ae.vdj$Score[ergo2ae.vdj$Peptide=='SLLQHLIGL' & ergo2ae.vdj$label==1]>quantile(ergo2ae.vdj$Score,1-1.5e-2),1,0),sum(ergo2ae.vdj$Score[ergo2ae.vdj$Peptide=='SLLQHLIGL' & ergo2ae.vdj$label==0]<quantile(ergo2ae.vdj$Score[ergo2ae.vdj$Peptide=='SLLQHLIGL'],1-1.5e-2))/750),
                      ergo2lsmac=c(ifelse(ergo2ls.mac$Score[ergo2ls.mac$Peptide=='SLLQHLIGL' & ergo2ls.mac$label==1]>quantile(ergo2ls.mac$Score,1-1.5e-2),1,0),sum(ergo2ls.mac$Score[ergo2ls.mac$Peptide=='SLLQHLIGL' & ergo2ls.mac$label==0]<quantile(ergo2ls.mac$Score[ergo2ls.mac$Peptide=='SLLQHLIGL'],1-1.5e-2))/750),
                      ergo2lsvdj=c(ifelse(ergo2ls.vdj$Score[ergo2ls.vdj$Peptide=='SLLQHLIGL' & ergo2ls.vdj$label==1]>quantile(ergo2ls.vdj$Score,1-1.5e-2),1,0),sum(ergo2ls.vdj$Score[ergo2ls.vdj$Peptide=='SLLQHLIGL' & ergo2ls.vdj$label==0]<quantile(ergo2ls.vdj$Score[ergo2ls.vdj$Peptide=='SLLQHLIGL'],1-1.5e-2))/750),
                      epact=c(rep(0,15),sum(epact.prame$Pred[epact.prame$label==0]<quantile(epact.prame$Pred,1-1.5e-2))/750),
                      nettcr=c(rep(0,15),sum(nettcr.prame$prediction[nettcr.prame$label==0]<quantile(nettcr.prame$prediction,1-1.5e-2))/750),
                      tulip=c(rep(0,15),sum(tulip.prame$rank[tulip.prame$label==0]>quantile(tulip.prame$rank,1.5e-2))/750)
)

plot.data<-plot.data[1:15,]
plot.data$tcr<-paste('tcr_',15:1,sep='_')

plot.data.reshape<-reshape(plot.data,idvar = 'tcr',varying = list(names(plot.data)[1:12]),timevar = 'method',direction = 'long',
                           times = c("pmtnetv2","pmtnetv1","teim","atm","imrex","ergo2aemac","ergo2aevdj","ergo2lsmac","ergo2lsvdj","epact","nettcr","tulip"),v.names = 'value')

plot.data.reshape$method<-factor(plot.data.reshape$method,levels = c("pmtnetv2","pmtnetv1",'nettcr',"teim","epact","atm","tulip","imrex","ergo2aevdj","ergo2aemac","ergo2lsvdj","ergo2lsmac"))

myPalette<-colorRampPalette(c("grey90","#BEE4A0", "#FDBE6F", "#9E0142"))
ggplot(plot.data.reshape,aes(x=method,y=tcr))+
  geom_point(aes(color=value),size=5)+ 
  xlab("")+theme_classic()+
  theme(axis.text.x = element_text(angle=20,vjust = 0.6))+
  scale_colour_gradientn(colours = myPalette(4), limits=c(0, 1))

plot.data<-plot.data[16,]
plot.data$tcr<-'neg'
plot.data.reshape<-reshape(plot.data,idvar = 'tcr',varying = list(names(plot.data)[1:12]),timevar = 'method',direction = 'long',
                           times = c("pmtnetv2","pmtnetv1","teim","atm","imrex","ergo2aemac","ergo2aevdj","ergo2lsmac","ergo2lsvdj","epact","nettcr","tulip"),v.names = 'value')
plot.data.reshape$method<-factor(plot.data.reshape$method,levels = c("pmtnetv2","pmtnetv1",'nettcr',"teim","epact","atm","tulip","imrex","ergo2aevdj","ergo2aemac","ergo2lsvdj","ergo2lsmac"))

ggplot(plot.data.reshape,aes(x=method,y=tcr,color=value))+
  geom_point(size=5)+
  xlab("")+theme_classic()+
  theme(axis.text.x = element_text(angle=20,vjust = 0.6))+
  scale_colour_gradientn(colours = myPalette(4), limits=c(0, 1))

# Fig. 4 KRAS
predata<-read.table('alex_kras_fewshot_G12VA11training_newPos_AlexValidationOnly_pre_seed100.csv',header = T,sep = ',')
pmtnetv1<-read.table('pMTnetv1_prediction_g12v_newPos.csv',header = T,sep=',')
teim<-read.table('teim_g12v_newPos_level_binding.csv',header = T,sep=',')
atm<-read.table('pred_original_input_g12v_negPos_atmtcr.csv',header = F,sep='\t')
imrex<-read.table('imrex_output_g12v_negPos.csv',header = T,sep=',')
ergo2ae.mac<-read.table('ergo2_g12v_newPos_ae_mc.csv',header = T,sep=',')
ergo2ae.vdj<-read.table('ergo2_g12v_newPos_ae_vdj.csv',header = T,sep=',')
ergo2ls.mac<-read.table('ergo2_g12v_newPos_lstm_mc.csv',header = T,sep=',')
ergo2ls.vdj<-read.table('ergo2_g12v_newPos_lstm_vdj.csv',header = T,sep=',')
epact.kras<-read.table('preds_epact_origin.csv',header = T,sep=',')
epact.kras<-epact.kras[epact.kras$Epitope.peptide!='VVGAVGVGK' | epact.kras$MHC!='A*11:01',]
nettcr.kras<-read.table('alex_nettcr_result.csv',header = T,sep=',')
nettcr.kras<-nettcr.kras[nettcr.kras$peptide!='VVGAVGVGK' | nettcr.kras$mhca!='A*11:01',]
lines <- readLines("ranking-tulip-origin.jsonl")
tulip.kras <- lapply(lines, fromJSON)
tulip.kras <- rbind(data.frame(peptide=rep(tulip.kras[[1]]$peptide),cdr3a=tulip.kras[[1]]$CDR3a,cdr3b=tulip.kras[[1]]$CDR3b,rank=tulip.kras[[1]]$rank),
                    data.frame(peptide=rep(tulip.kras[[2]]$peptide),cdr3a=tulip.kras[[2]]$CDR3a,cdr3b=tulip.kras[[2]]$CDR3b,rank=tulip.kras[[2]]$rank))
tulip.kras<-tulip.kras[(tulip.kras$cdr3a %in% epact.kras$CDR3.alpha.aa) & (tulip.kras$cdr3b %in% epact.kras$CDR3.beta.aa),]

sum(predata$cdr3b==pmtnetv1$CDR3)+sum(predata$peptide==pmtnetv1$Antigen)
sum(predata$cdr3b==teim$cdr3)+sum(predata$peptide==teim$epitope)
sum(predata$cdr3b==imrex$cdr3)+sum(predata$peptide==imrex$antigen.epitope)
sum(predata$cdr3b==ergo2ae.mac$TRB)+sum(predata$peptide==ergo2ae.mac$Peptide)+sum(predata$cdr3a==ergo2ae.mac$TRA)
sum(predata$cdr3b==ergo2ae.vdj$TRB)+sum(predata$peptide==ergo2ae.vdj$Peptide)+sum(predata$cdr3a==ergo2ae.vdj$TRA)
sum(predata$cdr3b==ergo2ls.mac$TRB)+sum(predata$peptide==ergo2ls.mac$Peptide)+sum(predata$cdr3a==ergo2ls.mac$TRA)
sum(predata$cdr3b==ergo2ls.vdj$TRB)+sum(predata$peptide==ergo2ls.vdj$Peptide)+sum(predata$cdr3a==ergo2ls.vdj$TRA)

pmtnetv1$label<-predata$label
teim$label<-predata$label
imrex$label<-predata$label
ergo2ae.mac$label<-predata$label
ergo2ae.vdj$label<-predata$label
ergo2ls.mac$label<-predata$label
ergo2ls.vdj$label<-predata$label

atm<-atm[order(atm$V1,atm$V2),]
predata<-predata[order(predata$peptide,predata$cdr3b),]
sum(predata$cdr3b==atm$V2)+sum(predata$peptide==atm$V1)
atm$label<-predata$label

predata<-predata[order(predata$peptide,predata$cdr3a,predata$cdr3b),]
epact.kras<-epact.kras[order(epact.kras$Epitope.peptide,epact.kras$CDR3.alpha.aa,epact.kras$CDR3.beta.aa),]
sum(predata$cdr3a==epact.kras$CDR3.alpha.aa)+sum(predata$cdr3b==epact.kras$CDR3.beta.aa)+sum(predata$peptide==epact.kras$Epitope.peptide)
epact.kras$label<-predata$label

nettcr.kras<-nettcr.kras[order(nettcr.kras$peptide,nettcr.kras$cdr3a,nettcr.kras$cdr3b),]
sum(predata$cdr3a==nettcr.kras$cdr3a)+sum(predata$cdr3b==nettcr.kras$cdr3b)+sum(predata$peptide==nettcr.kras$peptide)
nettcr.kras$label<-predata$label

tulip.kras<-tulip.kras[order(tulip.kras$peptide,tulip.kras$cdr3a,tulip.kras$cdr3b),]
sum(predata$cdr3a==tulip.kras$cdr3a)+sum(predata$cdr3b==tulip.kras$cdr3b)+sum(predata$peptide==tulip.kras$peptide)
tulip.kras$label<-predata$label

cutoff<-0.015

ifelse(pmtnetv1$Rank[pmtnetv1$label==1]<cutoff,1,0)
ifelse(teim$binding[teim$label==1]>quantile(teim$binding,1-cutoff),1,0)
ifelse(atm$V5[atm$label==1]>quantile(atm$V5,1-cutoff),1,0)
ifelse(imrex$prediction_score[imrex$label==1]>quantile(imrex$prediction_score,1-cutoff),1,0)
ifelse(ergo2ae.mac$Score[ergo2ae.mac$label==1]>quantile(ergo2ae.mac$Score,1-cutoff),1,0)
ifelse(ergo2ae.vdj$Score[ergo2ae.vdj$label==1]>quantile(ergo2ae.vdj$Score,1-cutoff),1,0)
ifelse(ergo2ls.mac$Score[ergo2ls.mac$label==1]>quantile(ergo2ls.mac$Score,1-cutoff),1,0)
ifelse(ergo2ls.vdj$Score[ergo2ls.vdj$label==1]>quantile(ergo2ls.vdj$Score,1-cutoff),1,0)
ifelse(epact.kras$Pred[epact.kras$label==1]>quantile(epact.kras$Pred,1-cutoff),1,0)
ifelse(nettcr.kras$prediction[nettcr.kras$label==1]>quantile(nettcr.kras$prediction,1-cutoff),1,0)
ifelse(tulip.kras$rank[tulip.kras$label==1 & tulip.kras$peptide=="VVGACGVGK"]<quantile(tulip.kras$rank[tulip.kras$peptide=="VVGACGVGK"],cutoff),1,0)
ifelse(tulip.kras$rank[tulip.kras$label==1 & tulip.kras$peptide=="VVGAVGVGK"]<quantile(tulip.kras$rank[tulip.kras$peptide=="VVGAVGVGK"],cutoff),1,0)

plot.data<-data.frame(pmtnetv2=c(ifelse(predata$avg_percentile_rank[predata$label==1]<cutoff,1,0),sum(predata$avg_percentile_rank[predata$label==0]>cutoff)/sum(predata$label==0)),
                      pmtnetv1=c(ifelse(pmtnetv1$Rank[pmtnetv1$label==1]<cutoff,1,0),sum(pmtnetv1$Rank[pmtnetv1$label==0]>cutoff)/sum(pmtnetv1$label==0)),
                      teim=c(ifelse(teim$binding[teim$label==1]>quantile(teim$binding,1-cutoff),1,0),sum(teim$binding[teim$label==0]<quantile(teim$binding,1-cutoff))/sum(teim$label==0)),
                      atm=c(ifelse(atm$V5[atm$label==1]>quantile(atm$V5,1-cutoff),1,0),sum(atm$V5[atm$label==0]<quantile(atm$V5,1-cutoff))/sum(atm$label==0)),
                      imrex=c(ifelse(imrex$prediction_score[imrex$label==1]>quantile(imrex$prediction_score,1-cutoff),1,0),sum(imrex$prediction_score[imrex$label==0]<quantile(imrex$prediction_score,1-cutoff))/sum(imrex$label==0)),
                      ergo2aemac=c(ifelse(ergo2ae.mac$Score[ergo2ae.mac$label==1]>quantile(ergo2ae.mac$Score,1-cutoff),1,0),sum(ergo2ae.mac$Score[ergo2ae.mac$label==0]<quantile(ergo2ae.mac$Score,1-cutoff))/sum(ergo2ae.mac$label==0)),
                      ergo2aevdj=c(ifelse(ergo2ae.vdj$Score[ergo2ae.vdj$label==1]>quantile(ergo2ae.vdj$Score,1-cutoff),1,0),sum(ergo2ae.vdj$Score[ergo2ae.vdj$label==0]<quantile(ergo2ae.vdj$Score,1-cutoff))/sum(ergo2ae.vdj$label==0)),
                      ergo2lsmac=c(ifelse(ergo2ls.mac$Score[ergo2ls.mac$label==1]>quantile(ergo2ls.mac$Score,1-cutoff),1,0),sum(ergo2ls.mac$Score[ergo2ls.mac$label==0]<quantile(ergo2ls.mac$Score,1-cutoff))/sum(ergo2ls.mac$label==0)),
                      ergo2lsvdj=c(ifelse(ergo2ls.vdj$Score[ergo2ls.vdj$label==1]>quantile(ergo2ls.vdj$Score,1-cutoff),1,0),sum(ergo2ls.vdj$Score[ergo2ls.vdj$label==0]<quantile(ergo2ls.vdj$Score,1-cutoff))/sum(ergo2ls.vdj$label==0)),
                      epact=c(ifelse(epact.kras$Pred[epact.kras$label==1]>quantile(epact.kras$Pred,1-cutoff),1,0),sum(epact.kras$Pred[epact.kras$label==0]<quantile(epact.kras$Pred,1-cutoff))/sum(epact.kras$label==0)),
                      nettcr=c(ifelse(nettcr.kras$prediction[nettcr.kras$label==1]>quantile(nettcr.kras$prediction,1-cutoff),1,0),sum(nettcr.kras$prediction[nettcr.kras$label==0]<quantile(nettcr.kras$prediction,1-cutoff))/sum(nettcr.kras$label==0)),
                      tulip=c(ifelse(tulip.kras$rank[tulip.kras$label==1 & tulip.kras$peptide=="VVGACGVGK"]<quantile(tulip.kras$rank[tulip.kras$peptide=="VVGACGVGK"],cutoff),1,0),
                              ifelse(tulip.kras$rank[tulip.kras$label==1 & tulip.kras$peptide=="VVGAVGVGK"]<quantile(tulip.kras$rank[tulip.kras$peptide=="VVGAVGVGK"],cutoff),1,0),
                              sum(c(tulip.kras$rank[tulip.kras$label==0 & tulip.kras$peptide=="VVGACGVGK"]>quantile(tulip.kras$rank[tulip.kras$peptide=="VVGACGVGK"],cutoff),
                                    tulip.kras$rank[tulip.kras$label==0 & tulip.kras$peptide=="VVGAVGVGK"]>quantile(tulip.kras$rank[tulip.kras$peptide=="VVGAVGVGK"],cutoff)))/sum(tulip.kras$label==0))
)
myPalette<-colorRampPalette(c("grey90","#BEE4A0", "#FDBE6F", "#9E0142"))
plot.data1<-plot.data[1:7,]
plot.data1$tcr<-paste('tcr_',7:1,sep='_')
plot.data1.reshape<-reshape(plot.data1,idvar = 'tcr',varying = list(names(plot.data1)[1:12]),timevar = 'method',direction = 'long',
                            times = c("pmtnetv2","pmtnetv1","teim","atm","imrex","ergo2aemac","ergo2aevdj","ergo2lsmac","ergo2lsvdj","epact","nettcr","tulip"),v.names = 'value')
plot.data1.reshape$method<-factor(plot.data1.reshape$method,levels = c("pmtnetv2","pmtnetv1",'nettcr',"teim","epact","atm","tulip","imrex","ergo2aevdj","ergo2aemac","ergo2lsvdj","ergo2lsmac"))

ggplot(plot.data1.reshape,aes(x=method,y=tcr))+
  geom_point(aes(color=value),size=5)+ 
  xlab("")+theme_classic()+
  theme(axis.text.x = element_text(angle=20,vjust = 0.6))+
  scale_colour_gradientn(colours = myPalette(4), limits=c(0, 1))

plot.data2<-plot.data[8,]
plot.data2$tcr<-'neg'
plot.data2.reshape<-reshape(plot.data2,idvar = 'tcr',varying = list(names(plot.data2)[1:12]),timevar = 'method',direction = 'long',
                            times = c("pmtnetv2","pmtnetv1","teim","atm","imrex","ergo2aemac","ergo2aevdj","ergo2lsmac","ergo2lsvdj","epact","nettcr","tulip"),v.names = 'value')
plot.data2.reshape$method<-factor(plot.data2.reshape$method,levels = c("pmtnetv2","pmtnetv1",'nettcr',"teim","epact","atm","tulip","imrex","ergo2aevdj","ergo2aemac","ergo2lsvdj","ergo2lsmac"))

ggplot(plot.data2.reshape,aes(x=method,y=tcr,color=value))+
  geom_point(size=5)+
  xlab("")+theme_classic()+
  theme(axis.text.x = element_text(angle=20,vjust = 0.6))+
  scale_colour_gradientn(colours = myPalette(4), limits=c(0, 1))

# Sup. File 2 Fig.4 unseen epitope
# prepare input data
training<-read.table("training_individual_1.txt",header = T,sep = "\t")
training_peptide<-unique(training$peptide)

validation<-read.table("feb20_pseudoVgene_upHc2Mouse_val5neg_train50neg_maxEpoch70/validation.txt",header = T,sep = "\t")
validation_neg<-validation[(validation$label==0) & (validation$TCR_SPECIES=="human"),c("vaseq","cdr3a","vbseq","cdr3b")]
validation_neg<-validation_neg[!duplicated(validation_neg),]

immrep25.iedb<-read.table("iedb_positives.csv",header = T,sep=',')
immrep25.iedb$dataset<-'immrep25iedb'
immrep25.vdjdb<-read.table("vdjdb_positives.csv",header = T,sep=',',comment.char = "", na.strings = c("", "*"), stringsAsFactors = FALSE)
immrep25.vdjdb$dataset<-'immrep25vdjdb'
immrep25<-rbind(immrep25.iedb,immrep25.vdjdb)

immrep25.unique.epitope<-immrep25[apply(stringdist::stringdistmatrix(immrep25$Peptide,training_peptide),1,min)>2,]
immrep25.unique.epitope<-immrep25.unique.epitope[immrep25.unique.epitope$just_10X=='False',]
immrep25.unique.epitope<-immrep25.unique.epitope[!grepl("HLA-E*",immrep25.unique.epitope$HLA,fixed = T),]
immrep25.unique.epitope$HLA[immrep25.unique.epitope$HLA=="HLA-A*02:01:48"]<-"HLA-A*02:01"
immrep25.unique.epitope$HLA[immrep25.unique.epitope$HLA=="HLA-B*44:05:01"]<-"HLA-B*44:05"
immrep25.unique.epitope$HLA2<-substr(immrep25.unique.epitope$HLA, 5, nchar(immrep25.unique.epitope$HLA))

input.pmtnetv2<-immrep25.unique.epitope[,c("Peptide","HLA2","Va","CDR3a_extended","Vb","CDR3b_extended")]
input.pmtnetv2$mhcb<-"human_microglobulin"
input.pmtnetv2$pMHC_SPECIES<-'human'
input.pmtnetv2$TCR_SPECIES<-'human'
input.pmtnetv2$label<-1

human_a<-Biostrings::readAAStringSet(filepath = "human_trav_pro.txt.completed.txt")
human_name_a <- names(human_a)
human_seq_a <- paste(human_a)
vgene_human_a<-data.frame(a_allele=sapply(strsplit(human_name_a,split = "|",fixed = T),function(x){x[2]}), human_seq_a)
rownames(vgene_human_a)<-vgene_human_a$a_allele

human_b<-Biostrings::readAAStringSet(filepath = "human_trbv_pro.txt.completed.txt")
human_name_b <- names(human_b)
human_seq_b <- paste(human_b)
vgene_human_b<-data.frame(b_allele=sapply(strsplit(human_name_b,split = "|",fixed = T),function(x){x[2]}), human_seq_b)
rownames(vgene_human_b)<-vgene_human_b$b_allele

stopifnot(length(setdiff(input.pmtnetv2$Va,vgene_human_a$a_allele))==0)
stopifnot(length(setdiff(input.pmtnetv2$Vb,vgene_human_b$b_allele))==0)

input.pmtnetv2$vaseq<-vgene_human_a[input.pmtnetv2$Va,"human_seq_a"]
input.pmtnetv2$vbseq<-vgene_human_b[input.pmtnetv2$Vb,"human_seq_b"]

stopifnot(sum(is.na(input.pmtnetv2))==0)

colnames(input.pmtnetv2)<-c("peptide","mhca","va","cdr3a","vb","cdr3b","mhcb","pMHC_SPECIES","TCR_SPECIES","label","vaseq","vbseq")
input.pmtnetv2<-input.pmtnetv2[,c("peptide","mhca","mhcb","pMHC_SPECIES","TCR_SPECIES","vaseq","cdr3a","vbseq","cdr3b","label")]

h_neg_back_a<-read.table("human_alpha.txt",header = T,sep="\t",quote = "") #dim=822103      3
h_neg_back_b<-read.table("human_beta.txt",header = T,sep="\t",quote = "")  #dim=5233069       3

input.pmtnetv2.neg<-cbind(input.pmtnetv2[rep(1:nrow(input.pmtnetv2),5),c("peptide","mhca","mhcb","pMHC_SPECIES","TCR_SPECIES")],
                          h_neg_back_a[sample(1:nrow(h_neg_back_a),5*nrow(input.pmtnetv2)),c("vaseq","cdr3a")],
                          h_neg_back_b[sample(1:nrow(h_neg_back_b),5*nrow(input.pmtnetv2)),c("vbseq","cdr3b")],
                          data.frame(label=rep(0,5*nrow(input.pmtnetv2))))

input.pmtnetv2<-rbind(input.pmtnetv2,input.pmtnetv2.neg)

write.table(input.pmtnetv2,"immrep25_pmtnetv2_input_new_random_neg.tsv",quote = F,sep = "\t",row.names = F,col.names = T)

#pmtnetv1
input.pmtnetv1.data<-input.pmtnetv2[,c("cdr3b","peptide","mhca","label")]
input.pmtnetv1.data<-input.pmtnetv1.data[!duplicated(input.pmtnetv1.data),]
input.pmtnetv1<-input.pmtnetv1.data[,c("cdr3b","peptide","mhca")]
write.table(input.pmtnetv1,"immrep25_pmtnetv1_input_new_random_neg.csv",quote = F,sep = ",",row.names = F,col.names = T)

#TEIM
input.teim.data<-input.pmtnetv2[,c("cdr3b","peptide","label")]
colnames(input.teim.data)<-c("cdr3","epitope","label")
input.teim.data<-input.teim.data[!duplicated(input.teim.data),]
input.teim.data<-input.teim.data[nchar(input.teim.data$cdr3)<=20,]
input.teim<-input.teim.data[,c("cdr3","epitope")]
write.table(input.teim,"immrep25_teim_input_new_random_neg.csv",quote = F,sep = ",",row.names = F,col.names = T)

# ATM-TCR
atmtcr_data<-data.frame(input.pmtnetv2$peptide,input.pmtnetv2$cdr3b,input.pmtnetv2$label) #no header
atmtcr_data<-atmtcr_data[!duplicated(atmtcr_data),]
atmtcr_data<-atmtcr_data[nchar(atmtcr_data$input.pmtnetv2.cdr3b)<=20,]
atmtcr<-cbind(atmtcr_data[,1:2],rep(0,nrow(atmtcr_data)))
write.table(atmtcr,"immrep25_atmtcr_input_new_random_neg.csv",quote = F,sep = ",",row.names = F,col.names = F)

#ImRex
imrex<-data.frame(cdr3=input.pmtnetv2$cdr3b,antigen.epitope=input.pmtnetv2$peptide) #sep=；
imrex<-imrex[!duplicated(imrex),]
write.table(imrex,"immrep25_imrex_input_new_random_neg.csv",col.names = T, row.names = F,sep = ';',quote = F)

# ERGOII
ergo2<-data.frame(TRA=input.pmtnetv2$cdr3a,TRB=input.pmtnetv2$cdr3b,
                  TRAV=rep(""),TRAJ=rep(""),TRBV=rep(""),TRBJ=rep(""),`T-Cell-Type`=rep(""),
                  Peptide=input.pmtnetv2$peptide,
                  MHC=paste('HLA-',unlist(lapply(strsplit(input.pmtnetv2$mhca,split = ":"),'[[',1)),sep = "")) # write  'TRA', 'TRB', 'TRAV', 'TRAJ', 'TRBV', 'TRBJ', 'T_Cell_Type', 'Peptide' and 'MHC', then change T_Cell_Type to T-Cell-Type
ergo2<-ergo2[!duplicated(ergo2),]
write.table(ergo2,"immrep25_ergo_input_new_random_neg.csv",col.names = T, row.names = F,sep = ',',quote = F)

#NetTCR
immrep25$vaseq<-vgene_human_a[immrep25$Va,"human_seq_a"]
immrep25$vbseq<-vgene_human_b[immrep25$Vb,"human_seq_b"]
cdr12a<-immrep25[,c("vaseq","CDR1a","CDR2a")]
cdr12a<-rbind(cdr12a,data.frame(vaseq=c("TQSVTQLDGHITVSEEAPLELKCNYSYSGVPSLFWYVQYSSQSLQLLLKDLTEATQVKGIRGFEAEFKKSETSFYLRKPSTHVSDAAEYFCAVG","AQSVTQLGSHVSVSEGALVLLRGATHYCCPPILFWYVQYPNQGLQLLLKYTSAATLVKGINGFEAEFKKSETSFHLTKPAAHMSDAAEYFCAVS"),
                                CDR1a=c("YSGVPS","YCCPPI"),
                                CDR2a=c("DLTEATQV","YTSAATLV")))
cdr12a<-cdr12a[!duplicated(cdr12a),]
cdr12a<-na.omit(cdr12a)
rownames(cdr12a)<-cdr12a$vaseq
cdr12b<-immrep25[,c("vbseq","CDR1b","CDR2b")]
cdr12b<-rbind(cdr12b,data.frame(vbseq=c("DTKVTQRPRLLVKASEQKAKMDCVPIKAHSYVYWYRKKLEEELKFLVYFQNEELIQKAEIINERFLAQCSKNSSCTLEIQSTESGDTALYFCASSK","DTGITQTPKYLVTAMGSKRTMKREHLGHDSMYWYRQKAKKSLEFMFYYNCKEFIENKTVPNHFTPECPDSSRLYLHVVALQQEDSAAYLCTSSQ",
                                        "GAVVSQHPSRVICKSGTSVKIECRSLDFQATTMFWYRQFPKKSLMLMATSNEGSKATYEQGVEKDKFLINHASLTLSTLTVTSAHPEDSSFYICSAS","GAGVSQSPRYEVTQRGQDVAPRCDPISGQVTLYWYRQTLGQGQEFLTSFQDETQQDKSGLLSDQFSTERSEDLSPPEDPAHRARATRLCISVPEA",
                                        "SQTIHQWPATLVQPVGSPLSLECTVEGTSNPNLYWYRQAAGRGLQLLFYSVGIGQISSEVPQNLSASRPQDRQFILSSKKLLLSDSGFYLCAWG","SAVISQKPSRDICQRGTSLTIQCQVDSQVTMIFWYRQQPGQSLTLIATANQGSEATYESGFVIDKFPISRPNLTFSTLTVSNMSPEDSSIYLCSAG",
                                        "GAGVSQSPRYKVAKRGQDVALRCDPISGHVSLFWYQQALGQGPEFLTYFQNEAQLDKSGLPSDRFFAERPEGSVSTLKIQRTQQEDSAVYLCASSR","DTGVSQNPRHKITKRGQNVTFRCDPISEHNRLYWYRQTLGQGPEFLTYFQNEAQLEKSRLLSDRFSAERPKGSLSTLEIQRTEQGDSAMYLCASTL",
                                        "NAGVTQTPKFRILKIGQSMTLQCTQDMNHEYMYWYRQDPGMGLKLIYYSVGAGITDKGEVPNGYNVSRSTTEDFPLRLELAAPSQTSVYFCASSR"),
                                CDR1b=c("KAHSY","LGHDS","DFQATT","SGQVT","GTSNPN","SQVTM","SGHVS","SEHNR","MNHEY"),
                                CDR2b=c("FQNEEL","YNCKEF","SNEGSKA","FQDETQ","SVGIG","ANQGSEA","FQNEAQ","FQNEAQ","SVGAGI")))
cdr12b<-cdr12b[!duplicated(cdr12b),]
cdr12b<-na.omit(cdr12b)
rownames(cdr12b)<-cdr12b$vbseq

unique(input.pmtnetv2$vaseq[!(input.pmtnetv2$vaseq %in% cdr12a$vaseq)])
unique(input.pmtnetv2$vbseq[!(input.pmtnetv2$vbseq %in% cdr12b$vbseq)])

nettcr<-data.frame(peptide=input.pmtnetv2$peptide,
                   A1=cdr12a[input.pmtnetv2$vaseq,"CDR1a"],A2=cdr12a[input.pmtnetv2$vaseq,"CDR2a"],A3=input.pmtnetv2$cdr3a,
                   B1=cdr12b[input.pmtnetv2$vbseq,"CDR1b"],B2=cdr12b[input.pmtnetv2$vbseq,"CDR2b"],B3=input.pmtnetv2$cdr3b,label=input.pmtnetv2$label)

rownames(nettcr)<-NULL
stopifnot(sum(is.na(nettcr))==0)
nettcr.lack<-nettcr[!((nchar(nettcr$A2)<=7) & (nchar(nettcr$B3)<=23)),]
nettcr<-nettcr[(nchar(nettcr$A2)<=7) & (nchar(nettcr$B3)<=23),]
write.table(nettcr,"immrep25_nettcr_input_new_random_neg.csv",col.names = T, row.names = F,sep = ',',quote = F)
write.table(nettcr.lack,"immrep25_nettcr_input_new_random_neg_lack.csv",col.names = T, row.names = F,sep = ',',quote = F)

#EPACT
epact<-data.frame(CDR3.alpha.aa=input.pmtnetv2$cdr3a,
                  CDR3.beta.aa=input.pmtnetv2$cdr3b,
                  Epitope.peptide=input.pmtnetv2$peptide,
                  MHC=input.pmtnetv2$mhca)
epact<-epact[!duplicated(epact),]
write.table(epact,"immrep25_epact_input_new_random_neg.csv",col.names = T, row.names = F,sep = ',',quote = F)

#TULIP
tulip<-data.frame(CDR3b=input.pmtnetv2$cdr3b,	CDR3a=input.pmtnetv2$cdr3a,	peptide=input.pmtnetv2$peptide,	MHC=paste0("HLA-",input.pmtnetv2$mhca))
tulip<-tulip[!duplicated(tulip),]
write.table(tulip,"immrep25_tulip_input_new_random_neg.csv",col.names = T, row.names = F,sep = ',',quote = F)

get_partial_auc01<-function(epitope,label,prediction){
  auc.par.all<-c()
  auc_min<-0.1*0.1/2
  auc_max<-0.1
  for (i in unique(epitope)) {
    auc.par.tmp<-pROC::auc(response=label[epitope==i],predictor=prediction[epitope==i],
                           direction="<",
                           partial.auc=c(1,0.9),partial.auc.focus="specificity", partial.auc.correct=F,quiet=T)
    auc.par.tmp<-((as.numeric(auc.par.tmp) - auc_min) / (auc_max - auc_min) + 1)*0.5
    auc.par.all<-c(auc.par.all,auc.par.tmp)
  }
  return(mean(auc.par.all))
}

pmtnetv2<-read.table("immrep25_pmtnetv2_output_all_mhc_new_random_neg.tsv",header = T,sep=",")
get_partial_auc01(epitope = pmtnetv2$peptide,label = pmtnetv2$label,prediction = pmtnetv2$avg_rank_asc)

pmtnetv1.out<-read.table('pMTnet_immrep25_new_random_neg.csv',header = T,sep=',')
input.pmtnetv2<-read.table("immrep25_pmtnetv2_input_new_random_neg.tsv",sep = '\t',header = T)
pmtnetv1.input<-input.pmtnetv2[,c("cdr3b","peptide","mhca","label")]
pmtnetv1.input<-pmtnetv1.input[!duplicated(pmtnetv1.input),]

stopifnot(identical(pmtnetv1.input$cdr3b,pmtnetv1.out$CDR3))
stopifnot(identical(pmtnetv1.input$peptide,pmtnetv1.out$Antigen))
stopifnot(identical(pmtnetv1.input$mhca,pmtnetv1.out$HLA))

get_partial_auc01(epitope = pmtnetv1.out$Antigen,label = pmtnetv1.input$label,prediction = 1-pmtnetv1.out$Rank)

#TEIM
input.pmtnetv2<-read.table("immrep25_pmtnetv2_input_new_random_neg.tsv",sep = '\t',header = T)
input.teim.data<-input.pmtnetv2[,c("cdr3b","peptide","label")]
colnames(input.teim.data)<-c("cdr3","epitope","label")
input.teim.data<-input.teim.data[!duplicated(input.teim.data),]
input.teim.data.lack<-input.teim.data[!(nchar(input.teim.data$cdr3)<=20),]
input.teim.data<-input.teim.data[nchar(input.teim.data$cdr3)<=20,]
teim.out<-read.table('teim_immrep25_new_random_neg_seq_level_binding.csv',header = T,sep = ',')
teim.out.lack<-cbind(input.teim.data.lack[,1:2],binding=0)
teim.out<-rbind(teim.out[,2:4],teim.out.lack)
input.teim.data<-rbind(input.teim.data,input.teim.data.lack)
stopifnot(sum(teim.out$cdr3!=input.teim.data$cdr3)==0)
stopifnot(sum(teim.out$epitope!=input.teim.data$epitope)==0)
get_partial_auc01(epitope = teim.out$epitope,label = input.teim.data$label,prediction = teim.out$binding)

#ATM-TCR
input.pmtnetv2<-read.table("immrep25_pmtnetv2_input_new_random_neg.tsv",sep = '\t',header = T)
atmtcr_data<-data.frame(input.pmtnetv2$peptide,input.pmtnetv2$cdr3b,input.pmtnetv2$label) #no header
atmtcr_data<-atmtcr_data[!duplicated(atmtcr_data),]
atmtcr_data.lack<-atmtcr_data[!(nchar(atmtcr_data$input.pmtnetv2.cdr3b)<=20),]
atmtcr_data<-atmtcr_data[nchar(atmtcr_data$input.pmtnetv2.cdr3b)<=20,]
atmtcr_data<-atmtcr_data[order(atmtcr_data$input.pmtnetv2.peptide,atmtcr_data$input.pmtnetv2.cdr3b),]

atmtcr.result<-read.table("pred_original_immrep25_atmtcr_input_new_random_neg.csv",header = F,sep = '\t')
atmtcr.result<-atmtcr.result[order(atmtcr.result$V1,atmtcr.result$V2),]
rownames(atmtcr.result)<-NULL

atmtcr.result<-rbind(atmtcr.result[,c(1,2,5)],data.frame(V1=atmtcr_data.lack$input.pmtnetv2.peptide,V2=atmtcr_data.lack$input.pmtnetv2.cdr3b,V5=rep(0,nrow(atmtcr_data.lack))))
atmtcr_data<-rbind(atmtcr_data,atmtcr_data.lack)

stopifnot(sum(atmtcr.result$V1!=atmtcr_data$input.pmtnetv2.peptide)==0)
stopifnot(sum(atmtcr.result$V2!=atmtcr_data$input.pmtnetv2.cdr3b)==0)

get_partial_auc01(epitope = atmtcr.result$V1,label = atmtcr_data$input.pmtnetv2.label,prediction = atmtcr.result$V5)

#Imrex
input.pmtnetv2<-read.table("immrep25_pmtnetv2_input_new_random_neg.tsv",header = T,sep='\t')
imrex<-data.frame(cdr3=input.pmtnetv2$cdr3b,antigen.epitope=input.pmtnetv2$peptide,label=input.pmtnetv2$label) #sep=；
imrex<-imrex[!duplicated(imrex[,1:2]),]
imrex.result<-read.table("immrep25_imrex_new_random_neg_output.csv",header = T,sep = ',')
imrex.lack<-imrex[!((imrex$antigen.epitope %in% imrex.result$antigen.epitope) & (imrex$cdr3 %in% imrex.result$cdr3)),]
imrex<-imrex[imrex$antigen.epitope %in% imrex.result$antigen.epitope,]
imrex<-imrex[imrex$cdr3 %in% imrex.result$cdr3,]
imrex<-imrex[order(imrex$antigen.epitope,imrex$antigen.epitope),]
imrex.result<-imrex.result[order(imrex.result$antigen.epitope,imrex.result$antigen.epitope),]
imrex<-rbind(imrex,imrex.lack)
imrex.result<-rbind(imrex.result,cbind(imrex.lack[,1:2],prediction_score=rep(0,nrow(imrex.lack))))
stopifnot(sum(unlist(imrex[,1:2])!=unlist(imrex.result[,1:2]))==0)
get_partial_auc01(epitope = imrex.result$antigen.epitope,label = imrex$label,prediction = imrex.result$prediction_score)

#ERGO
ergo.1<-read.table("ergo_au_mc_new_random_neg_results.csv",header = T,sep = ',')
ergo.2<-read.table("ergo_au_vd_new_random_neg_results.csv",header = T,sep = ',')
ergo.3<-read.table("ergo_ls_mc_new_random_neg_results.csv",header = T,sep = ',')
ergo.4<-read.table("ergo_ls_vd_new_random_neg_results.csv",header = T,sep = ',')
ergo.1<-ergo.1[,c(1:2,8:10)] |> arrange(TRA,TRB,Peptide,MHC)
ergo.2<-ergo.2[,c(1:2,8:10)] |> arrange(TRA,TRB,Peptide,MHC)
ergo.3<-ergo.3[,c(1:2,8:10)] |> arrange(TRA,TRB,Peptide,MHC)
ergo.4<-ergo.4[,c(1:2,8:10)] |> arrange(TRA,TRB,Peptide,MHC)
stopifnot(sum(unlist(ergo.1[,1:4])!=unlist(ergo.2[,1:4]))==0)
stopifnot(sum(unlist(ergo.1[,1:4])!=unlist(ergo.3[,1:4]))==0)
stopifnot(sum(unlist(ergo.1[,1:4])!=unlist(ergo.4[,1:4]))==0)
input.pmtnetv2<-read.table("immrep25_pmtnetv2_input_new_random_neg.tsv",header = T,sep='\t')
ergo2<-data.frame(TRA=input.pmtnetv2$cdr3a,TRB=input.pmtnetv2$cdr3b,
                  Peptide=input.pmtnetv2$peptide,
                  MHC=paste('HLA-',unlist(lapply(strsplit(input.pmtnetv2$mhca,split = ":"),'[[',1)),sep = ""),
                  label=input.pmtnetv2$label)
ergo2<-ergo2[!duplicated(ergo2[,1:4]),]
ergo2<- ergo2 |> arrange(TRA,TRB,Peptide,MHC)
stopifnot(sum(unlist(ergo2[,1:4])!=unlist(ergo.1[,1:4]))==0)
get_partial_auc01(epitope = ergo.1$Peptide,label = ergo2$label,prediction = ergo.1$Score)
get_partial_auc01(epitope = ergo.2$Peptide,label = ergo2$label,prediction = ergo.2$Score)
get_partial_auc01(epitope = ergo.3$Peptide,label = ergo2$label,prediction = ergo.3$Score)
get_partial_auc01(epitope = ergo.4$Peptide,label = ergo2$label,prediction = ergo.4$Score)

#NetTCR
get_partial_auc01_nettcr<-function(epitope,label,prediction){
  auc.par.all<-c()
  auc_min<-0.1*0.1/2
  auc_max<-0.1
  for (i in unique(epitope)) {
    if(length(unique(label[epitope==i]))==1){next}
    auc.par.tmp<-pROC::auc(response=label[epitope==i],predictor=prediction[epitope==i],
                           direction="<",
                           partial.auc=c(1,0.9),partial.auc.focus="specificity", partial.auc.correct=F,quiet=T)
    auc.par.tmp<-((as.numeric(auc.par.tmp) - auc_min) / (auc_max - auc_min) + 1)*0.5
    auc.par.all<-c(auc.par.all,auc.par.tmp)
  }
  return(mean(auc.par.all))
}
nettcr.input<-read.table("immrep25_nettcr_input_new_random_neg.csv",header = T,sep = ',')
nettcr<-read.table("nettcr_predictions_new_random_neg.csv",header = T,sep = ',')
nettcr<-nettcr[,1:8]

stopifnot(sum(nettcr.input$peptide!=nettcr$peptide)==0)
stopifnot(sum(unlist(nettcr.input[,1:7])!=unlist(nettcr[,1:7]))==0)

get_partial_auc01_nettcr(epitope = nettcr$peptide,label = nettcr.input$label,prediction = nettcr$prediction)

nettcr.lack<-read.table("immrep25_nettcr_input_new_random_neg_lack.csv",header = T,sep = ',')
nettcr.lack$prediction<-0
nettcr.input.lack<-nettcr.lack[,1:8]
nettcr.lack<-nettcr.lack[,c(1:7,9)]
nettcr.full<-rbind(nettcr,nettcr.lack)
nettcr.input.full<-rbind(nettcr.input,nettcr.input.lack)
stopifnot(sum(nettcr.input.full$peptide!=nettcr.full$peptide)==0)
stopifnot(sum(unlist(nettcr.input.full[,1:7])!=unlist(nettcr.full[,1:7]))==0)

get_partial_auc01(epitope = nettcr.full$peptide,label = nettcr.input.full$label,prediction = nettcr.full$prediction)

#EPACT
input.pmtnetv2<-read.table("immrep25_pmtnetv2_input_new_random_neg.tsv",header = T,sep='\t')
epact<-data.frame(CDR3.alpha.aa=input.pmtnetv2$cdr3a,
                  CDR3.beta.aa=input.pmtnetv2$cdr3b,
                  Epitope.peptide=input.pmtnetv2$peptide,
                  MHC=input.pmtnetv2$mhca,
                  label=input.pmtnetv2$label)

epact<-epact[!duplicated(epact[,1:4]),]
epact<- epact |> arrange(CDR3.alpha.aa,    CDR3.beta.aa, Epitope.peptide,     MHC)
epact.ori.result<-read.table("epact-origin-wenqi-v2.csv",header = T,sep=',')
epact.ori.result<- epact.ori.result |> arrange(CDR3.alpha.aa,    CDR3.beta.aa, Epitope.peptide,     MHC)
get_partial_auc01(epitope = epact.ori.result$Epitope.peptide,label = epact$label,prediction = epact.ori.result$Pred)

#TULIP
tmp<- stream_in(file("tulip-origin-wenqi-v2.jsonl"))
tulip.ori.result <- tmp |> unnest(cols = c( CDR3a ,      CDR3b, peptide,    rank ))
tulip.ori.result$pcdr<-paste(tulip.ori.result$peptide,tulip.ori.result$CDR3a,tulip.ori.result$CDR3b,sep='_')
tulip.ori.result<-tulip.ori.result[!duplicated(tulip.ori.result$pcdr),]

input.pmtnetv2<-read.table("immrep25_pmtnetv2_input_new_random_neg.tsv",header = T,sep='\t')
tulip<-data.frame(CDR3b=input.pmtnetv2$cdr3b,	CDR3a=input.pmtnetv2$cdr3a,	peptide=input.pmtnetv2$peptide,	MHC=paste0("HLA-",input.pmtnetv2$mhca),label=input.pmtnetv2$label)
tulip<-tulip[!duplicated(tulip[,1:3]),]
tulip$pcdr<-paste(tulip$peptide,tulip$CDR3a,tulip$CDR3b,sep='_')

tulip.ori.result<-tulip.ori.result[tulip.ori.result$pcdr %in% tulip$pcdr,]
tulip<-tulip[tulip$pcdr %in% tulip.ori.result$pcdr,]

tulip.ori.result<- tulip.ori.result |> arrange(peptide, CDR3a, CDR3b)
tulip.tao.result<- tulip.tao.result |> arrange(peptide, CDR3a, CDR3b)
tulip<- tulip |> arrange(peptide, CDR3a, CDR3b)

stopifnot(identical(tulip$peptide,tulip.ori.result$peptide))
stopifnot(identical(tulip$peptide,tulip.tao.result$peptide))
stopifnot(identical(tulip$CDR3a,tulip.ori.result$CDR3a))
stopifnot(identical(tulip$CDR3a,tulip.tao.result$CDR3a))
stopifnot(identical(tulip$CDR3b,tulip.ori.result$CDR3b))
stopifnot(identical(tulip$CDR3b,tulip.tao.result$CDR3b))

get_partial_auc01(epitope = tulip$peptide,label = tulip$label,prediction = 2000-tulip.ori.result$rank)
get_partial_auc01(epitope = tulip$peptide,label = tulip$label,prediction = 2000-tulip.tao.result$rank)

#unseen epitope and HLA
#prepare input data
training<-read.table("training_individual_1.txt",header = T,sep = "\t")
training_peptide<-unique(training$peptide)

immrep25.iedb<-read.table("iedb_positives.csv",header = T,sep=',')
immrep25.iedb$dataset<-'immrep25iedb'
immrep25.vdjdb<-read.table("vdjdb_positives.csv",header = T,sep=',',comment.char = "", na.strings = c("", "*"), stringsAsFactors = FALSE)
immrep25.vdjdb$dataset<-'immrep25vdjdb'
immrep25<-rbind(immrep25.iedb,immrep25.vdjdb)

immrep25.unique.epitope<-immrep25[apply(stringdist::stringdistmatrix(immrep25$Peptide,training_peptide),1,min)>2,]
immrep25.unique.epitope<-immrep25.unique.epitope[immrep25.unique.epitope$just_10X=='False',]
immrep25.unique.epitope<-immrep25.unique.epitope[!grepl("HLA-E*",immrep25.unique.epitope$HLA,fixed = T),]
immrep25.unique.epitope$HLA[immrep25.unique.epitope$HLA=="HLA-A*02:01:48"]<-"HLA-A*02:01"
immrep25.unique.epitope$HLA[immrep25.unique.epitope$HLA=="HLA-B*44:05:01"]<-"HLA-B*44:05"
immrep25.unique.epitope$HLA2<-substr(immrep25.unique.epitope$HLA, 5, nchar(immrep25.unique.epitope$HLA))
immrep25.unique.epitope<-immrep25.unique.epitope[!(immrep25.unique.epitope$HLA2 %in% training$mhca),]

get_partial_auc01<-function(epitope,label,prediction){
  auc.par.all<-c()
  auc_min<-0.1*0.1/2
  auc_max<-0.1
  for (i in unique(epitope)) {
    auc.par.tmp<-pROC::auc(response=label[epitope==i],predictor=prediction[epitope==i],
                           direction="<",
                           partial.auc=c(1,0.9),partial.auc.focus="specificity", partial.auc.correct=F,quiet=T)
    auc.par.tmp<-((as.numeric(auc.par.tmp) - auc_min) / (auc_max - auc_min) + 1)*0.5
    auc.par.all<-c(auc.par.all,auc.par.tmp)
  }
  return(mean(auc.par.all))
}

#pmtnetv2
pmtnetv2<-read.table("immrep25_pmtnetv2_output_all_mhc_new_random_neg.tsv",header = T,sep=",")
pmtnetv2<-pmtnetv2[pmtnetv2$mhca %in% immrep25.unique.epitope$HLA2,]
pmtnetv2$pcdr3b<-paste(pmtnetv2$peptide,pmtnetv2$cdr3b,sep='_')
pmtnetv2$pcdr3ab<-paste(pmtnetv2$peptide,pmtnetv2$cdr3a,pmtnetv2$cdr3b,sep='_')

get_partial_auc01(epitope = pmtnetv2$peptide,label = pmtnetv2$label,prediction = pmtnetv2$avg_rank_asc)

#pmtnetv1
pmtnetv1.out<-read.table('pMTnet_immrep25_new_random_neg.csv',header = T,sep=',')
input.pmtnetv2<-read.table("immrep25_pmtnetv2_input_new_random_neg.tsv",sep = '\t',header = T)
pmtnetv1.input<-input.pmtnetv2[,c("cdr3b","peptide","mhca","label")]
pmtnetv1.input<-pmtnetv1.input[!duplicated(pmtnetv1.input),]

stopifnot(identical(pmtnetv1.input$cdr3b,pmtnetv1.out$CDR3))
stopifnot(identical(pmtnetv1.input$peptide,pmtnetv1.out$Antigen))
stopifnot(identical(pmtnetv1.input$mhca,pmtnetv1.out$HLA))

pmtnetv1.input<-pmtnetv1.input[pmtnetv1.input$mhca %in% immrep25.unique.epitope$HLA2,]
pmtnetv1.out<-pmtnetv1.out[pmtnetv1.out$HLA %in% immrep25.unique.epitope$HLA2,]

get_partial_auc01(epitope = pmtnetv1.out$Antigen,label = pmtnetv1.input$label,prediction = 1-pmtnetv1.out$Rank)

#TEIM
input.pmtnetv2<-read.table("immrep25_pmtnetv2_input_new_random_neg.tsv",sep = '\t',header = T)
input.teim.data<-input.pmtnetv2[,c("cdr3b","peptide","label")]
colnames(input.teim.data)<-c("cdr3","epitope","label")
input.teim.data<-input.teim.data[!duplicated(input.teim.data),]
input.teim.data.lack<-input.teim.data[!(nchar(input.teim.data$cdr3)<=20),]
input.teim.data<-input.teim.data[nchar(input.teim.data$cdr3)<=20,]
teim.out<-read.table('teim_immrep25_new_random_neg_seq_level_binding.csv',header = T,sep = ',')

teim.out.lack<-cbind(input.teim.data.lack[,1:2],binding=0)
teim.out<-rbind(teim.out[,2:4],teim.out.lack)
input.teim.data<-rbind(input.teim.data,input.teim.data.lack)
teim.out$pcdr3b<-paste(teim.out$epitope,teim.out$cdr3,sep='_')
input.teim.data$pcdr3b<-paste(input.teim.data$epitope,input.teim.data$cdr3,sep='_')
teim.out<-teim.out[(teim.out$pcdr3b %in% pmtnetv2$pcdr3b),]
input.teim.data<-input.teim.data[(input.teim.data$pcdr3b %in% pmtnetv2$pcdr3b),]

stopifnot(sum(teim.out$cdr3!=input.teim.data$cdr3)==0)
stopifnot(sum(teim.out$epitope!=input.teim.data$epitope)==0)

get_partial_auc01(epitope = teim.out$epitope,label = input.teim.data$label,prediction = teim.out$binding)

rm(list=setdiff(ls(),c("get_partial_auc01","training","pmtnetv2")))

#ATM-TCR
input.pmtnetv2<-read.table("immrep25_pmtnetv2_input_new_random_neg.tsv",sep = '\t',header = T)
atmtcr_data<-data.frame(input.pmtnetv2$peptide,input.pmtnetv2$cdr3b,input.pmtnetv2$label) #no header
atmtcr_data<-atmtcr_data[!duplicated(atmtcr_data),]
atmtcr_data.lack<-atmtcr_data[!(nchar(atmtcr_data$input.pmtnetv2.cdr3b)<=20),]
atmtcr_data<-atmtcr_data[nchar(atmtcr_data$input.pmtnetv2.cdr3b)<=20,]
atmtcr_data<-atmtcr_data[order(atmtcr_data$input.pmtnetv2.peptide,atmtcr_data$input.pmtnetv2.cdr3b),]

atmtcr.result<-read.table("pred_original_immrep25_atmtcr_input_new_random_neg.csv",header = F,sep = '\t')
atmtcr.result<-atmtcr.result[order(atmtcr.result$V1,atmtcr.result$V2),]
rownames(atmtcr.result)<-NULL

atmtcr.result<-rbind(atmtcr.result[,c(1,2,5)],data.frame(V1=atmtcr_data.lack$input.pmtnetv2.peptide,V2=atmtcr_data.lack$input.pmtnetv2.cdr3b,V5=rep(0,nrow(atmtcr_data.lack))))
atmtcr_data<-rbind(atmtcr_data,atmtcr_data.lack)

atmtcr.result$pcdr3b<-paste(atmtcr.result$V1,atmtcr.result$V2,sep='_')
atmtcr_data$pcdr3b<-paste(atmtcr_data$input.pmtnetv2.peptide,atmtcr_data$input.pmtnetv2.cdr3b,sep='_')

atmtcr.result<-atmtcr.result[atmtcr.result$pcdr3b %in% pmtnetv2$pcdr3b,]
atmtcr_data<-atmtcr_data[atmtcr_data$pcdr3b %in% pmtnetv2$pcdr3b,]

stopifnot(sum(atmtcr.result$V1!=atmtcr_data$input.pmtnetv2.peptide)==0)
stopifnot(sum(atmtcr.result$V2!=atmtcr_data$input.pmtnetv2.cdr3b)==0)

get_partial_auc01(epitope = atmtcr.result$V1,label = atmtcr_data$input.pmtnetv2.label,prediction = atmtcr.result$V5)

rm(list=setdiff(ls(),c("get_partial_auc01","training","pmtnetv2")))

#Imrex
input.pmtnetv2<-read.table("immrep25_pmtnetv2_input_new_random_neg.tsv",header = T,sep='\t')
imrex<-data.frame(cdr3=input.pmtnetv2$cdr3b,antigen.epitope=input.pmtnetv2$peptide,label=input.pmtnetv2$label) #sep=；
imrex<-imrex[!duplicated(imrex[,1:2]),]

imrex.result<-read.table("immrep25_imrex_new_random_neg_output.csv",header = T,sep = ',')

imrex.lack<-imrex[!((imrex$antigen.epitope %in% imrex.result$antigen.epitope) & (imrex$cdr3 %in% imrex.result$cdr3)),]
imrex<-imrex[imrex$antigen.epitope %in% imrex.result$antigen.epitope,]
imrex<-imrex[imrex$cdr3 %in% imrex.result$cdr3,]

imrex<-imrex[order(imrex$antigen.epitope,imrex$antigen.epitope),]
imrex.result<-imrex.result[order(imrex.result$antigen.epitope,imrex.result$antigen.epitope),]

imrex<-rbind(imrex,imrex.lack)
imrex.result<-rbind(imrex.result,cbind(imrex.lack[,1:2],prediction_score=rep(0,nrow(imrex.lack))))

imrex$pcdr3<-paste(imrex$antigen.epitope,imrex$cdr3,sep='_')
imrex.result$pcdr3<-paste(imrex.result$antigen.epitope,imrex.result$cdr3,sep='_')

imrex<-imrex[imrex$pcdr3 %in% pmtnetv2$pcdr3b,]
imrex.result<-imrex.result[imrex.result$pcdr3 %in% pmtnetv2$pcdr3b,]

stopifnot(sum(unlist(imrex[,1:2])!=unlist(imrex.result[,1:2]))==0)

get_partial_auc01(epitope = imrex.result$antigen.epitope,label = imrex$label,prediction = imrex.result$prediction_score)

rm(list=setdiff(ls(),c("get_partial_auc01","training","pmtnetv2")))

#ERGO
ergo.1<-read.table("ergo_au_mc_new_random_neg_results.csv",header = T,sep = ',')
ergo.2<-read.table("ergo_au_vd_new_random_neg_results.csv",header = T,sep = ',')
ergo.3<-read.table("ergo_ls_mc_new_random_neg_results.csv",header = T,sep = ',')
ergo.4<-read.table("ergo_ls_vd_new_random_neg_results.csv",header = T,sep = ',')

ergo.1<-ergo.1[,c(1:2,8:10)] |> arrange(TRA,TRB,Peptide,MHC)
ergo.2<-ergo.2[,c(1:2,8:10)] |> arrange(TRA,TRB,Peptide,MHC)
ergo.3<-ergo.3[,c(1:2,8:10)] |> arrange(TRA,TRB,Peptide,MHC)
ergo.4<-ergo.4[,c(1:2,8:10)] |> arrange(TRA,TRB,Peptide,MHC)

stopifnot(sum(unlist(ergo.1[,1:4])!=unlist(ergo.2[,1:4]))==0)
stopifnot(sum(unlist(ergo.1[,1:4])!=unlist(ergo.3[,1:4]))==0)
stopifnot(sum(unlist(ergo.1[,1:4])!=unlist(ergo.4[,1:4]))==0)

input.pmtnetv2<-read.table("immrep25_pmtnetv2_input_new_random_neg.tsv",header = T,sep='\t')
ergo2<-data.frame(TRA=input.pmtnetv2$cdr3a,TRB=input.pmtnetv2$cdr3b,
                  Peptide=input.pmtnetv2$peptide,
                  MHC=paste('HLA-',unlist(lapply(strsplit(input.pmtnetv2$mhca,split = ":"),'[[',1)),sep = ""),
                  label=input.pmtnetv2$label)
ergo2<-ergo2[!duplicated(ergo2[,1:4]),]

ergo2<- ergo2 |> arrange(TRA,TRB,Peptide,MHC)

ergo.1$pcdr3<-paste(ergo.1$Peptide,ergo.1$TRA,ergo.1$TRB,sep='_')
ergo.2$pcdr3<-paste(ergo.2$Peptide,ergo.2$TRA,ergo.2$TRB,sep='_')
ergo.3$pcdr3<-paste(ergo.3$Peptide,ergo.3$TRA,ergo.3$TRB,sep='_')
ergo.4$pcdr3<-paste(ergo.4$Peptide,ergo.4$TRA,ergo.4$TRB,sep='_')

ergo.1<-ergo.1[ergo.1$pcdr3 %in% pmtnetv2$pcdr3ab,]
ergo.2<-ergo.2[ergo.2$pcdr3 %in% pmtnetv2$pcdr3ab,]
ergo.3<-ergo.3[ergo.3$pcdr3 %in% pmtnetv2$pcdr3ab,]
ergo.4<-ergo.4[ergo.4$pcdr3 %in% pmtnetv2$pcdr3ab,]

stopifnot(sum(unlist(ergo.1[,1:4])!=unlist(ergo.2[,1:4]))==0)
stopifnot(sum(unlist(ergo.1[,1:4])!=unlist(ergo.3[,1:4]))==0)
stopifnot(sum(unlist(ergo.1[,1:4])!=unlist(ergo.4[,1:4]))==0)

ergo2$pcdr3<-paste(ergo2$Peptide,ergo2$TRA,ergo2$TRB,sep='_')
ergo2<-ergo2[ergo2$pcdr3 %in% pmtnetv2$pcdr3ab,]

stopifnot(sum(unlist(ergo2[,1:4])!=unlist(ergo.1[,1:4]))==0)

get_partial_auc01(epitope = ergo.1$Peptide,label = ergo2$label,prediction = ergo.1$Score)
get_partial_auc01(epitope = ergo.2$Peptide,label = ergo2$label,prediction = ergo.2$Score)
get_partial_auc01(epitope = ergo.3$Peptide,label = ergo2$label,prediction = ergo.3$Score)
get_partial_auc01(epitope = ergo.4$Peptide,label = ergo2$label,prediction = ergo.4$Score)

rm(list=setdiff(ls(),c("get_partial_auc01","training","pmtnetv2")))

#NetTCR
get_partial_auc01_nettcr<-function(epitope,label,prediction){
  auc.par.all<-c()
  auc_min<-0.1*0.1/2
  auc_max<-0.1
  for (i in unique(epitope)) {
    if(length(unique(label[epitope==i]))==1){next}
    auc.par.tmp<-pROC::auc(response=label[epitope==i],predictor=prediction[epitope==i],
                           direction="<",
                           partial.auc=c(1,0.9),partial.auc.focus="specificity", partial.auc.correct=F,quiet=T)
    auc.par.tmp<-((as.numeric(auc.par.tmp) - auc_min) / (auc_max - auc_min) + 1)*0.5
    auc.par.all<-c(auc.par.all,auc.par.tmp)
  }
  return(mean(auc.par.all))
}

nettcr.input<-read.table("immrep25_nettcr_input_new_random_neg.csv",header = T,sep = ',')
nettcr<-read.table("nettcr_predictions_new_random_neg.csv",header = T,sep = ',')
nettcr<-nettcr[,1:8]

stopifnot(sum(nettcr.input$peptide!=nettcr$peptide)==0)
stopifnot(sum(unlist(nettcr.input[,1:7])!=unlist(nettcr[,1:7]))==0)

get_partial_auc01_nettcr(epitope = nettcr$peptide,label = nettcr.input$label,prediction = nettcr$prediction)

nettcr.lack<-read.table("immrep25_nettcr_input_new_random_neg_lack.csv",header = T,sep = ',')
nettcr.lack$prediction<-0

nettcr.input.lack<-nettcr.lack[,1:8]
nettcr.lack<-nettcr.lack[,c(1:7,9)]

nettcr.full<-rbind(nettcr,nettcr.lack)
nettcr.input.full<-rbind(nettcr.input,nettcr.input.lack)

nettcr.full$pcdr3ab<-paste(nettcr.full$peptide,nettcr.full$A3,nettcr.full$B3,sep='_')
nettcr.input.full$pcdr3ab<-paste(nettcr.input.full$peptide,nettcr.input.full$A3,nettcr.input.full$B3,sep='_')

nettcr.full<-nettcr.full[nettcr.full$pcdr3ab %in% pmtnetv2$pcdr3ab,]
nettcr.input.full<-nettcr.input.full[nettcr.input.full$pcdr3ab %in% pmtnetv2$pcdr3ab,]

nettcr.full<-nettcr.full[!duplicated(nettcr.full),]
nettcr.input.full<-nettcr.input.full[!duplicated(nettcr.input.full),]

stopifnot(sum(nettcr.input.full$peptide!=nettcr.full$peptide)==0)
stopifnot(sum(unlist(nettcr.input.full[,1:7])!=unlist(nettcr.full[,1:7]))==0)

get_partial_auc01(epitope = nettcr.full$peptide,label = nettcr.input.full$label,prediction = nettcr.full$prediction)

rm(list=setdiff(ls(),c("get_partial_auc01","training","pmtnetv2")))

#EPACT
input.pmtnetv2<-read.table("immrep25_pmtnetv2_input_new_random_neg.tsv",header = T,sep='\t')
epact<-data.frame(CDR3.alpha.aa=input.pmtnetv2$cdr3a,
                  CDR3.beta.aa=input.pmtnetv2$cdr3b,
                  Epitope.peptide=input.pmtnetv2$peptide,
                  MHC=input.pmtnetv2$mhca,
                  label=input.pmtnetv2$label)

epact<-epact[!duplicated(epact[,1:4]),]
epact<- epact |> arrange(CDR3.alpha.aa,    CDR3.beta.aa, Epitope.peptide,     MHC)
epact.ori.result<-read.table("epact-origin-wenqi-v2.csv",header = T,sep=',')
epact.ori.result<- epact.ori.result |> arrange(CDR3.alpha.aa,    CDR3.beta.aa, Epitope.peptide,     MHC)
epact.ori.result<-epact.ori.result[epact.ori.result$MHC %in% pmtnetv2$mhca,]
epact<-epact[epact$MHC %in% pmtnetv2$mhca,]
stopifnot(identical(epact$Epitope.peptide,epact.ori.result$Epitope.peptide))
stopifnot(identical(epact$CDR3.beta.aa,epact.ori.result$CDR3.beta.aa))
stopifnot(identical(epact$CDR3.alpha.aa,epact.ori.result$CDR3.alpha.aa))

get_partial_auc01(epitope = epact.ori.result$Epitope.peptide,label = epact$label,prediction = epact.ori.result$Pred)

rm(list=setdiff(ls(),c("get_partial_auc01","training","pmtnetv2")))

#TULIP
tmp<- stream_in(file("tulip-origin-wenqi-v2.jsonl"))
tulip.ori.result <- tmp |> unnest(cols = c( CDR3a ,      CDR3b, peptide,    rank ))
tulip.ori.result$pcdr<-paste(tulip.ori.result$peptide,tulip.ori.result$CDR3a,tulip.ori.result$CDR3b,sep='_')
tulip.ori.result<-tulip.ori.result[!duplicated(tulip.ori.result$pcdr),]

input.pmtnetv2<-read.table("immrep25_pmtnetv2_input_new_random_neg.tsv",header = T,sep='\t')
tulip<-data.frame(CDR3b=input.pmtnetv2$cdr3b,	CDR3a=input.pmtnetv2$cdr3a,	peptide=input.pmtnetv2$peptide,	MHC=paste0("HLA-",input.pmtnetv2$mhca),label=input.pmtnetv2$label)
tulip<-tulip[!duplicated(tulip[,1:3]),]
tulip$pcdr<-paste(tulip$peptide,tulip$CDR3a,tulip$CDR3b,sep='_')

tulip.ori.result<-tulip.ori.result[tulip.ori.result$pcdr %in% tulip$pcdr,]
tulip<-tulip[tulip$pcdr %in% tulip.ori.result$pcdr,]

tulip.ori.result<- tulip.ori.result |> arrange(peptide, CDR3a, CDR3b)
tulip<- tulip |> arrange(peptide, CDR3a, CDR3b)

stopifnot(identical(tulip$peptide,tulip.ori.result$peptide))
stopifnot(identical(tulip$CDR3a,tulip.ori.result$CDR3a))
stopifnot(identical(tulip$CDR3b,tulip.ori.result$CDR3b))

tulip<-tulip[tulip$pcdr %in% pmtnetv2$pcdr3ab,]
tulip.ori.result<-tulip.ori.result[tulip.ori.result$pcdr %in% pmtnetv2$pcdr3ab,]

stopifnot(identical(tulip$peptide,tulip.ori.result$peptide))
stopifnot(identical(tulip$CDR3a,tulip.ori.result$CDR3a))
stopifnot(identical(tulip$CDR3b,tulip.ori.result$CDR3b))

get_partial_auc01(epitope = tulip$peptide,label = tulip$label,prediction = 2000-tulip.ori.result$rank)

# Sup. File 2 Fig. 7
pmtnetv1<-read.table('pMTnet_v1_prame_rootpath_validated_prediction.csv',header = T,sep=',')
teim<-read.table('teim_prame_rootpath_validated_sequence_level_binding.csv',header = T,sep=',')
atm<-read.table('pred_original_input_rootpath_validated_atmtcr.csv',header = F,sep='\t')
imrex<-read.table('output_rootpath_validated_imrex.csv',header = T,sep=',')
ergo2ae.mac<-read.table('ergo2-ae-mc-prame-rootpath-validated.csv',header = T,sep=',')
ergo2ae.vdj<-read.table('ergo2-ae-vdj-prame-rootpath-validated.csv',header = T,sep=',')
ergo2ls.mac<-read.table('ergo2-ls-mc-prame-rootpath-validated.csv',header = T,sep=',')
ergo2ls.vdj<-read.table('ergo2-ls-vdj-prame-rootpath-validated.csv',header = T,sep=',')

epact.immunocore<-read.table('preds-epact-origin_immunocore.csv',header = T,sep=',')
nettcr.immunocore<-read.table('immunocore_result.csv',header = T,sep=',')
lines <- readLines("ranking-tulip-origin_immunocore.jsonl")
tulip.immunocore <- lapply(lines, fromJSON)
tulip.immunocore <- data.frame(peptide=rep(tulip.immunocore[[1]]$peptide),cdr3a=tulip.immunocore[[1]]$CDR3a,cdr3b=tulip.immunocore[[1]]$CDR3b,rank=tulip.immunocore[[1]]$rank)

candidates<-read.table("immunocore_a79b71_predictions_after_5shot_learning.csv",header = T,sep = ',')
candidates<-candidates[order(candidates$cdr3a,candidates$cdr3b),]
epact.immunocore<-epact.immunocore[order(epact.immunocore$CDR3.alpha.aa,epact.immunocore$CDR3.beta.aa),]
nettcr.immunocore<-nettcr.immunocore[order(nettcr.immunocore$cdr3a,nettcr.immunocore$cdr3b),]
tulip.immunocore<-tulip.immunocore[order(tulip.immunocore$cdr3a,tulip.immunocore$cdr3b),]

stopifnot(sum(paste(epact.immunocore$CDR3.alpha.aa,epact.immunocore$CDR3.beta.aa,sep = '_')!=paste(candidates$cdr3a,candidates$cdr3b,sep='_'))==0)
stopifnot(sum(paste(nettcr.immunocore$cdr3a,nettcr.immunocore$cdr3b,sep = '_')!=paste(candidates$cdr3a,candidates$cdr3b,sep='_'))==0)
stopifnot(sum(paste(tulip.immunocore$cdr3a,tulip.immunocore$cdr3b,sep = '_')!=paste(candidates$cdr3a,candidates$cdr3b,sep='_'))==0)

epact.immunocore$spot<-candidates$elispot
nettcr.immunocore$spot<-candidates$elispot
tulip.immunocore$spot<-candidates$elispot

candidates<-read.table("immunocore_a79b71_predictions_after_5shot_learning.csv",header = T,sep = ',')
candidates<-candidates[order(candidates$cdr3b),]

stopifnot(sum(pmtnetv1$CDR3!=candidates$cdr3b)==0)
stopifnot(sum(teim$cdr3!=candidates$cdr3b)==0)
stopifnot(sum(imrex$cdr3!=candidates$cdr3b)==0)
stopifnot(sum(ergo2ae.mac$TRB!=candidates$cdr3b)==0)
stopifnot(sum(ergo2ae.mac$TRA!=candidates$cdr3a)==0)

stopifnot(sum(ergo2ae.vdj$TRB!=candidates$cdr3b)==0)
stopifnot(sum(ergo2ae.vdj$TRA!=candidates$cdr3a)==0)

stopifnot(sum(ergo2ls.mac$TRB!=candidates$cdr3b)==0)
stopifnot(sum(ergo2ls.mac$TRA!=candidates$cdr3a)==0)

stopifnot(sum(ergo2ls.vdj$TRB!=candidates$cdr3b)==0)
stopifnot(sum(ergo2ls.vdj$TRA!=candidates$cdr3a)==0)

atm<-atm[order(atm$V2),]
stopifnot(sum(atm$V2!=candidates$cdr3b)==0)

pmtnetv1$spot<-candidates$elispot
teim$spot<-candidates$elispot
imrex$spot<-candidates$elispot
ergo2ae.mac$spot<-candidates$elispot
ergo2ae.vdj$spot<-candidates$elispot
ergo2ls.mac$spot<-candidates$elispot
ergo2ls.vdj$spot<-candidates$elispot
atm$spot<-candidates$elispot

data.plot<-data.frame(group=c(rep('pmtnetv1_better',sum(pmtnetv1$Rank<pmtnetv1$Rank[pmtnetv1$CDR3=="CASSWWTGGSAPIYF"])),rep('pmtnetv1_worse',sum(pmtnetv1$Rank>pmtnetv1$Rank[pmtnetv1$CDR3=="CASSWWTGGSAPIYF"])),
                              rep('NetTCR_better',sum(nettcr.immunocore$prediction>nettcr.immunocore$prediction[nettcr.immunocore$spot==122])),rep('NetTCR_worse',sum(nettcr.immunocore$prediction<nettcr.immunocore$prediction[nettcr.immunocore$spot==122])),
                              rep('teim_better',sum(teim$binding>teim$binding[teim$cdr3=="CASSWWTGGSAPIYF"])),rep('teim_worse',sum(teim$binding<teim$binding[teim$cdr3=="CASSWWTGGSAPIYF"])),
                              rep('epact_better',sum(epact.immunocore$Pred>epact.immunocore$Pred[epact.immunocore$spot==122])),rep('epact_worse',sum(epact.immunocore$Pred<epact.immunocore$Pred[epact.immunocore$spot==122])),
                              rep('tulip_better',sum(tulip.immunocore$rank<tulip.immunocore$rank[tulip.immunocore$spot==122])),rep('tulip_worse',sum(tulip.immunocore$rank>tulip.immunocore$rank[epact.immunocore$spot==122])),
                              rep('imrex_better',1),rep('imrex_worse',sum(imrex$prediction_score<imrex$prediction_score[imrex$cdr3=="CASSWWTGGSAPIYF"])),
                              rep('ergo2-ae-mac_better',sum(ergo2ae.mac$Score>ergo2ae.mac$Score[ergo2ae.mac$TRB=="CASSWWTGGSAPIYF"])),rep('ergo2-ae-mac_worse',sum(ergo2ae.mac$Score<ergo2ae.mac$Score[ergo2ae.mac$TRB=="CASSWWTGGSAPIYF"])),
                              rep('ergo2-ae-vdj_better',sum(ergo2ae.vdj$Score>ergo2ae.vdj$Score[ergo2ae.vdj$TRB=="CASSWWTGGSAPIYF"])),rep('ergo2-ae-vdj_worse',sum(ergo2ae.vdj$Score<ergo2ae.vdj$Score[ergo2ae.vdj$TRB=="CASSWWTGGSAPIYF"])),
                              rep('ergo2-ls-mac_better',sum(ergo2ls.mac$Score>ergo2ls.mac$Score[ergo2ls.mac$TRB=="CASSWWTGGSAPIYF"])),rep('ergo2-ls-mac_worse',sum(ergo2ls.mac$Score<ergo2ls.mac$Score[ergo2ls.mac$TRB=="CASSWWTGGSAPIYF"])),
                              rep('ergo2-ls-vdj_better',sum(ergo2ls.vdj$Score>ergo2ls.vdj$Score[ergo2ls.vdj$TRB=="CASSWWTGGSAPIYF"])),rep('ergo2-ls-vdj_worse',sum(ergo2ls.vdj$Score<ergo2ls.vdj$Score[ergo2ls.vdj$TRB=="CASSWWTGGSAPIYF"])),
                              rep('atm_better',sum(atm$V5>atm$V5[atm$V2=="CASSWWTGGSAPIYF"])),rep('atm_worse',sum(atm$V5<atm$V5[atm$V2=="CASSWWTGGSAPIYF"]))
),
spot=c(pmtnetv1$spot[pmtnetv1$Rank<pmtnetv1$Rank[pmtnetv1$CDR3=="CASSWWTGGSAPIYF"]],pmtnetv1$spot[pmtnetv1$Rank>pmtnetv1$Rank[pmtnetv1$CDR3=="CASSWWTGGSAPIYF"]],
       nettcr.immunocore$spot[nettcr.immunocore$prediction>nettcr.immunocore$prediction[nettcr.immunocore$spot==122]],nettcr.immunocore$spot[nettcr.immunocore$prediction<nettcr.immunocore$prediction[nettcr.immunocore$spot==122]],
       teim$spot[teim$binding>teim$binding[teim$cdr3=="CASSWWTGGSAPIYF"]],teim$spot[teim$binding<teim$binding[teim$cdr3=="CASSWWTGGSAPIYF"]],
       epact.immunocore$spot[epact.immunocore$Pred>epact.immunocore$Pred[epact.immunocore$spot==122]],epact.immunocore$spot[epact.immunocore$Pred<epact.immunocore$Pred[epact.immunocore$spot==122]],
       tulip.immunocore$spot[tulip.immunocore$rank<tulip.immunocore$rank[tulip.immunocore$spot==122]],tulip.immunocore$spot[tulip.immunocore$rank>tulip.immunocore$rank[tulip.immunocore$spot==122]],
       0,imrex$spot[imrex$prediction_score<imrex$prediction_score[imrex$cdr3=="CASSWWTGGSAPIYF"]],
       ergo2ae.mac$spot[ergo2ae.mac$Score>ergo2ae.mac$Score[ergo2ae.mac$TRB=="CASSWWTGGSAPIYF"]],ergo2ae.mac$spot[ergo2ae.mac$Score<ergo2ae.mac$Score[ergo2ae.mac$TRB=="CASSWWTGGSAPIYF"]],
       ergo2ae.vdj$spot[ergo2ae.vdj$Score>ergo2ae.vdj$Score[ergo2ae.vdj$TRB=="CASSWWTGGSAPIYF"]],ergo2ae.vdj$spot[ergo2ae.vdj$Score<ergo2ae.vdj$Score[ergo2ae.vdj$TRB=="CASSWWTGGSAPIYF"]],
       ergo2ls.mac$spot[ergo2ls.mac$Score>ergo2ls.mac$Score[ergo2ls.mac$TRB=="CASSWWTGGSAPIYF"]],ergo2ls.mac$spot[ergo2ls.mac$Score<ergo2ls.mac$Score[ergo2ls.mac$TRB=="CASSWWTGGSAPIYF"]],
       ergo2ls.vdj$spot[ergo2ls.vdj$Score>ergo2ls.vdj$Score[ergo2ls.vdj$TRB=="CASSWWTGGSAPIYF"]],ergo2ls.vdj$spot[ergo2ls.vdj$Score<ergo2ls.vdj$Score[ergo2ls.vdj$TRB=="CASSWWTGGSAPIYF"]],
       atm$spot[atm$V5>atm$V5[atm$V2=="CASSWWTGGSAPIYF"]],atm$spot[atm$V5<atm$V5[atm$V2=="CASSWWTGGSAPIYF"]]))

data.plot$group<-factor(data.plot$group,levels = c("pmtnetv1_better","pmtnetv1_worse",
                                                   "NetTCR_better","NetTCR_worse",
                                                   "teim_better","teim_worse",
                                                   "epact_better","epact_worse",
                                                   "atm_better","atm_worse",
                                                   "tulip_better","tulip_worse",
                                                   "imrex_better","imrex_worse",
                                                   "ergo2-ae-vdj_better","ergo2-ae-vdj_worse",
                                                   "ergo2-ae-mac_better","ergo2-ae-mac_worse",
                                                   "ergo2-ls-vdj_better","ergo2-ls-vdj_worse",
                                                   "ergo2-ls-mac_better","ergo2-ls-mac_worse"))

data.plot$group2<-ifelse(grepl('better',data.plot$group),'better','worse')

ggplot(data.plot, aes(x=group, y=spot,fill=group2,color=group2)) +
  geom_bar(stat="summary", fun=median, color=NA,alpha=0.6) +
  stat_summary(fun = median,
               geom = "errorbar",width=0.4,
               fun.max = function(x) median(x) + sd(x) / sqrt(length(x)) * 1.2533,
               fun.min = function(x) median(x) - sd(x) / sqrt(length(x)) * 1.2533)+
  geom_point(position = position_jitter(width = 0.1),size=3,alpha=0.5)+
  geom_hline(yintercept=122, linetype="dashed",alpha=0.7)+
  labs(x = "",y = "# of IFN gamma ELISPOT") +
  theme_classic() +
  theme(axis.text=element_text(size=12), 
        axis.text.x = element_text(angle=30,hjust = 0.9,vjust = 0.9),
        axis.title=element_text(size=14,face="bold"),
        legend.position="none")+
  scale_fill_manual(values = c("#00A087B2","#E64B35B2")) +
  scale_color_manual(values = c("#00A087B2","#E64B35B2"))
