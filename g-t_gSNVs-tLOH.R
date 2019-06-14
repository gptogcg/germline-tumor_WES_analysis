
setwd("/home/mdiaz/Desktop/PROYECTOS/FAMCOLON_11/Germline-Tumor_Analysis_FAMCOLON_11")
library(data.table)
library(openxlsx)
library(readxl)
library(WriteXLS)

samples<-fread("g-t_samples_FAMCOLON_11.txt",sep="\t",header=T,data.table=F)
families<-levels(factor(samples$FAMILIES))[-33]
families_mispaired<-levels(factor(samples$FAMILIES))[33]

#Pipeline for paired samples

##Building of dataframe of samples
df_samples<-data.frame(families,matrix(nrow=length(families),ncol=4))
colnames(df_samples)<-c("Family","Germline","Tumor","Study","Tumor_project")
for (i in 1:length(families)){
    fam<-families[i]
    aux<-which(samples$FAMILIES==fam & samples$MOSTRA=="germinal")
    df_samples$Germline[i]<-samples$SAMPLES[aux]
    aux_2<-which(samples$FAMILIES==fam & samples$MOSTRA=="tumoral")
    df_samples$Tumor[i]<-samples$SAMPLES[aux_2]
    df_samples$Study[i]<-samples$ESTUDI[aux_2]
}

df_samples$Tumor_project[c(1:8,10,12,14,16:17,28:32)]<-"Marcos"
df_samples$Tumor_project[-c(1:8,10,12,14,16:17,28:32)]<-"CNAG"

##Germline variants dataframe loading

###CRC
gSNVs_CRC<-fread("/home/mdiaz/Desktop/PROYECTOS/FAMCOLON_11/germ_FAMCOLON_11/CRC_Famcolon_01.08__EXOMES_6.0.txt",sep="\t",header=T,data.table=F)

SAMPLES_2<-fread("/home/mdiaz/Desktop/PROYECTOS/FAMCOLON_11/germ_FAMCOLON_11/CRC_LIST_samples.txt",sep="\t",header=T,data.table=F)
EXOMES_6.0<-gSNVs_CRC

first_exome_column <- as.numeric(grep("_GT", colnames(EXOMES_6.0))[1])
colnames_tail_exomes <- as.numeric(length(as.numeric(grep("_PL", colnames(EXOMES_6.0)))))
last_exome_column <- as.numeric(grep("_PL", colnames(EXOMES_6.0))[colnames_tail_exomes])
NONE_EXOME_COLUMNS <- 1:(first_exome_column-1)
NONE_EXOME_COLUMNS_TAIL <- (last_exome_column+1):length(colnames(EXOMES_6.0))
n_cols_x_ind <- as.numeric(grep("_GT", colnames(EXOMES_6.0))[2])-first_exome_column
all_fams <- list()
for (ss in 1:length(unique(SAMPLES_2[,"FAMILIES"]))){
    fam <- as.character(unique(SAMPLES_2[,"FAMILIES"])[ss])
    inds <- as.character(SAMPLES_2[which(SAMPLES_2[,"FAMILIES"]==fam),1])
    if (length(inds)>6){
        print("Només està preparat per famíles amb un màxim de 6 persones!!!!")
    } else if (length(inds)<=6){
        all_cols <- list()
        for (xx in 1:length(inds)){
            name <- inds[xx]
            cols <- grep(name, colnames(EXOMES_6.0))
            all_cols[[as.character(name)]] <- cols
        }
        all_cols2 <- do.call(cbind, all_cols)
        INDS_COLS <- as.vector(all_cols2)
        FAM <- EXOMES_6.0[,c(NONE_EXOME_COLUMNS, INDS_COLS, NONE_EXOME_COLUMNS_TAIL)]
        first_exome_column <- as.numeric(grep("_GT", colnames(FAM))[1])
        nonexome_columns <- 1:(first_exome_column-1)
        first_GT <- first_exome_column
        first_DP <- as.numeric(grep("_DP", colnames(FAM))[1])
        if (length(inds)==2){
            col_GT_1 <- first_GT
            col_GT_2 <- first_GT + n_cols_x_ind
            sel <- which( (FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0") & (FAM[,col_GT_2]!="./." & FAM[,col_GT_2]!="0/0"))
        } else if (length(inds)==3){
            col_GT_1 <- first_GT
            col_GT_2 <- first_GT + n_cols_x_ind
            col_GT_3 <- first_GT + (n_cols_x_ind*2)
            sel <- which( (FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0") & (FAM[,col_GT_2]!="./." & FAM[,col_GT_2]!="0/0") & (FAM[,col_GT_3]!="./." & FAM[,col_GT_3]!="0/0"))
        } else if (length(inds)==4){
            col_GT_1 <- first_GT
            col_GT_2 <- first_GT + n_cols_x_ind
            col_GT_3 <- first_GT + (n_cols_x_ind*2)
            col_GT_4 <- first_GT + (n_cols_x_ind*3)
            sel <- which( (FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0") & (FAM[,col_GT_2]!="./." & FAM[,col_GT_2]!="0/0") & (FAM[,col_GT_3]!="./." & FAM[,col_GT_3]!="0/0") & (FAM[,col_GT_4]!="./." & FAM[,col_GT_4]!="0/0"))
        } else if (length(inds)==5){
            col_GT_1 <- first_GT
            col_GT_2 <- first_GT + n_cols_x_ind
            col_GT_3 <- first_GT + (n_cols_x_ind*2)
            col_GT_4 <- first_GT + (n_cols_x_ind*3)
            col_GT_5 <- first_GT + (n_cols_x_ind*4)
            sel <- which( (FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0") & (FAM[,col_GT_2]!="./." & FAM[,col_GT_2]!="0/0") & (FAM[,col_GT_3]!="./." & FAM[,col_GT_3]!="0/0") & (FAM[,col_GT_4]!="./." & FAM[,col_GT_4]!="0/0") & (FAM[,col_GT_5]!="./." & FAM[,col_GT_5]!="0/0"))
        } else if (length(inds)==6){
            col_GT_1 <- first_GT
            col_GT_2 <- first_GT + n_cols_x_ind
            col_GT_3 <- first_GT + (n_cols_x_ind*2)
            col_GT_4 <- first_GT + (n_cols_x_ind*3)
            col_GT_5 <- first_GT + (n_cols_x_ind*4)
            col_GT_6 <- first_GT + (n_cols_x_ind*5)
            sel <- which( (FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0") & (FAM[,col_GT_2]!="./." & FAM[,col_GT_2]!="0/0") & (FAM[,col_GT_3]!="./." & FAM[,col_GT_3]!="0/0") & (FAM[,col_GT_4]!="./." & FAM[,col_GT_4]!="0/0") & (FAM[,col_GT_5]!="./." & FAM[,col_GT_5]!="0/0") & (FAM[,col_GT_6]!="./." & FAM[,col_GT_6]!="0/0"))
        } else if (length(inds)==1){
            col_GT_1 <- first_GT
            sel <- which( (FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0"))
        }
        all_fams[[fam]] <- FAM[sel,] #keeping variants rows as dataframe
    }
}

gSNVs_CRC<-all_fams



###SPS
gSNVs_SPS<-fread("/home/mdiaz/Desktop/PROYECTOS/FAMCOLON_11/germ_FAMCOLON_11/SPS_Famcolon_01.08__EXOMES_6.0.txt",sep="\t",header=T,data.table=F)

SAMPLES_2<-fread("/home/mdiaz/Desktop/PROYECTOS/FAMCOLON_11/germ_FAMCOLON_11/SPS_LIST_samples.txt",sep="\t",header=T,data.table=F)
EXOMES_6.0<-gSNVs_SPS

first_exome_column <- as.numeric(grep("_GT", colnames(EXOMES_6.0))[1])
colnames_tail_exomes <- as.numeric(length(as.numeric(grep("_PL", colnames(EXOMES_6.0)))))
last_exome_column <- as.numeric(grep("_PL", colnames(EXOMES_6.0))[colnames_tail_exomes])
NONE_EXOME_COLUMNS <- 1:(first_exome_column-1)
NONE_EXOME_COLUMNS_TAIL <- (last_exome_column+1):length(colnames(EXOMES_6.0))
n_cols_x_ind <- as.numeric(grep("_GT", colnames(EXOMES_6.0))[2])-first_exome_column
all_fams <- list()
for (ss in 1:length(unique(SAMPLES_2[,"FAMILIES"]))){
    fam <- as.character(unique(SAMPLES_2[,"FAMILIES"])[ss])
    inds <- as.character(SAMPLES_2[which(SAMPLES_2[,"FAMILIES"]==fam),1])
    if (length(inds)>6){
        print("Només està preparat per famíles amb un màxim de 6 persones!!!!")
    } else if (length(inds)<=6){
        all_cols <- list()
        for (xx in 1:length(inds)){
            name <- inds[xx]
            cols <- grep(name, colnames(EXOMES_6.0))
            all_cols[[as.character(name)]] <- cols
        }
        all_cols2 <- do.call(cbind, all_cols)
        INDS_COLS <- as.vector(all_cols2)
        FAM <- EXOMES_6.0[,c(NONE_EXOME_COLUMNS, INDS_COLS, NONE_EXOME_COLUMNS_TAIL)]
        first_exome_column <- as.numeric(grep("_GT", colnames(FAM))[1])
        nonexome_columns <- 1:(first_exome_column-1)
        first_GT <- first_exome_column
        first_DP <- as.numeric(grep("_DP", colnames(FAM))[1])
        if (length(inds)==2){
            col_GT_1 <- first_GT
            col_GT_2 <- first_GT + n_cols_x_ind
            sel <- which( (FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0") & (FAM[,col_GT_2]!="./." & FAM[,col_GT_2]!="0/0"))
        } else if (length(inds)==3){
            col_GT_1 <- first_GT
            col_GT_2 <- first_GT + n_cols_x_ind
            col_GT_3 <- first_GT + (n_cols_x_ind*2)
            sel <- which( (FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0") & (FAM[,col_GT_2]!="./." & FAM[,col_GT_2]!="0/0") & (FAM[,col_GT_3]!="./." & FAM[,col_GT_3]!="0/0"))
        } else if (length(inds)==4){
            col_GT_1 <- first_GT
            col_GT_2 <- first_GT + n_cols_x_ind
            col_GT_3 <- first_GT + (n_cols_x_ind*2)
            col_GT_4 <- first_GT + (n_cols_x_ind*3)
            sel <- which( (FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0") & (FAM[,col_GT_2]!="./." & FAM[,col_GT_2]!="0/0") & (FAM[,col_GT_3]!="./." & FAM[,col_GT_3]!="0/0") & (FAM[,col_GT_4]!="./." & FAM[,col_GT_4]!="0/0"))
        } else if (length(inds)==5){
            col_GT_1 <- first_GT
            col_GT_2 <- first_GT + n_cols_x_ind
            col_GT_3 <- first_GT + (n_cols_x_ind*2)
            col_GT_4 <- first_GT + (n_cols_x_ind*3)
            col_GT_5 <- first_GT + (n_cols_x_ind*4)
            sel <- which( (FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0") & (FAM[,col_GT_2]!="./." & FAM[,col_GT_2]!="0/0") & (FAM[,col_GT_3]!="./." & FAM[,col_GT_3]!="0/0") & (FAM[,col_GT_4]!="./." & FAM[,col_GT_4]!="0/0") & (FAM[,col_GT_5]!="./." & FAM[,col_GT_5]!="0/0"))
        } else if (length(inds)==6){
            col_GT_1 <- first_GT
            col_GT_2 <- first_GT + n_cols_x_ind
            col_GT_3 <- first_GT + (n_cols_x_ind*2)
            col_GT_4 <- first_GT + (n_cols_x_ind*3)
            col_GT_5 <- first_GT + (n_cols_x_ind*4)
            col_GT_6 <- first_GT + (n_cols_x_ind*5)
            sel <- which( (FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0") & (FAM[,col_GT_2]!="./." & FAM[,col_GT_2]!="0/0") & (FAM[,col_GT_3]!="./." & FAM[,col_GT_3]!="0/0") & (FAM[,col_GT_4]!="./." & FAM[,col_GT_4]!="0/0") & (FAM[,col_GT_5]!="./." & FAM[,col_GT_5]!="0/0") & (FAM[,col_GT_6]!="./." & FAM[,col_GT_6]!="0/0"))
        } else if (length(inds)==1){
            col_GT_1 <- first_GT
            sel <- which( (FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0"))
        }
        all_fams[[fam]] <- FAM[sel,] #keeping variants rows as dataframe
    }
}


gSNVs_SPS<-all_fams


###All
gSNVs<-list()
for (i in 1:length(gSNVs_CRC)){
    gSNVs[[i]]<-gSNVs_CRC[[i]]
    names(gSNVs)[i]<-names(gSNVs_CRC)[i]
}
for (i in 1:length(gSNVs_SPS)){
    gSNVs[[(length(gSNVs_CRC)+i)]]<-gSNVs_SPS[[i]]
    names(gSNVs)[(length(gSNVs_CRC)+i)]<-names(gSNVs_SPS)[i]
}


##LOH matrix loading

ALFRED<-fread("/home/mdiaz/Desktop/PROYECTOS/Fran_Supek_Collaboration/LOH_ALFRED/Results_FAMCOLON_11/ALFRED_input_matrix_ver2_transposed.txt",sep="\t",header=T,data.table=F)

ALFRED<-ALFRED[grep(":LOH",ALFRED$Gene),]
ALFRED$Gene<-sapply(strsplit(ALFRED$Gene,":"),"[",1)
new_colnames<-c("Gene","MOS_SPS_14","Ind5","MOS_SPS_11","MOS_SPS_12","MOS_SPS_15","Ind9","CAR_SPS.8","CAR_SPS.9","CAR_SPS.4","CAR_SPS.6","CAR_SPS.7","CAR_SPS.1","CAR_SPS.2","CAR_SPS.3","Ind4","Ind3","fam1","fam3","fam4","Ind15","DONOS_SPS.5","Ind13","Ind12","Ind11","MOS_SPS_13","fam23","fam22","fam20","fam_nova_3","fam_nova_1","fam_nova_4","fam19")
colnames(ALFRED)<-new_colnames


##Pipeline SNV-LOH
gt<-list()
unknown<-integer(0)
unknown_fam<-integer(0)
for (i in 1:length(families)){
    
    #Patient definition
    family<-families[i]
    germline<-df_samples$Germline[i]
    tumor<-df_samples$Tumor[i]
    study<-df_samples$Study[i]
    tumor_project<-df_samples$Tumor_project[i]
    
    #Germline Variants
    g_df<-which(names(gSNVs)==family)
    g_df<-gSNVs[[g_df]]
    
    #Tumor LOHs
    t_loh<-which(colnames(ALFRED)==family)
    t_loh<-ALFRED[,c(1,t_loh)]

    #Cross Variants
    
    g_genes<-as.character(g_df$Gene_Name)
    
    g_catch<-integer(0)
    t_catch<-integer(0)
    for (j in 1:length(g_genes)){
        gene<-g_genes[j]
        
        catch<-integer(0)
        catch<-grep(paste0("\\b",gene,"\\b"),t_loh$Gene)
        
        if (length(catch)==0){
            print(paste0("No está el gen ",gene," en ALFRED!!!!"))
            unknown<-rbind(unknown,gene)
            unknown_fam<-rbind(unknown_fam,family)
        }
        if (length(catch)==1){
            if (t_loh[catch,2]==1){
            g_catch<-rbind(g_catch,g_df[j,])
            t_catch<-rbind(t_catch,t_loh[catch,])
            }
        }
        
    }
    
    if (length(g_catch)>0){
        gt[[paste0(family,"_G__",germline)]]<-g_catch
        gt[[paste0(family,"_T__",tumor)]]<-t_catch
    } else {
        gt[[paste0(family,"_G__",germline)]]<-data.frame()
        gt[[paste0(family,"_T__",tumor)]]<-data.frame()
    }
    
    
    print(i)
    print(family)
}

gt[["Unknown Genes"]]<-data.frame(cbind(unknown,unknown_fam))

WriteXLS(gt,"./Results_gSNV-tLOH/g-t_gSNV-tLOH_FAMCOLON_11.xlsx")



