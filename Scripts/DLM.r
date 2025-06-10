library(ggplot2)
library(margins)
library(sandwich)
library(stargazer)
library(lfe)
library(dplyr)
library(lemon) 
library(lspline)
library(fixest)
library(cowplot)
library(texreg)
library(magrittr)


# location
footsteps=4   # 3,4,6 footsteps is lag widths
thresholds=100          
loc_panel <-  paste0("Data/Panel/")
disease="Scarletfever"                    
### working directory
setwd("D://code//Data")
loc_save <- paste0("Data/RegressionResults/",disease,"/E",footsteps,"/",thresholds,"/")
loc_save2<- paste0("Data/significance/",disease,"/E",footsteps,"/",thresholds,"/")
loc_save3 <- paste0("Data/Coefficients/Future/",disease,"/E",footsteps,"/",thresholds,"/")
loc_save4 <- paste0("Data/Coefficients/")
# helper functions

source("HelperFunctions.R")

# read in data
y1 <- 2004                                       
y2 <- 2022
y1_in<-2002
y2_in<-2022
panel <- read.csv(paste(loc_panel,disease,"_ENSO_Growth_Panel_",y1_in,"-",y2_in,".csv",sep=""),fileEncoding="GBK")
panel$pop<-log(panel$pop*0.01)  
panel$year<-as.numeric(panel$year)
scale=""
###filter provinces
#panel<-panel<-panel %>%filter(iso != 140000 & iso != 150000 & iso != 220000 & iso != 310000 & iso != 620000& iso != 630000& iso != 640000& iso != 650000& iso != 230000& iso != 540000& iso != 130000& iso != 120000)
panel %>% filter(year>=y1&year<= y2)


# add vars
panel$t2 <- panel$t**2
panel$p2 <- panel$p**2
##################

############# 1000 bootstraps for heterogeneous effects by teleconnection strength
n_config <-1
configurations <- data.frame("enso_var"=numeric(n_config),  
                             "response_var"=numeric(n_config),
                             "trends"=numeric(n_config),
                             "teleconnection"=numeric(n_config))
configurations[1,] <- c("e_and_c","log_rate","none","t_p_corr_running")
# configurations[2,] <- c("e","log_rate","none","t_p_corr_running")
# configurations[3,] <- c("c","log_rate","none","t_p_corr_running")


## Heterogeneous effects by teleconnection strength
  if (!dir.exists(loc_save)) {
     dir.create(loc_save, recursive = TRUE)} 
  if (!dir.exists(loc_save4)) {
     dir.create(loc_save4, recursive = TRUE)} 
  nlag=24/footsteps
  fe <- "iso+month"            
  cl <- "0"
  nboot <- 1000
  autocorr_lags <- 0    

  
  for (c in c(1:n_config)){                                      
      set.seed(100)                                                 # make sure it's the same bootstraps across variables
      print(paste("configuration: ",c,sep=""))
      enso_var <- configurations[c,"enso_var"]                     
      if ((enso_var=="nino3")|(enso_var=="nino34")){enso_tc<-"e"}else{enso_tc<-enso_var}
      response_var <- configurations[c,"response_var"]
      trnds <- configurations[c,"trends"]
      tc <- configurations[c,"teleconnection"]
      if (trnds=="none"){trends=""}else if(trnds=="linear"){trends=trend_lin}else{trends=""}

     
      panel %>% 
        filter(year>=y1&year<= y2)-> dat
      

      ######## set up dataframe,
      coefs_dist<-data.frame("boot"=c(1:nboot))                   # all coefs
      fix_coefs_dist<-data.frame("boot"=c(1:nboot))               # all fixed effects

      enso_coefs_dist <- data.frame("boot"=c(1:nboot))                  
      exp_interact_coefs_dist <- data.frame("boot"=c(1:nboot))
      for (ll in c(1:nlag)){
        enso_coefs_dist[paste0("coef_lag",ll)] <- numeric(nboot)       
        exp_interact_coefs_dist[paste0("coef_lag",ll)] <- numeric(nboot)  
      }

      if (enso_var=="e_and_c"){
          # secondary set of ENSO effects to get C index as well
          enso_coefs_dist2 <- data.frame("boot"=c(1:nboot))
          exp_interact_coefs_dist2 <- data.frame("boot"=c(1:nboot))
          for (ll in c(1:nlag)){
          enso_coefs_dist2[paste0("coef_lag",ll)] <- numeric(nboot)
          exp_interact_coefs_dist2[paste0("coef_lag",ll)] <- numeric(nboot)
          #####
          coefs_dist[paste0("e_lag",ll)] <- numeric(nboot)       
          coefs_dist[paste0("int_e_lag",ll)] <- numeric(nboot) 
          coefs_dist[paste0("c_lag",ll)] <- numeric(nboot)       
          coefs_dist[paste0("int_c_lag",ll)] <- numeric(nboot)
          }
        }else{
          for (ll in c(1:nlag)){
            coefs_dist[paste0(enso_var,"_lag",ll)] <- numeric(nboot)
            coefs_dist[paste0("int_",enso_var,"_lag",ll)] <- numeric(nboot)
          }
        }

       coefs_dist["pop"] <- numeric(nboot)
       coefs_dist["StringencyIndex_Average"] <- numeric(nboot)

      results_list <- list()  
      all_possible_idx <- c(as.character(unique(panel$iso)),1:12)

      #################loop   bootstrap 
      for (n in c(1:nboot)){
        print(n)
        isos <- as.vector(unique(dat$iso))
        iso_boot <- sample(isos,size=length(isos),replace=T)    
        df_boot <- sapply(iso_boot, function(x) which(dat[,'iso']==x))  
        data_boot <- dat[unlist(df_boot),] 
       
        
        # run model for each set of lags  
        form_dl=""
        if (enso_var=="e_and_c"){
              form_dl <- paste0(response_var,"~","+pop+StringencyIndex_Average")    
              for (j in c(1:nlag)){
                form_dl <- paste0(form_dl," + e_lag",j," + ","e_lag",j,"*",tc,"_e",
                                  " + c_lag",j," + ","c_lag",j,"*",tc,"_c")
                                  }
              }else{
              form_dl <- paste0(response_var," ~ ","+pop+StringencyIndex_Average")
              for (j in c(1:nlag)){
                form_dl <- paste0(form_dl," + ",enso_var,"_lag",j," + ",
                                  enso_var,"_lag",j,"*",tc,"_",enso_tc)
              }
              }
    
        if (autocorr_lags>0){  
          for (j in c(1:autocorr_lags)){form_dl <- paste0(form_dl," + ",response_var,"_lag",j)}
        }
        
        form <- as.formula(paste0(form_dl,trends,"| ",fe," | 0 | ",cl))  #
        mdl <- felm(form,data=data_boot)    
        df = getfe(mdl)  

        # extract fix effects
        df <- df[, c("idx", "effect")]    
        missing_idx <- setdiff(all_possible_idx, df$idx)
       if(length(missing_idx) > 0) {
        df <- rbind(df, data.frame(idx = missing_idx, effect = NA))  
        }
        df_wide <- data.frame(t(df$effect))
        colnames(df_wide) <- df$idx 
        df_wide <- df_wide[, sort(colnames(df_wide))]       
        results_list[[n]] <- df_wide
        


        # extract coefficients
        if (enso_var=="e_and_c"){
            enso_coefs_dist[n,"coef_lag1"] <- as.numeric(coef(mdl)["e_lag1"])
            exp_interact_coefs_dist[n,"coef_lag1"] <- as.numeric(coef(mdl)[paste0("e_lag1:",tc,"_e")])
            enso_coefs_dist2[n,"coef_lag1"] <- as.numeric(coef(mdl)["c_lag1"])
            exp_interact_coefs_dist2[n,"coef_lag1"] <- as.numeric(coef(mdl)[paste0("c_lag1:",tc,"_c")])
            #####
            coefs_dist[n,"e_lag1"] <- as.numeric(coef(mdl)["e_lag1"])
            coefs_dist[n,"int_e_lag1"] <- as.numeric(coef(mdl)[paste0("e_lag1:",tc,"_e")])
            coefs_dist[n,"c_lag1"] <- as.numeric(coef(mdl)["c_lag1"])
            coefs_dist[n,"int_c_lag1"] <- as.numeric(coef(mdl)[paste0("c_lag1:",tc,"_c")])         
            
            for (ll in c(2:nlag)){
            enso_coefs_dist[n,paste0("coef_lag",ll)] <- as.numeric(coef(mdl)[paste0("e_lag",ll)])
            exp_interact_coefs_dist[n,paste0("coef_lag",ll)] <- as.numeric(coef(mdl)[paste0(tc,"_e:e_lag",ll)])
            enso_coefs_dist2[n,paste0("coef_lag",ll)] <- as.numeric(coef(mdl)[paste0("c_lag",ll)])
            exp_interact_coefs_dist2[n,paste0("coef_lag",ll)] <- as.numeric(coef(mdl)[paste0(tc,"_c:c_lag",ll)])
            #####
            coefs_dist[n,paste0("e_lag",ll)] <- as.numeric(coef(mdl)[paste0("e_lag",ll)])
            coefs_dist[n,paste0("int_e_lag",ll)] <- as.numeric(coef(mdl)[paste0(tc,"_e:e_lag",ll)])
            coefs_dist[n,paste0("c_lag",ll)] <- as.numeric(coef(mdl)[paste0("c_lag",ll)])
            coefs_dist[n,paste0("int_c_lag",ll)] <- as.numeric(coef(mdl)[paste0(tc,"_c:c_lag",ll)])           
                                }
                  }

        else {
            enso_coefs_dist[n,"coef_lag1"] <- as.numeric(coef(mdl)[paste0(enso_var,"_lag1")])
            exp_interact_coefs_dist[n,"coef_lag1"] <- as.numeric(coef(mdl)[paste0(enso_var,"_lag1",":",tc,"_",enso_tc)])
            ##### 
            coefs_dist[paste0(enso_var,"_lag1")] <- as.numeric(coef(mdl)[paste0(enso_var,"_lag1")])
            coefs_dist[paste0("int_",enso_var,"_lag1")] <- as.numeric(coef(mdl)[paste0(enso_var,"_lag1",":",tc,"_",enso_tc)])
          for (ll in c(2:nlag)){
            enso_coefs_dist[n,paste0("coef_lag",ll)] <- as.numeric(coef(mdl)[paste0(enso_var,"_lag",ll)])
            exp_interact_coefs_dist[n,paste0("coef_lag",ll)] <- as.numeric(coef(mdl)[paste0(tc,"_",enso_tc,":",enso_var,"_lag",ll)])
            #####
            coefs_dist[paste0(enso_var,"_lag",ll)] <- as.numeric(coef(mdl)[paste0(enso_var,"_lag",ll)])
            coefs_dist[paste0("int_",enso_var,"_lag",ll)] <- as.numeric(coef(mdl)[paste0(tc,"_",enso_tc,":",enso_var,"_lag",ll)])
          }
         }
        #enso_coefs_dist[n,"gdp"] <- as.numeric(coef(mdl)["gdp"])
        coefs_dist[n,"pop"] <- as.numeric(coef(mdl)["pop"])
        coefs_dist[n,"StringencyIndex_Average"] <- as.numeric(coef(mdl)["StringencyIndex_Average"])
      }



      if (enso_var=="e_and_c"){
        # write out ENSO coefs
        fname_enso_dist1 <- paste0(disease,y1,y2,scale,"ENSO_teleconnection-interacted_coefs_bootstrap-dist_lag",nlag,"_e_e-and-c_",
                                  response_var,"_trend",trnds,".csv")
        write.csv(enso_coefs_dist,paste(loc_save,fname_enso_dist1,sep=""))
        print(fname_enso_dist1)
        fname_enso_dist2 <- paste0(disease,y1,y2,scale,"ENSO_teleconnection-interacted_coefs_bootstrap-dist_lag",nlag,"_c_e-and-c_",
                                    response_var,"_trend",trnds,".csv")
        write.csv(enso_coefs_dist2,paste(loc_save,fname_enso_dist2,sep=""))
        print(fname_enso_dist2)
        
        # write out interaction coefs
        fname_int_dist1 <- paste0(disease,y1,y2,scale,"ENSO_teleconnection-interaction_coefs_bootstrap-dist_lag",nlag,"_e_e-and-c_",
                                  response_var,"_trend",trnds,".csv")
        write.csv(exp_interact_coefs_dist,paste(loc_save,fname_int_dist1,sep=""))
        print(fname_int_dist1)
        fname_int_dist2 <- paste0(disease,y1,y2,scale,"ENSO_teleconnection-interaction_coefs_bootstrap-dist_lag",nlag,"_c_e-and-c_",
                                  response_var,"_trend",trnds,".csv")
        write.csv(exp_interact_coefs_dist2,paste(loc_save,fname_int_dist2,sep=""))
        print(fname_int_dist2)
        
        
      }else {
        # write out ENSO coefs
        fname_enso_dist <- paste0(disease,y1,y2,scale,"ENSO_teleconnection-interacted_coefs_bootstrap-dist_lag",nlag,"_",enso_var,"_",
                                  response_var,"_trend",trnds,".csv")
        write.csv(enso_coefs_dist,paste(loc_save,fname_enso_dist,sep=""))
        print(fname_enso_dist)
        
        # write out interaction coefs
        fname_int_dist <- paste0(disease,y1,y2,scale,"ENSO_teleconnection-interaction_coefs_bootstrap-dist_lag",nlag,"_",enso_var,"_",
                                  response_var,"_trend",trnds,".csv")
        write.csv(exp_interact_coefs_dist,paste(loc_save,fname_int_dist,sep=""))
        print(fname_int_dist)
      }
        # write out all coefs
        all_coefs_dist <- paste0(disease,"_bootstrap_coefs_lag",nlag,"_",enso_var,"_",
                                  response_var,".csv")
        write.csv(coefs_dist,paste(loc_save3,all_coefs_dist,sep=""))
        print(all_coefs_dist)

        # write out fixed effects
        final_df <- do.call(rbind, results_list)
        fname_fix_dist <- paste0(disease,"_bootstrap_fix_lag",nlag,"_",enso_var,"_",
                                  response_var,".csv")
        write.csv(final_df,paste(loc_save3,fname_fix_dist,sep=""))
        print(fname_fix_dist)
      }

  ###
  ######################################################## main model and significance test #####################################################
  boot=1
  for (c in c(1:n_config)){                                   #
      print(paste("configuration: ",c,sep=""))
      response_var <- configurations[c,"response_var"]
      print(response_var)
      enso_var <- configurations[c,"enso_var"]
      trnds <- configurations[c,"trends"]
      tc <- configurations[c,"teleconnection"]  
      if ((enso_var=="nino3")|(enso_var=="nino34")){enso_tc<-"e"}else{enso_tc<-enso_var}
      if (trnds=="none"){trends=""}else if(trnds=="linear"){trends=trend_lin}else{trends=""}  

      form_dl_test=""
      if (enso_var=="e_and_c"){
          form_dl_test <- paste0(response_var," ~ ","+pop+StringencyIndex_Average")       
          for (j in c(1:nlag)){
            form_dl_test <- paste0(form_dl_test," + e_lag",j," + ","e_lag",j,"*",tc,"_e",
                              " + c_lag",j," + ","c_lag",j,"*",tc,"_c")
                              }
          }else{
          form_dl_test <- paste0(response_var," ~ ","+pop+StringencyIndex_Average")
          for (j in c(1:nlag)){
            form_dl_test <- paste0(form_dl_test," + ",enso_var,"_lag",j," + ",
                              enso_var,"_lag",j,"*",tc,"_",enso_tc)
          }

          }
    
        form_test <- as.formula(paste0(form_dl_test,trends,"| ",fe," | 0 | ",cl))  #
        mdl_test <- felm(form_test,data=panel)   

        cov_matrix <- vcov(mdl_test)
        coefs <- coef(mdl_test)
        t_critical <- qt(0.975, mdl_test$df.residual)

        model_s<-summary(mdl_test)
        aic_value <- AIC(mdl_test)
        bic_value <- BIC(mdl_test)
        residual_se <- model_s$sigma
        
        if (!dir.exists(loc_save2)) {
          dir.create(loc_save2, recursive = TRUE)
        }

      sink(paste0(loc_save2,y1,"_",y2,"_",response_var,"_",enso_var,"_",nlag,"lag-regression_table.txt"))
      cat("Model Formula:\n",as.character(form_test), "\n\n")
      cat(capture.output(model_s), sep = "\n")
      cat("AIC: ", aic_value, "\n")
      cat("BIC: ", bic_value, "\n")
      cat("Residual Standard Error: ", residual_se, "\n")
      sink()
      

### output model cofficient for caculating current cumlative effects
  enso_coefs_dist <- data.frame("boot"=c(1:boot))
  enso_coefs_dist2 <- data.frame("boot"=c(1:boot))
  exp_interact_coefs_dist <- data.frame("boot"=c(1:boot))
  exp_interact_coefs_dist2 <- data.frame("boot"=c(1:boot))
 if (enso_var=="e_and_c"){
      enso_coefs_dist[n,"coef_lag1"] <- as.numeric(coef(mdl_test)["e_lag1"])
      exp_interact_coefs_dist[n,"coef_lag1"] <- as.numeric(coef(mdl_test)[paste0("e_lag1:",tc,"_e")])
      enso_coefs_dist2[n,"coef_lag1"] <- as.numeric(coef(mdl_test)["c_lag1"])
      exp_interact_coefs_dist2[n,"coef_lag1"] <- as.numeric(coef(mdl_test)[paste0("c_lag1:",tc,"_c")])
    
            
      for (ll in c(2:nlag)){
      enso_coefs_dist[n,paste0("coef_lag",ll)] <- as.numeric(coef(mdl_test)[paste0("e_lag",ll)])
      exp_interact_coefs_dist[n,paste0("coef_lag",ll)] <- as.numeric(coef(mdl_test)[paste0(tc,"_e:e_lag",ll)])
      enso_coefs_dist2[n,paste0("coef_lag",ll)] <- as.numeric(coef(mdl_test)[paste0("c_lag",ll)])
      exp_interact_coefs_dist2[n,paste0("coef_lag",ll)] <- as.numeric(coef(mdl_test)[paste0(tc,"_c:c_lag",ll)])          
                          }
            }

  else {
      enso_coefs_dist[n,"coef_lag1"] <- as.numeric(coef(mdl_test)[paste0(enso_var,"_lag1")])
      exp_interact_coefs_dist[n,"coef_lag1"] <- as.numeric(coef(mdl_test)[paste0(enso_var,"_lag1",":",tc,"_",enso_tc)])

    for (ll in c(2:nlag)){
      enso_coefs_dist[n,paste0("coef_lag",ll)] <- as.numeric(coef(mdl_test)[paste0(enso_var,"_lag",ll)])
      exp_interact_coefs_dist[n,paste0("coef_lag",ll)] <- as.numeric(coef(mdl_test)[paste0(tc,"_",enso_tc,":",enso_var,"_lag",ll)])
    }
    }

##
  if (enso_var=="e_and_c"){
    # write out ENSO coefs
    fname_enso_dist1 <- paste0(disease,y1,y2,scale,"ENSO_teleconnection-interacted_coefs_bootstrap-dist_lag",nlag,"_e_e-and-c_",
                              response_var,"_trend",trnds,".csv")
    write.csv(enso_coefs_dist,paste(loc_save4,fname_enso_dist1,sep=""))
    print(fname_enso_dist1)
    fname_enso_dist2 <- paste0(disease,y1,y2,scale,"ENSO_teleconnection-interacted_coefs_bootstrap-dist_lag",nlag,"_c_e-and-c_",
                                response_var,"_trend",trnds,".csv")
    write.csv(enso_coefs_dist2,paste(loc_save4,fname_enso_dist2,sep=""))
    print(fname_enso_dist2)
    
    # write out interaction coefs
    fname_int_dist1 <- paste0(disease,y1,y2,scale,"ENSO_teleconnection-interaction_coefs_bootstrap-dist_lag",nlag,"_e_e-and-c_",
                              response_var,"_trend",trnds,".csv")
    write.csv(exp_interact_coefs_dist,paste(loc_save4,fname_int_dist1,sep=""))
    print(fname_int_dist1)
    fname_int_dist2 <- paste0(disease,y1,y2,scale,"ENSO_teleconnection-interaction_coefs_bootstrap-dist_lag",nlag,"_c_e-and-c_",
                              response_var,"_trend",trnds,".csv")
    write.csv(exp_interact_coefs_dist2,paste(loc_save4,fname_int_dist2,sep=""))
    print(fname_int_dist2)
    write.csv(cov_matrix,paste(loc_save4,"cov_matrix_",disease,".csv",sep=""),row.names = T)
    write.csv(coefs,paste(loc_save4,"coefs_",disease,".csv",sep=""),row.names = T)
    
  }else {
    # write out ENSO coefs
    fname_enso_dist <- paste0(disease,y1,y2,scale,"ENSO_teleconnection-interacted_coefs_bootstrap-dist_lag",nlag,"_",enso_var,"_",
                              response_var,"_trend",trnds,".csv")
    write.csv(enso_coefs_dist,paste(loc_save4,fname_enso_dist,sep=""))
    print(fname_enso_dist)
    
    # write out interaction coefs
    fname_int_dist <- paste0(disease,y1,y2,scale,"ENSO_teleconnection-interaction_coefs_bootstrap-dist_lag",nlag,"_",enso_var,"_",
                              response_var,"_trend",trnds,".csv")
    write.csv(exp_interact_coefs_dist,paste(loc_save4,fname_int_dist,sep=""))
    print(fname_int_dist)

    write.csv(cov_matrix,paste(loc_save4,"cov_matrix_",disease,".csv",sep=""),row.names = T)
    write.csv(coefs,paste(loc_save4,"coefs_",disease,".csv",sep=""),row.names = T)
  }
  }