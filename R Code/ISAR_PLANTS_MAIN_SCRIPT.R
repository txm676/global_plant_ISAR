
##Note is fixed to use log10

source("R Code\\ISAR_PLANTS_SOURCE.R")

##SETTINGS FOR THRESHOLD MODELS
##set the interval
intz_set <- 0.1
##number of cores; found more than this makes it slower
Ncore <- 8
#point size
PS <- 1.6

##Select dataset: all islands ("All"), oceanic ("Oce"), 
#archipelago ("Arch")
dd <- "All"

if (dd == "All"){
  datM <- datAll
  FN <- c("All_islands")
} else if (dd == "Oce"){
  datM <- datO
  FN <- c("Oceanic_islands")
} else {
  datM <- arch
  FN <- c("Archipelagos")
}

datMEndZer <- filter(datM, endemic_count > 0)#endemics > 0 richness

####################################################
####MODEL COMPARISON AND THRESHOLD MODELS#########
####################################################

#######################################
###All islands, all species
##########################################
ResAll <- fit_thr(datz = select(datM, area, native_count),
                      logBase = log10,
                      parallelz = TRUE,
                      coresz = Ncore,
                      intz = intz_set)

saveRDS(ResAll, file = paste0("Results/ResAll_",FN,".RDS"))

if (dd == "All"){
##Fit of the power model in semi-log space
jpeg("Results/semi_log_power.jpeg", width =18, height = 13,
     res = 300, units = "cm")
pow_plot(ResAll)
dev.off()
}

##Make plot with continents included
#yscale warning is fine; it just relates
#to extrap expanding ylim to fit the extrapolated
#linear model in.
if (dd == "Arch"){
  dco <- NULL
  xro <- dd
  lo <- FALSE
  pf <- FALSE
} else if (dd == "All") {
  dco <- dat_cont
  xro <- TRUE
  lo <- TRUE
  pf <- FALSE
} else {
  dco <- dat_cont
  xro <- TRUE
  lo <- FALSE
  pf <- TRUE
}

RAP3 <- plot_thr(x = ResAll[[2]][[2]],
                 datz = datM,
                 type = "both",
                 pow_fit = NULL,
                 lin_fit = pf,
                 point_end = 3,
                 point_cont = dco,
                 extrap = TRUE,
                 endemic_cont = FALSE,
                 point_size = PS,
                 vir_pal = "inferno",
                 largest = lo,
                 xRaw = xro)

#RAP3[[2]]
AC_LL <- RAP3[[2]]  + ggtitle("a)")

##Correlation between residuals and percentage endemism
M6 <- ResAll[[2]][[2]][[1]][[5]]
M6d <- datM[order(datM$area),]
M6d$area <- log10(M6d$area)
if (!identical(M6d$area, M6$model$x)) stop("Gorgeous")
C6 <- cor.test(residuals(M6), M6d$PercEnd, 
         method = "spearman")

####################################################
###All islands, endemic species (>0 richness islands)
#########################################################

ResAllEndZer <- fit_thr(datz = select(datMEndZer, 
                                      area, endemic_count),
                     logBase = log10,
                     parallelz = TRUE,
                     coresz = Ncore,
                     intz = intz_set)

saveRDS(ResAllEndZer, file = paste0("Results/ResAllEndZer_",FN,".RDS"))

RAEZ3 <- plot_thr(x = ResAllEndZer[[2]][[2]],
                 datz = datMEndZer,
                 type = "both",
                 pow_fit = NULL,
                 lin_fit = TRUE,
                 point_end = 3,
                 point_cont = dco,
                 extrap = TRUE,
                 endemic_cont = TRUE,
                 point_size = PS,
                 vir_pal = "inferno",
                 largest = FALSE,
                 xRaw = xro)

EZC_LL <- RAEZ3[[2]] + ggtitle("b)")

##Correlation between residuals and percentage endemism
M7 <- ResAllEndZer[[2]][[2]][[1]][[5]]
M7d <- datMEndZer[order(datMEndZer$area),]
M7d$area <- log10(M7d$area)
if (!identical(M7d$area, M7$model$x)) stop("Gorgeous")
C7 <- cor.test(residuals(M7), M7d$PercEnd, 
         method = "spearman")

#Export correlation results
CTb <- as.data.frame(round(matrix(c(C6$estimate, C6$p.value,
                C7$estimate, C7$p.value),
              ncol = 2, byrow = TRUE), 2))
colnames(CTb) <- c("Rho", "P")
rownames(CTb) <- c("All", "End")
write.csv(CTb, file = paste0("Results/Corr_tab_",FN,".csv"))

###############################################
############FIGURE 1 & TABLE 2#######################
######################################

jpeg(paste0("Results/Figure_1_",FN,".jpeg"), 
     width =30, height = 13,
     res = 300, units = "cm")
ggpubr::ggarrange(AC_LL, EZC_LL, nrow=1, 
                               common.legend = T,
                               legend="right")
dev.off()

#Both - Table 2
Tab2 <- bind_rows(ResAll[[1]][[2]][,1:4],
ResAllEndZer[[1]][[2]][,1:4])
Tab2$Th1 <- round(10^Tab2$Th1,5)
Tab2$Th2 <- round(10^Tab2$Th2,5)
write.csv(Tab2, file = paste0("Results/Table2_",FN,".csv"))

######################################################
#######FIGURE 2: SEMI-LOG############################
####################################################

if (dd == "All"){
  
  SLNB <- semi_log_NB(datM, datMEndZer)

  jpeg(paste0("Results/Figure_2_",FN,".jpeg"), width =16, height = 11,
       res = 300, units = "cm")
  print(SLNB[[1]])
  dev.off()
  
  jpeg(paste0("Results/Figure_2_SIE_",FN,".jpeg"), width =16, height = 11,
       res = 300, units = "cm")
  print(SLNB[[2]])
  dev.off()
  
  write.csv(SLNB[[3]], "Results/NB_model_output_All.csv")
  
}#eo if dd == All

######################################################
#######Isolation Analyses#####################
##################################################

#00C19A

gi1 <- parreg(datAll3 = datM, Title = "a)", 
              S = "logS", lcol = "red",
              point_size = PS)

gi2 <- parreg(datAll3 = datMEndZer,
              Title = "b)", S = "logES", lcol = "red",
              point_size = PS)

gi_distri <- ggpubr::ggarrange(gi1[[1]], gi2[[1]], nrow=1, 
                               common.legend = T,
                               legend="right")
ggsave(paste0("Results/Figure_3_",FN,".jpeg"), gi_distri, 
       height = 5.5, width = 11)

#Info for figure caption
f3c <- c(round(summary(gi1[[2]])$coefficients[2,c(1,4)],2),
  round(summary(gi2[[2]])$coefficients[2,c(1,4)],2)) 
write.csv(f3c, file = paste0("Results/Figure_3_cap_",FN,".csv"))

###############################
#########Figure 5#################
###################################
if (dd == "Arch"){
  gi1b <- parreg(datAll3 = datM, Title = "b)", 
                S = "logS", lcol = "red",
                point_size = PS)
  AC_LLb <- AC_LL + ggtitle("a)")
  gi_arch <- ggpubr::ggarrange(AC_LLb, gi1b[[1]],
                               nrow=1, 
                    common.legend = T,
                    legend="right")
  jpeg(paste0("Results/Figure_5_",FN,".jpeg"), 
       width = 22, height = 11,
       res = 300, units = "cm")
  print(gi_arch)
  dev.off()
  
  #lm results for caption
  laa <- lm(logS ~ LogArea, data = datM)
  f4c <- c(length(laa$residuals),
           round(summary(laa)$coefficients[2,c(1,4)],2),
           round(summary(laa)$r.squared,2),
           round(summary(gi1b[[2]])$coefficients[2,c(1,4)],2),
           round(summary(gi1b[[2]])$r.squared,2))
  write.csv(f4c, file = paste0("Results/Figure_5_cap_",FN,".csv"))
}

#########################################################
################FIGURE 4##############################
#######################################################

if (dd %in% c("All", "Oce")){

PEF <- Perc_end(datM)
  
ggsave(paste0("Results/Figure_4_",FN,".jpeg"), 
       PEF[[1]], height = 12, width = 8)

#z-values
write.csv(round(PEF[[2]], 2), 
          file = paste0("Results/Figure_4_caption",FN,".csv"))
}
#################################################################
########Raw, Spatial and Mixed Effect Models: log-log##################
######################################################################

#This fits spatial error models (log-log), and 
#mixed effect models (with random intercept and slopes).
#Both model types are fitted four times, once for all species
#and once for endemics, and for both models with and without isolation
#are fitted (i.e., 8 models in total)
if (dd == "Arch"){
  TT <- spaMM(datM, datMEndZer, arch = TRUE)
} else {
  TT <- spaMM(datM, datMEndZer, arch = FALSE)
}

write.csv(TT[[1]], 
          file = paste0("Results/Raw_mod_output_", FN,".csv"))#Raw model output
write.csv(TT[[2]], 
          file = paste0("Results/Spatial_mod_output_", FN,".csv"))#spatial model output
if (dd != "Arch"){
write.csv(TT[[3]], 
          file = paste0("Results/Mixed_mod_output_", FN,".csv"))#Mixed model output
}

#correlograms
jpeg(paste0("Results/Correlograms_", FN,".jpeg"),
     width = 25, height = 25,
     res = 300, units = "cm")
ggpubr::ggarrange(TT[[4]][[1]], TT[[4]][[2]],
                        TT[[4]][[3]], TT[[4]][[4]], 
                  ncol = 2, nrow = 2,
                  common.legend = T,
                  legend="top")
dev.off()

##Residual plots
NN <- c("Resid_All_area.jpeg", "Resid_All_areaIso.jpeg",
        "Resid_End_area.jpeg", "Resid_End_areaIso.jpeg")
for (i in 1:length(TT[[5]])){
jpeg(paste0("Results/",NN[i], "_", FN,".jpeg"), width = 30, height = 30,
     res = 300, units = "cm")
print(TT[[5]][[i]])
dev.off()
}

#########################################################
##########HISTOGRAM in SI###########################
##############################################################

if (dd == "All"){

g1 <- ggplot(datM, aes(x=area, fill=category)) +
  
  geom_histogram(binwidth=0.3) + scale_x_log10() + theme_classic() +
  
  scale_y_continuous(expand = c(0,0)) + labs(title = "Area distribution",
                                             
                                             x=expression(paste("Area - log"[10], "(km"^2, ")")), y="Count") +
  
  scale_fill_manual(name="Island Category", values=c("#999999", "#E69F00", "#56B4E9"))



g2 <- ggplot(datM, aes(x=native_count, fill=category)) +
  
  geom_histogram(binwidth=0.3) + scale_x_log10() + theme_classic() +
  
  scale_y_continuous(expand = c(0,0)) + labs(title = "Native richness distribution",
                                             
                                             x = expression(paste("Species richness - log"[10])),
                                             
                                             y="Count") +
  
  scale_fill_manual(name="Island Category", values=c("#999999", "#E69F00", "#56B4E9"))


g3 <- ggplot(filter(datM, endemic_count > 0), aes(x=endemic_count, fill=category)) +
  
  geom_histogram(binwidth=0.3) + scale_x_log10() + theme_classic() +
  
  scale_y_continuous(expand = c(0,0)) + labs(title = "Endemic richness distribution",
                                             
                                             x = expression(paste("Species richness - log"[10])),
                                             
                                             y="Count") +
  
  scale_fill_manual(name="Island Category", values=c("#999999", "#E69F00", "#56B4E9"))


g_distri <- ggpubr::ggarrange(g1, g2, g3, ncol=1, 
                              common.legend = T,
                              legend="top")


ggsave("Results/histogram_distri.jpg", g_distri, height = 10, width = 5)

}

###############################################################


# #############################################################
# ##########Check against segmented package results#############
# #########################################################
# 
# library(segmented)
# 
# #semi-log
# dd1 <- data.frame("s" = datAll$native_count,
#                   "a" = log10(datAll$area))
# 
# #log-log
# dd2 <- data.frame("s" = log10(datAll$native_count),
#                   "a" = log10(datAll$area))
# 
# ##For sar_threshold fits, run:
# ResAll <- fit_thr(datz = datAll,
#                   logBase = log10,
#                   parallelz = TRUE,
#                   coresz = Ncore,
#                   intz = intz_set)
# 
# 
# ##NOTE - I have checked the number of parameters used,
# #and they match between segmented and sars (5 for contOne
# #and 7 for contTwo)
# 
# ######TWO THRESHOLD CONTINUOUS
# ###Semi-log
# ##sar_threshold()
# r1 <- ResAll[[2]][[1]][[1]][[3]]#model
# s1 <- ResAll[[1]][[1]]#summary table
# 
# ##Segmented
# #First do a free search of parameter space
# l <- lm(s ~ a, data = dd1)
# o <- segmented(l, seg.Z=~a, npsi=2, 
#                control=seg.control(display=FALSE))
# plot(o)
# o
# BIC(o)
# logLik(o)
# s1["ContTwo", "BIC"]#BIC from SARs
# 
# #We see that the breakpoints are very different to those
# #generated by sar_threshold (s1) but the BIC is much higher (s1),
# #i.e., it looks like the segmented fit is not optimal.
# #We can check we get the same results with segmented, by
# #fixing the parameter values to those from sar_threshold fit,
# #which are stored in s1
# o2 <- segmented(l, seg.Z=~a, npsi=2, 
#                control=seg.control(display=FALSE,
#                                    it.max = 0),
#                psi = c(3.455, 5.455))
# o2
# BIC(o2)
# s1["ContTwo", "BIC"]#BIC from SARs
# 
# #However, what happens here is that there are no returned
# #estimated breakpoints in o2 and BIC is lower than in
# #sar_threshold; but this is because it no longer adds in the
# #two extra parameters for breakpoint searching (as they are
# #fixed), so we instead need to compare the log-likelihood values,
# #which are the same:
# logLik(o2)
# logLik(r1)#this is the stored sars model object
# 
# #So it looks like in this case, segmented just doesn't do a good
# #search of parameter space (with the default settings at least),
# #and so doesn't find the best parameters
# 
# ###log-log
# ##sar_threshold()
# r2 <- ResAll[[2]][[2]][[1]][[3]]#model
# s2 <- ResAll[[1]][[2]]#summary table
# 
# ##Segmented
# l2 <- lm(s ~ a, data = dd2)
# o3 <- segmented(l2, seg.Z=~a, npsi=2, 
#                control=seg.control(display=FALSE))
# o3
# BIC(o3)
# s2["ContTwo", "BIC"]#BIC from SARs
# 
# #BIC is similar to sar_threshold, but again different 
# #parameter estimates. So again, we can fix them in segmented:
# o4 <- segmented(l2, seg.Z=~a, npsi=2, 
#                 control=seg.control(display=FALSE,
#                                     it.max = 0),
#                 psi = c(3.355, 3.655))
# logLik(o4)
# logLik(r2)
# 
# #Now, we have identical log-likelihoods
# 
# ######One THRESHOLD CONTINUOUS
# ###Semi-log
# ##sar_threshold()
# r3 <- ResAll[[2]][[1]][[1]][[1]]#model
# s1
# 
# ##Segmented
# #First do a free search of parameter space
# #l <- lm(s ~ a, data = dd1)
# o5 <- segmented(l, seg.Z=~a, npsi=1, 
#                control=seg.control(display=FALSE))
# o5
# BIC(o5)
# s1["ContOne", "BIC"]#BIC from SARs
# 
# #Fix in segmented
# o6 <- segmented(l, seg.Z=~a, npsi=1, 
#                 control=seg.control(display=FALSE,
#                                     it.max = 0),
#                 psi = c(4.955))
# logLik(o6)
# logLik(r3)
# 
# ###log-log
# ##sar_threshold()
# r4 <- ResAll[[2]][[2]][[1]][[1]]#model
# s2
# 
# ##Segmented
# #l2 <- lm(s ~ a, data = dd2)
# o7 <- segmented(l2, seg.Z=~a, npsi=1, 
#                 control=seg.control(display=FALSE))
# o7
# BIC(o7)
# s2["ContOne", "BIC"]#BIC from SARs
# 
# #BIC is similar to sar_threshold, but again different 
# #parameter estimates. So again, we can fix them in segmented:
# o7 <- segmented(l2, seg.Z=~a, npsi=1, 
#                 control=seg.control(display=FALSE,
#                                     it.max = 0),
#                 psi = c(5.355))
# logLik(o7)
# logLik(r4)
