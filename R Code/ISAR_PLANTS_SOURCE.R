
library(dplyr)
library(sars)
library(ggplot2)
library(ggnewscale)
library(RColorBrewer)
#devtools::install_github("jaredhuling/jcolors")
library(jcolors)
library(nlme)
library(mgcv)
library(lme4)
library(glmmfields)
library(glmmTMB)
library(car)
library(spdep)
library(tidyr)
library(maps)

#set ggplot theme
theme_set(theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(), 
                axis.line = element_line(colour = "black"),
                plot.title = element_text(size=22, hjust = 0,
                                          vjust = 2),
                axis.title = element_text(size = 15),
                legend.text=element_text(size=12),
                legend.title=element_text(size=13)))

dat <- read.csv("Data\\data_global_SAR.csv")


##############################################################################
#####Dataset Filtering################
####################################

#No longer undertaken as datasets have been 
#pre-filtered

# #just remove Greenland
# dat <- dat[-which(dat$geo_entity == "Greenland"),]
# 
# #filter out badly sampled islands
# dat <- filter(dat, suit_geo == 1)
# #remove the single island with 0 species, as this
# #causes errors for percentage endemism.
# dat <- filter(dat, native_count > 0)
# #remove the islands with <0.8 in area prop column
# dat <- filter(dat, area_prop >=0.80)
##############################################################

#create percentage endemic column
dat$PercEnd <- dat$endemic_count / dat$native_count

dat$LogArea <- log10(dat$area)
dat$SqrtDist <- sqrt(dat$dist)
dat$LogIso <- log10(dat$dist)
dat$logS <- log10(dat$native_count)
dat$logES <- log10(dat$endemic_count)

#remove continents but save a version for plot
dat_cont <- filter(dat, entity_class == "continent")
#save version with continents for endemic analysis
dat_cont2 <- dat
#Remove continents for main analyses
dat <- filter(dat, entity_class != "continent")


datAll <- dat
datAllEndZer <- filter(datAll, endemic_count > 0)#endemics > 0 richness

#Oceanic
datO <- filter(dat, category == "oceanic")

##load archipelago dataset
arch <- read.csv("Data\\data_archipelago.csv")

##Filter archipelago dataset: no longer needed as pre-filtered
# arch <- filter(arch, suit_geo_arch == 1) %>%
#   dplyr::select("archipelago_name","suit_geo_arch",
#          "area","dist", "island_num",
#          "island_cat","lat","long",
#          "native_count","endemic_count")
# arch <- filter(arch, archipelago_name != "Greenland")
# colnames(arch)[7] <- "latitude"
# colnames(arch)[8] <- "longitude"

arch$LogArea <- log10(arch$area)
arch$SqrtDist <- sqrt(arch$dist)
arch$LogIso <- log10(arch$dist)
arch$logS <- log10(arch$native_count)
arch$logES <- log10(arch$endemic_count)
arch$PercEnd <- arch$endemic_count / arch$native_count


######################################################################
#########THRESHOLD MODELS and ISAR MODEL COMPARISON##############
#################################################################

####INTERNAL THRESHOLD FITTING FUNCTION

fit_thr_int <- function(datz, logAxesz = "none",
                        logBase = log10,
                        conz = 0.1,
                        parallelz = FALSE,
                        coresz = NULL,
                        intz = NULL, ...){

  if (ncol(datz) != 2) stop("incorrect col N")

  modz <- c("ContOne", "ZslopeOne",
            "ContTwo", "ZslopeTwo")

  s <- sar_threshold(data = datz,
                     mod = modz,
                     interval = intz,
                     logAxes = logAxesz,
                     con = conz,
                     logT = logBase,
                     parallel = parallelz,
                     cores = coresz,
                     nisl = 5, ...)
  return(s)
}

##fit the threshold models to both semi-log and
#log-log data, and return both the trimmed model
#summary tables and the sar_threshold fit objects.
#Returns two lists, the first has the trimmed
#summary tables for semi-log and log-log, and the
#second is a list of the two fit objects.
fit_thr <- function(datz,
                    logBase = log10,
                    conz = 0.1,
                    parallelz = FALSE,
                    coresz = NULL,
                    intz = NULL,
                    type = "both", ...){
  #type = "area" or "both"; if former it just does
  #the semi-log
  m1 <- fit_thr_int(datz, logAxesz = "area",
                    logBase = logBase,
                    conz = conz,
                    parallelz = parallelz,
                    coresz = coresz,
                    intz = intz, ...)

  sm1 <- summary(m1)

  rsm1 <- sm1[[2]] %>%
    dplyr::select(BIC, R2, Th1, Th2, seg1, seg2, seg3)

  ##add in power model
  mpow <- sar_power(datz)
  rsm1 <- add_row(rsm1, BIC = round(mpow$BIC, 2),
                  R2 = round(mpow$R2, 2),
          Th1 = NA, Th2 = NA, seg1 = NA,
          seg2 = NA, seg3 = NA)
  if (rownames(rsm1)[7] != "...7") stop("ethi")
  rownames(rsm1)[7] <- "Power"
  rsm1 <- rsm1[order(rsm1$BIC),]

  # #get power values for plotting
  # mpow_vals <- data.frame("A"= mpow$data$A,
  #                         "S" = mpow$data$S,
  #                         "F" = mpow$calculated)
  # mpow_vals <- mpow_vals[order(mpow_vals$A),]

  mpow_vals <- mpow

  if (type == "area"){
    snow <- list(list(rsm1), list(m1), list(mpow_vals))
    return(snow)
  }

  m2 <- fit_thr_int(datz, logAxesz = "both",
                    logBase = logBase,
                    conz = conz,
                    parallelz = parallelz,
                    coresz = coresz,
                    intz = intz, ...)


  sm2 <- summary(m2)

  rsm2 <- sm2[[2]] %>%
    dplyr::select(BIC, R2, Th1, Th2, seg1, seg2, seg3)
  
  #add z-value of log-log power to model table
  lpz <- lin_pow(datz, logT = log10)
  lpz2 <- round(lpz$Model$coefficients[2], 2)
  rsm2$z <- rep(NA, nrow(rsm2))
  rsm2$z[which(rownames(rsm2) == "Linear")] <- lpz2

  snow <- list(list(rsm1, rsm2), list(m1, m2), list(mpow_vals))
  return(snow)
}

###############################################
###Plot all model fits on same plot#########
##########################################

#also option to colour points by % endemism.

#Returns two plots, one with all model fits, and
#one with just best model, linear and, for semi-log,
#the power model

#x = a fit object from our fit_thr function,
#but relating to a fit in a specific log-space,
#datz = Main dataset
#type = whether in semilog (area) or log-log (both) space
#pow_fit = show the power model fit (only for semilog)
#lin_fit = show the fit of the linear model as well as the
#best model (if the linear is not the best model).
#point_end = should points be coloured by %endemism,
#using data provided in datz argument. 0 = no, 1 = yes
#(standard method), 2 = yes (0% points are hollow points),
#3 = yes (same as 2, but logged endemism points)
#point_cont = should the continents also be plotted, with
#names. Either NULL or the continent dataset
#extrap = should the linear model be extrapolated to
#the continents (only works if point_cont included)
#endemic_cont = whether the response variable is native
#richness or endemic richness. Only relevant if continent
#data are being plotted.
#largest = should the nine largest island be labelled;
#only should be run if using all native species and log-log.
#xRaw = should the log10 x-axis area values be plotted on the 
#untransformed scale. Should either be TRUE, FALSE or Arch,
#if archipelago dataset used.
#yraw = should the log10 y-axis richness values be plotted on the 
#untransformed scale. Should either be TRUE, FALSE

plot_thr <- function(x, datz = NULL,
                     type = "area",
                     pow_fit = NULL,
                     lin_fit = FALSE,
                     point_end = 1,
                     point_size = 1.8,
                     point_cont = NULL,
                     extrap = FALSE,
                     endemic_cont = FALSE,
                     vir_pal = "inferno",
                     largest = FALSE,
                     xRaw = FALSE,
                     yRaw = FALSE,
                     ...){

  #if plotting continents and using endemic richness as response,
  #replace native richness with endemic richness for the continent
  #dataset
  if (!is.null(point_cont) & endemic_cont){
    point_cont$native_count <- point_cont$endemic_count
  }

  data <- x[[4]]#nb. this will already be log-transformed if user has selected
  colnames(data) <- c("A", "S")
  #both data and datz have to be in same row order
  data <- data[order(data$A),]
  datz <- datz[order(datz$area),]
  if (data$A[10] != log10(datz$area[10])){
    stop("Seb")
  }
  mods <- x[[1]]
  names <- x[[2]]
  th <- x[[3]]
  names(mods) <- names
  names(th) <- names
  yy <- data$S
  xx <- data$A

  #create new binary % column: 0 = 0%, 1 = >0%
  datz$PercEnd2 <- ifelse(datz$PercEnd == 0, "0", "1")

  cont <- c("ContOne","ZslopeOne","ContTwo","ZslopeTwo","Linear","Intercept")
  xypred.cont <- data.frame("x" = seq(min(data$A), max(data$A),
                                      (max(data$A) - min(data$A))/1000))
  for (i in which(names %in% cont)) {
    xypred.cont[, names[i]] <- stats::predict(mods[[i]], xypred.cont)
  }

  #if semi-log space, also plot the power model fitted
  #in untransformed space
  if (type == "area"){
    #convert the logged area to raw areas
    rawArea <- 10^xypred.cont$x
    #predict using the power model
    xypred.cont$Power <- pow_fit$par[1] *
     (rawArea ^ pow_fit$par[2])
    #Make a long version for ggplot
    xyLong <- tidyr::pivot_longer(xypred.cont,
                                  cols = ContOne:Power,
                                  names_to = "Models")
  } else {
    xyLong <- tidyr::pivot_longer(xypred.cont,
                                  cols = ContOne:Intercept,
                                  names_to = "Models")
  }

  ##get max y-value across model fits and observed to
  #set y-axis range (unless provided)
    all_values <- data$S #start with observed richness values
    #then add in any predicted values from cont or disc models (if they are included)
    #The >1 argument is because the first column is just the area values
    all_values <- c(all_values,
                    unlist(xypred.cont[,2:ncol(xypred.cont)]))
    #if adding continents, need to extend y-axis
    #to their richness values
    if (!is.null(point_cont)){
      if (type == "both"){
        point_cont$native_count <- log10(point_cont$native_count)
      }
      all_values <- c(all_values,
                      point_cont$native_count)
    }

    ff <- range(all_values)
    yRange <- c(ff[1], ff[2])

    #if point_end, then colour points by % endemism,
    #which is provided through the datz argument
    if (point_end == 0){
      gTh0 <- ggplot() +
        geom_point(data = data, aes(x = A, y = S),
                   size = point_size)
    } else if (point_end == 1){
      gTh0 <- ggplot() +
        geom_point(data = data,
                   aes(x = A, y = S,
                       fill = datz$PercEnd),
                   shape = 21,
                   size = point_size)
    } else if  (point_end == 2){
      gTh0 <- ggplot() +
        geom_point(data = data,
                   aes(x = A, y = S,
                   shape = datz$PercEnd2,
                   fill = datz$PercEnd),
                   size = point_size) +
               scale_shape_manual(values = c(1,21)) +
        guides(shape = "none")
    } else if  (point_end == 3){
      
      #first create two datasets, one with zero end
      #and with > 0 end
      data0 <- data[which(datz$PercEnd == 0),]
      dataG0 <- data[which(datz$PercEnd > 0),]
      datzG0 <- datz[which(datz$PercEnd > 0),]
      
      #Then plot them separately, using new_scale_color
      #to create a new plotting scale
      
      gTh0 <- ggplot() +
        geom_point(data = data0,
                   aes(x = A, y = S),
                       shape = 16,
                 # shape = 1,
                       alpha = 0.6,
                 color = "lightgrey",
                   size = point_size) +
                # new_scale_color() +
               #  new_scale_fill() +
                geom_point(data = dataG0,
                 aes(x = A, y = S,
                     fill= log10(datzG0$PercEnd)),
                     shape = 21,
                     alpha = 0.8,
                 size = point_size) +
        scale_fill_viridis_c(option = vir_pal,
                              breaks = c(-0.30, -1, -2, -2.9),
                              labels = c("50", "10",
                                         "1", "0.001"),
                              name = "log10(% End.)",
                              ...)
    }

    gTh1 <- gTh0 +
      ylim(yRange) +
      new_scale_color() +
      geom_line(data = xyLong,
                aes(x = x, y = value,
                    colour = Models)) #+
     # theme_bw() 
  #    scale_color_brewer(palette = "Dark2") +
    #  scale_fill_jcolors_contin(palette = "pal12")

    
    
  #make second version with just the best model,
  #power and linear models
  sData <- summary(x)
  besTMod <- (sData$Model_table %>% rownames())[1]

  # if (type == "area"){
  # xyLong2 <- filter(xyLong,
  #                   Models %in% c("Power", "Linear",
  #                                 besTMod))
  # } else {
  #   xyLong2 <- filter(xyLong,
  #                     Models %in% c("Linear",
  #                                   besTMod))
  # }
  
  if (lin_fit){
    if ("Linear" %in% besTMod){
      xyLong2 <- filter(xyLong,
                        Models %in% besTMod)
    } else {
      xyLong2 <- filter(xyLong,
                        Models %in% c("Linear",
                                      besTMod))
      xyLong2$Models <- factor(xyLong2$Models, 
                               levels = c("Linear", besTMod))
      
    }
  } else {

  #Or just show the best model, not linear and power (unless best)
  if (type == "area"){
    xyLong2 <- filter(xyLong,
                      Models %in% besTMod)
  } else {
    xyLong2 <- filter(xyLong,
                      Models %in% besTMod)
  }
  }#eo if lin_fit

  gTh2 <- gTh0  + 
    new_scale_color() +
    ylim(yRange) +
    geom_line(data = xyLong2,
              aes(x = x, y = value,
                  color = Models),
              linewidth = 1.1) +
   # theme_bw() +
    scale_color_brewer(palette = "Dark2") 
  
  #if adding continents
  if (!is.null(point_cont)){
    point_cont$geo_entity <- stringr::str_to_title(point_cont$geo_entity)

    xRange <- range(data$A)
    xRange[2] <- max(log10(point_cont$area)) + 2
    gTh2 <- gTh2  +
    #  new_scale("shape") +
      xlim(xRange) +
      geom_point(data = point_cont,
                 aes(x = log10(area),
                     y = native_count,
                     shape = geo_entity),
                 color = "darkgreen",
                 size = 2.5) +
      labs(shape = "Continent")
    #if extrapolating the linear model to the continents,
    #we need to extract the largest of the island areas,
    #and the largest of the continent areas, then predict the
    #richness of these two end points
    if (extrap){
      extrapDF2 <- NULL
      maxObsA <- max(data$A)
      maxContA <- max(log10(point_cont$area))
      
      #if native richness, extrapolate linear model, but for
      #endemic richness, extrapolate the best model as this will
      #be a threshold model (and now also the linear model)
        # if (!endemic_cont){
        #   bm2 <- "Linear"
        # } else 
        if (besTMod == "Linear" | lin_fit == FALSE){
          bm2 <- besTMod
        } else {
          bm2 <- c(besTMod, "Linear")
          if (length(bm2) > 2) stop("far away eyes")
        }
        
        extrapR <- sapply(bm2, function(z){
          stats::predict(mods[[z]],
                         data.frame("x" = c(maxObsA,maxContA)))
        })
        
        if (lin_fit & besTMod != "Linear"){
        #Need to put linear model column first, so it matches
        #the correct colour
        extrapR <- extrapR[,c("Linear", besTMod)]
        extrapDF2 <- data.frame("x1" = maxObsA, "x2" = maxContA,
                                "y1" =  extrapR[1,2], 
                                "y2" =  extrapR[2,2])
        }#eo if lin_fit
        
        extrapDF <- data.frame("x1" = maxObsA, "x2" = maxContA,
                               "y1" =  extrapR[1,1], 
                               "y2" =  extrapR[2,1])

      #extend the y-range if predicted continent richness
      #is large
        if (yRange[2] < max(extrapR)){
        yRange[2] = max(extrapR)
      }
      if (yRange[1] > min(extrapR)){
          yRange[1] = min(extrapR)
        } 
        
      #first 2 colours from used colour-pallete
      extrCol <- brewer.pal(3,"Dark2")[1:2]
      
      gTh2 <- gTh2 + geom_segment(data = extrapDF,
                          aes(x = x1, xend = x2,
                           y = y1, yend = y2),
                          linetype = "dashed",
                          colour = extrCol[1],
                          linewidth = 1.1) +
        ylim(yRange)
      
      if (length(extrapDF2) > 1){
        gTh2 <- gTh2 + geom_segment(data = extrapDF2,
                                    aes(x = x1, xend = x2,
                                        y = y1, yend = y2),
                                    linetype = "dashed",
                                    colour = extrCol[2],
                                    linewidth = 1.1) +
          ylim(yRange)
      }

      }#eo if extrap

  }
  gTh1 <- gTh1 + xlab("Area") + ylab("Species richness")
  gTh2 <- gTh2 + xlab("Area") + ylab("Species richness")
  
  #if adding on the 9 largest islands
  if (largest){
    if (nrow(datz) > 500){
       sr <- c("native_count","logS")
    } else {
      sr <- c("endemic_count","logES")
    }
    
    datz9 <- datz[order(datz[[sr[1]]], decreasing = TRUE),]
    largest <-  datz9[1:9,]
    largest$rank <- (max(rank(largest[[sr[1]]]))+1) -
      rank(largest[[sr[1]]])
    
    gTh2 <- gTh2 +
      ggrepel::geom_text_repel(data=largest, 
                               mapping=aes(x=LogArea, 
                                           y=!!sym(sr[2]),
                                           label=rank),
                               force= 3,
                               nudge_x = -1.5,
                               hjust= 0,
                               segment.size = 0.2,
                               box.padding = 0.35,
                               size=3)
  }#eo if largest
  
  ##Convert x-axis values to raw scale
  if (!is.null(xRaw)){
    xra <- ifelse(xRaw == "Arch", TRUE, FALSE)
    gTh2 <- x2r(gTh2, cont = point_cont,
                arch = xra)
    gTh2 <- gTh2 + xlab(expression(paste("Area (km"^2,")")))
  }
  ##Convert y-axis values to raw scale
  if (!is.null(yRaw)){
    gTh2 <- y2r(gTh2)
    gTh2 <- gTh2 + ylab(expression(paste("Species richness")))
  }
  gTh4 <- list(gTh1, gTh2)
  return(gTh4)
}#eo function

# ###quick power plot function
pow_plot <- function(ResAll){
  plot(log10(ResAll[[3]][[1]]$data$A),
       ResAll[[3]][[1]]$data$S,
       xlab = "Log(Area)",
       ylab = "Species richness")
  lines(log10(ResAll[[3]][[1]]$data$A),
        ResAll[[3]][[1]]$calculated,
        col = "red", lwd = 3)
}

##############################################
#########Function to convert log10 x-axis to ####
#########raw scale##############################
#######################################################

#NB. only works properly with ISAR (not isolation) plots
#cont = continents included on the plot or not (if null, it
#is assumed continents are not on plot)
#arch = is it plotting the archipelago dataset
x2r <- function(gobj, cont = NULL, arch = FALSE){
  
  #  xr <- layer_scales(gTh2)$x$range$range
  
  if (!arch){
    
    if (!is.null(cont)){
      
      gobj <- gobj + scale_x_continuous(breaks = c(-5,-2,1,4,7,10,13),
                                        labels = c(expression(paste("10"^-5)), 
                                                   expression(paste("10"^-2)),
                                                   expression(paste("10"^1)),
                                                   expression(paste("10"^4)),
                                                   expression(paste("10"^7)),
                                                   expression(paste("10"^10)),
                                                   expression(paste("10"^13))))
    } else {
      gobj <- gobj + 
        scale_x_continuous(breaks = c(-5,-2,1,4,7),
                           labels = c(expression(paste("10"^-5)), 
                                      expression(paste("10"^-2)),
                                      expression(paste("10"^1)),
                                      expression(paste("10"^4)),
                                      expression(paste("10"^7))))
    }
  } else {
    gobj <- gobj + 
      scale_x_continuous(breaks = c(-2,0,2,4,6),
                         labels = c(expression(paste("10"^-2)), 
                                    expression(paste("10"^0)),
                                    expression(paste("10"^2)),
                                    expression(paste("10"^4)),
                                    expression(paste("10"^6))))
  }
  return(gobj)
}

##Equivalent version for y-axis
y2r <- function(gobj, cont = NULL){
  
      gobj <- gobj + scale_y_continuous(breaks = c(-1,0,1,2,3,4,5),
                                        labels = c(expression(paste("10"^-1)),
                                                   "0", 
                                                   expression(paste("10"^1)),
                                                   expression(paste("10"^2)),
                                                   expression(paste("10"^3)),
                                                   expression(paste("10"^4)),
                                                   expression(paste("10"^5))))
                                                   

  return(gobj)
}


######################################################################
#########PARTIAL REGRESSION: ISOLATION##############
#################################################################

parreg <- function(datAll3, Title, S, lcol = "red",
                   point_size){
  
  S9 <- datAll3[,S]
  
  m2 <- lm(S9 ~ LogArea + LogIso, data=datAll3)
  
  #avoid printing the plot
  ff <- tempfile()
  png(filename=ff)
  partReg <- car::avPlot(m2, "LogIso")
  dev.off()
  
  datAll3$partLogS <- partReg[,2]
  datAll3$partLogIso <- partReg[,1]
  
  data0 <- datAll3[which(datAll3$PercEnd == 0),]
  dataG0 <- datAll3[which(datAll3$PercEnd > 0),]
  
  g_iso <- ggplot() + geom_point(data = data0, 
                                 aes(x = partLogIso, y = partLogS),
                                 shape = 16,
                                 alpha = 0.8, 
                                 color = "lightgrey", 
                                 size = point_size)  +
    geom_point(data = dataG0, aes(x = partLogIso, 
                                  y = partLogS, 
                                  fill= log10(PercEnd)),
               shape = 21, alpha = 0.8, size = point_size) +
    scale_fill_viridis_c(option = "inferno", breaks = c(-0.30, -1, -2, -2.9),
                         labels = c("50", "10","1", "0.001"),
                         name = "log10(% End.)") + 
    labs(title = Title, 
         x=expression(paste("Isolation (km) - log"[10], " | Area")), 
         y = expression(paste("Species richness - log"[10], " | Area"))) + 
    geom_smooth(method="lm", se=F, size=1.1, color=lcol, 
                data = datAll3, mapping=aes(x = partLogIso, 
                                            y = partLogS)) +
    theme(plot.title = element_text(size=22, hjust = 0,
                                    vjust = 2)) + 
    theme(axis.title = element_text(size = 15)) +
    theme(legend.text=element_text(size=12),
          legend.title=element_text(size=13))
  
  #lm results for caption
  li <- lm(partLogS ~ partLogIso, data = datAll3)
  
  ri <- list(g_iso, li)
  return(ri)
  
}#eo function

######################################################################
#########SPATIAL AUTOREGRESSIVE MODEL##############
#################################################################

# function for spatial autoregressive for area #


spatial_model <- function(x, y, area, dist=NULL, S)
{
  nb_list <- nb_object(x, y)
  
  if(is.null(dist))
    
  {mod <- spatialreg::errorsarlm(S ~ area, listw=nb_list, zero.policy=T) }
  
  else 
    
  {mod <- spatialreg::errorsarlm(S ~ area + dist, listw=nb_list, zero.policy=T) }
  
  return(list("nb" = nb_list, "model" = mod))
}


nb_object <- function(x, y){
  coords <- data.frame(x=x, y=y)
  dd <- spdep::dnearneigh(coords, 0, 1000, longlat=T)
  nb <- spdep::nb2listw(dd, zero.policy = TRUE, 
                        glist = lapply(spdep::nbdists(dd, coords), 
                                       function(x) 1 - x/max(dist(coords))))
  return(nb)
  
}

######################################################################
#########MIXED MODELs##############
#################################################################

mixed_model <- function(data, dist=NULL)
{
  
  # S9 <- data[,S]
  
  if(is.null(dist)){
    m <- lme(S ~ LogArea, 
            random = ~ 1+LogArea|archipelago_name, 
            correlation = corExp(form=~longitude+latitude),
                                 data = data)
    } else {
    ##Iterations increased to aid in convergence for the endemics
    #model
  m <- lme(S ~ LogArea+LogIso, 
            random = ~ 1+LogArea+LogIso|archipelago_name, 
            correlation = corExp(form=~longitude+latitude),
            control = list("maxIter" = 50000, "msMaxIter" = 50000,
                           "niterEM" = 50000,
                           "msMaxEval" = 50000),
            data=data)
  }
  
  return(m)
}

###########################################################
########FIT SPATIAL AND MIXED EFFECT MODELS AND FORMAT######
########OUTPUT INTO TABLES#########################################
################################################################

##This fits standard, spatial error models (log-log), and 
#mixed effect models (with random intercept and slopes).
#All three model types are fitted four times, once for all species
#and once for endemics, and for both models with and without isolation
#are fitted (i.e., 12 models in total)
#datAll = datM; datAllEndZer = datMEndZer
spaMM <- function(datAll, datAllEndZer, arch = FALSE){
  
  ##Raw models
  # nb_list_all <- nb_object(datAll$longitude, datAll$latitude)
  # nb_list_end <- nb_object(datAllEndZer$longitude, datAllEndZer$latitude)
  # 
  modAr <- lm(logS ~ LogArea, data=datAll)
  RwS1 <- summary(modAr)
  
  modArIso <- lm(logS ~ LogArea + LogIso, data=datAll)
  RwS2 <- summary(modArIso)
  
  modAr_end <- lm(logES ~ LogArea, data=datAllEndZer)
  RwS3 <- summary(modAr_end)
  
  modArIso_end <- lm(logES ~ LogArea + LogIso, data=datAllEndZer)
  RwS4 <- summary(modArIso_end)
  
  ##Build Spatial results table
  RwM <- matrix(as.vector(c(RwS1$coefficients[2,],
                            rep(NA,4),
                            RwS1$adj.r.squared,
                            RwS2$coefficients[2,],
                            RwS2$coefficients[3,],
                            RwS2$adj.r.squared,
                            RwS3$coefficients[2,],
                            rep(NA,4),
                            RwS3$adj.r.squared,
                            RwS4$coefficients[2,],
                            RwS4$coefficients[3,],
                            RwS4$adj.r.squared)),
                ncol = 9, byrow = TRUE)
  RwM <- as.data.frame(round(RwM, 2))
  colnames(RwM) <- c("Area_coef.", "Area_SE", 
                     "Area_t", "Area_p",
                     "Iso_coef.", "Iso_SE", 
                     "Iso_t", "Iso_p",
                     "Adjusted_R2")
  rownames(RwM) <- c("AllSp._Area",
                     "AllSp._AreaIso",
                     "EndSp._Area",
                     "EndSp._AreaIso")
  
  ##Spatial models
  spatial_loglog_isar <- spatial_model(x=datAll$longitude, 
                                       y=datAll$latitude,
                                       area=datAll$LogArea,
                                       S=datAll$logS)
  
  S1 <- summary(spatial_loglog_isar$model, Nagelkerke = TRUE)
  
  # model area + isolation
  spatial_loglog_ArIso <- spatial_model(x=datAll$longitude, 
                                        y=datAll$latitude,
                                        area=datAll$LogArea,
                                        S=datAll$logS,
                                        dist=datAll$LogIso)
  S2 <- summary(spatial_loglog_ArIso$model, Nagelkerke = TRUE)
  
  ##Endemic versions
  spatial_loglog_isar_end <- spatial_model(x=datAllEndZer$longitude, 
                                           y=datAllEndZer$latitude,
                                           area=datAllEndZer$LogArea,
                                           S=datAllEndZer$logES)
  
  S3 <- summary(spatial_loglog_isar_end$model, Nagelkerke = TRUE)
  
  spatial_loglog_ArIso_end <- spatial_model(x=datAllEndZer$longitude, 
                                            y=datAllEndZer$latitude,
                                            area=datAllEndZer$LogArea,
                                            S=datAllEndZer$logES,
                                            dist=datAllEndZer$LogIso)
  S4 <- summary(spatial_loglog_ArIso_end$model, Nagelkerke = TRUE)
  
  ##Build Spatial results table
  RSM <- matrix(as.vector(c(S1$Coef[2,],rep(NA,4),
                            S1$NK,
                            S2$Coef[2,],S2$Coef[3,],
                            S2$NK,
                            S3$Coef[2,],rep(NA,4),
                            S3$NK,
                            S4$Coef[2,],S4$Coef[3,],
                            S4$NK)),
                ncol = 9, byrow = TRUE)
  RSM <- as.data.frame(round(RSM, 2))
  colnames(RSM) <- c("Area_coef.", "Area_SE", 
                     "Area_z", "Area_p",
                     "Iso_coef.", "Iso_SE", 
                     "Iso_z", "Iso_p",
                     "Pseudo_R2")
  rownames(RSM) <- c("AllSp._Area",
                     "AllSp._AreaIso",
                     "EndSp._Area",
                     "EndSp._AreaIso")
  
  if (!arch){
  ##Mixed effect models
  dat_MM <- datAll
  dat_MM$S <- dat_MM$logS
  mixmod_isar <- mixed_model(data=dat_MM, dist=NULL)
  mixmod_ArIso <- mixed_model(data=dat_MM, dist=TRUE)

  dat_MM2 <- datAllEndZer
  dat_MM2$S <- dat_MM2$logES
  mixmod_isar_end <- mixed_model(data=dat_MM2, dist=NULL)
  mixmod_ArIso_end <- mixed_model(data=dat_MM2, dist=TRUE)
  
  MM1 <- summary(mixmod_isar)
  MM2 <- summary(mixmod_ArIso)
  MM3 <- summary(mixmod_isar_end)
  MM4 <- summary(mixmod_ArIso_end)
  
  MR1 <- piecewiseSEM::rsquared(mixmod_isar)
  MR2 <- piecewiseSEM::rsquared(mixmod_ArIso)
  MR3 <- piecewiseSEM::rsquared(mixmod_isar_end)
  MR4 <- piecewiseSEM::rsquared(mixmod_ArIso_end)
  
  ##Build MM results table
  RMM <- matrix(as.vector(unlist(c(MM1$tTable[2,c(1,2,4,5)], 
                                   rep(NA,4),
                                   MR1[5:6],
                                   MM2$tTable[2,c(1,2,4,5)],
                                   MM2$tTable[3,c(1,2,4,5)],
                                   MR2[5:6],
                                   MM3$tTable[2,c(1,2,4,5)], 
                                   rep(NA,4),
                                   MR3[5:6],
                                   MM4$tTable[2,c(1,2,4,5)],
                                   MM4$tTable[3,c(1,2,4,5)],
                                   MR4[5:6]))), 
                ncol = 10, byrow = TRUE)
  RMM <- as.data.frame(round(RMM, 2))
  colnames(RMM) <- c("Area_coef.", "Area_SE", 
                     "Area_t", "Area_p",
                     "Iso_coef.", "Iso_SE", 
                     "Iso_t", "Iso_p",
                     "Marginal_R2", "Conditional_R2")
  rownames(RMM) <- c("AllSp._Area",
                     "AllSp._AreaIso",
                     "EndSp._Area",
                     "EndSp._AreaIso")
  
  } else {
    mixmod_isar <- NULL
    mixmod_ArIso <- NULL
    mixmod_isar_end <- NULL
    mixmod_ArIso_end <-NULL
    RMM <- NULL
  }
  ##Correlograms
  CO1 <- correl(dat = datM, modLM = modAr, 
         modSpat = spatial_loglog_isar, 
         modMM = mixmod_isar, 
         title = "All - area")
  CO2 <- correl(dat = datM, modLM = modArIso, 
         modSpat = spatial_loglog_ArIso, 
         modMM = mixmod_ArIso, 
         title = "All - area + isolation")
  CO3 <- correl(dat = datAllEndZer, modLM = modAr_end, 
         modSpat = spatial_loglog_isar_end, 
         modMM = mixmod_isar_end, 
         title = "Endemics - area")
  CO4 <- correl(dat = datAllEndZer, modLM = modArIso_end, 
         modSpat = spatial_loglog_ArIso_end, 
         modMM = mixmod_ArIso_end, 
         title = "Endemics - area + isolation")
  
  CO <- list(CO1, CO2, CO3, CO4)
  
  ##Residual plots
  RP1 <- RFP2(modLM = modAr, 
              modSpat = spatial_loglog_isar, 
              modMM = mixmod_isar, 
              title = "All - area")
  RP2 <- RFP2(modLM = modArIso, 
              modSpat = spatial_loglog_ArIso, 
              modMM = mixmod_ArIso, 
              title = "All - area + isolation")
  RP3 <- RFP2(modLM = modAr_end, 
              modSpat = spatial_loglog_isar_end, 
              modMM = mixmod_isar_end, 
              title = "Endemics - area")
  RP4 <- RFP2(modLM = modArIso_end, 
               modSpat = spatial_loglog_ArIso_end, 
               modMM = mixmod_ArIso_end, 
               title = "Endemics - area + isolation")
  
  RPO <- list(RP1, RP2, RP3, RP4)
  
    Rres <- list(RwM, RSM, RMM, CO, RPO)
    
  return(Rres)
}#eo function

########################################################
####Fit Negative binomial GLM: semi-log##############
####Generate plots (fig2) and table###############
##########################################################

semi_log_NB <- function(datM, datMEndZer){
  
  ####NEGATIVE BINOMIAL GLM AND FIGURE 2
  model_nb <- MASS::glm.nb(native_count~LogArea, data=datM)
  
  nbP1 <- performance::r2_efron(model_nb)
  nbP2 <- performance::r2_nagelkerke(model_nb)
  nbP3 <- performance::r2_coxsnell(model_nb)
  nbS1 <- summary(model_nb)
  
  pred.nb <- predict(model_nb, type="response")
  datAll2 <- datM
  datAll2$pred.nb <- pred.nb
  
  datAll3 <- datAll2[order(datAll2$native_count, decreasing = T),]
  datAll3 <- datAll3[1:9,] 
  
  cbb <- c("#999999", "#E69F00", "#56B4E9")
  
  gnb1 <- ggplot() + geom_point(data=datAll2, 
                                mapping=aes(x=LogArea, 
                                            y=native_count, 
                                            fill = category), 
                                size=1.5,
                                pch = 21) + 
    geom_line(data=datAll2, 
              mapping=aes(x=LogArea, y=pred.nb), 
              color="black", linewidth=1) + 
    #theme_bw() +
    xlab(expression(paste("Area (km"^2,")"))) + 
    ylab("Species richness") +
    scale_fill_manual(values = cbb,
                      labels = c("Continental",
                                 "Fragment",
                                 "Oceanic")) +
    labs(fill = "Island Type") +
    geom_text(data = datAll3, 
              label = datAll3$geo_entity,
              nudge_x = -1.5,
              size = 3,
              aes(x = LogArea, y = native_count)) +
    # scale_x_continuous(breaks = -5:6,
    #                    labels = -5:6) + 
    annotate(geom="text", x=-2, y=10000, 
             label="Taiwan   Hainan",
             color="black", size = 3)+ 
    theme(axis.title = element_text(size = 13)) +
    theme(legend.text=element_text(size=11),
          legend.title=element_text(size=12)) +
    #For extending the x-axis for the 10^7 tick mark
    geom_point(aes(x = 7, y=5000), 
               color="white")
  
  gnb1 <- x2r(gnb1, arch = FALSE)

  ####FIGURE S3: Endemics################
  
  model_nb2 <- MASS::glm.nb(endemic_count ~ LogArea, 
                            data=datMEndZer)
  
  nb2P1 <- performance::r2_efron(model_nb2)
  nb2P2 <- performance::r2_nagelkerke(model_nb2)
  nb2P3 <- performance::r2_coxsnell(model_nb2)
  nb2S1 <- summary(model_nb2)
  
  pred.nb2 <- predict(model_nb2, type="response")
  
  datAll5 <- datMEndZer
  datAll5$pred.nb2 <- pred.nb2
  
  gnb2 <- ggplot() + geom_point(data=datAll5, 
                                mapping=aes(x=LogArea, 
                                            y=endemic_count, 
                                            fill = category), 
                                size=1.5,
                                pch = 21) + 
    geom_line(data=datAll5, 
              mapping=aes(x=LogArea, y=pred.nb2), 
              color="black", linewidth=1) + 
    #theme_bw() +
    xlab(expression(paste("Area (km"^2,")"))) +
    ylab("Species richness") +
    scale_fill_manual(values = cbb,
                      labels = c("Continental",
                                 "Fragment",
                                 "Oceanic")) +
    labs(fill = "Island Type")+ 
    theme(axis.title = element_text(size = 13)) +
    theme(legend.text=element_text(size=11),
          legend.title=element_text(size=12))
  gnb2 <- x2r(gnb2, arch = FALSE)

  ##NB summary table
  
  NB_ST <- matrix(as.vector(unlist(c(nbS1$coefficients[1,1],
                                     nbS1$coefficients[2,],
                                     nbP1,nbP2,nbP3,
                                     nb2S1$coefficients[1,1],
                                     nb2S1$coefficients[2,],
                                     nb2P1,nb2P2,nb2P3))), 
                  ncol = 8, byrow = TRUE)
  
  NB_ST <- as.data.frame(round(NB_ST, 2))
  colnames(NB_ST) <- c("Intercept", "Area_coef.", "Area_SE", 
                       "Area_z", "Area_p",
                       "Efron-R2", "Nagelkerke-R2",
                       "CoxSnell-R2")
  rownames(NB_ST) <- c("AllSp",
                       "EndSp")

  return(list(gnb1, gnb2, NB_ST))
  
}#eo semilog function


#########################################################
################FIGURE 4##############################
#######################################################

####Continuous endemism window function##################
#calculate log-log power z, c, R2 and mean area for moving
#window based on %endemism. Starts with all islands, then
#removes islands with < 0.02% endemism, then 0.5%, then 1%,
#then in units of 1% up to 25%.

end_cont <- function(dat){
  
  dat2 <- dat
  x0 <- c(0, 0.002, 0.005)
  x <- c(x0, seq(0.01, 0.25, 0.01))
  res <- matrix(ncol = 10, nrow = length(x))
  colnames(res) <- c("end", "n", "c", "z", "R2", 
                     "p","z_SE", "Median_area", "Mean_area",
                     "Mean_log_area")
  
  for (i in 1:length(x)){
    if (i == 1){
      f <- dat2
    } else {
      f <- filter(dat2, PercEnd > x[i])
    } 
    #   if (nrow(f) < 15) break
    modDum <- lm(logS ~ LogArea, data = f)
    res[i,1] <- x[i] * 100
    res[i,2] <- nrow(f)
    res[i,3:4] <- round(as.vector(modDum$coefficients),3)
    res[i,5] <- round(summary(modDum)$r.squared,3)
    res[i,6] <- summary(modDum)$coefficients[8]
    res[i,7] <- round(summary(modDum)$coefficients[4],3)
    res[i,8] <- round(median(f$area), 2)
    res[i,9] <- round(mean(f$area), 2)
    res[i,10] <- round(mean(f$LogArea), 2)
  }
  
  #  res <- res[1:(i-1),]
  res <- as.data.frame(res)
  
  #add FDR correct p-value column
  res$p_fdr <- p.adjust(p = res$p,
                        method = "fdr")
  
  return(res)
}

###Run end_cont (with and without continents), and save
#results table and build Figure 4. Note that two versions
#of the figure are produced, the smaller version is used to
#extract the legends.
figure4 <- function(dat, dat_cont2,
                    oce = FALSE){
  
  ec1 <- end_cont(dat)
  #with continents included
  ec2 <- end_cont(dat_cont2)
  
  rr <- bind_rows(ec1, ec2)

  rr$Type <- c(rep("Without", nrow(ec1)), 
               rep("With", nrow(ec2)))
  
  # range(ec1$z)
  # range(ec2$z)
  
  if (oce){
    yl4 <- c(0.25, 0.55)
    li4 <- c(10, 710)
    br4 <- c(10, 50, 100, 
             500)
    vl4 <- c("red", "#7DBDE4")
    lb4 <- c("Continent",
             "Oceanic isl.")
  } else {
    yl4 <- c(0.29, 0.55)
    li4 <- c(20, 1270)
    br4 <- c(20, 50, 100, 
             500, 1000)
    vl4 <- c("red","#A4A4A4",
             "#E39F11", "#7DBDE4")
    lb4 <- c("Continent",
             "Continental isl.",
             "Fragment isl.",
             "Oceanic isl.")
  }
  

  
  ge1 <- ggplot(data = rr) + 
    geom_point(aes(x = end, y = z, 
                   col = Type, size = n),
               alpha = 0.5) +
    scale_color_manual(values = c("#009E73", 
                                  "#000000")) +
    ylim(yl4) + 
    xlab("Endemism % cut-off") +
    labs(colour="Continents") +
    labs(size = "No. Isl.") +
    scale_size_continuous(limits =  li4, 
                          breaks = br4) +
    ggtitle("a)") + guides(col="none", size = "none") #+
  #  geom_errorbar(aes(x = end, ymin=z-z_SE, 
  #                   ymax=z+z_SE, col = Type))
  
  ge2 <- ggplot(data = rr) + 
    geom_point(aes(x = end, y = R2, 
                   col = Type, size = n),
               alpha = 0.5) +
    scale_color_manual(values = c("#009E73", "#000000")) +
    xlab("Endemism % cut-off") +
    labs(colour="Continents") +
    labs(size = "No. Isl.") +
    scale_size_continuous(limits =  li4, 
                          breaks = br4) +
    ggtitle(bquote('b)')) +
    ylab(bquote(~R^2)) + guides(col="none", size = "none")
  
  #Mean log area
  ge3 <- ggplot(data = rr) + 
    geom_point(aes(x = end, y = Mean_log_area, 
                   col = Type, size = n),
               alpha = 0.5) +
    scale_color_manual(values = c("#009E73",
                                  "#000000")) +
    scale_size_continuous(limits =  li4, 
                          breaks = br4) +
    xlab("Endemism % cut-off") +
    ylab("Mean log(Area)") +
    guides(col="none", size = "none") +
    ggtitle("c)")
  
  
  #Adding z and R2 to plot. bquote used to make R2 
  #superscript. Note the is.na() warning can be ignored, just
  #relates to use of annotate() with bquote expression
  ge4lm <- lm(logS~LogArea, data = dat_cont2)
  g4tex <- bquote(
    z == .(round(ge4lm$coefficients[2], 2)) * "," ~ R^2 == .(round(summary(ge4lm)$r.squared, 2)))
  
  ge4 <- ggplot(data = dat_cont2) + 
    geom_point(aes(x = LogArea, y = logS, 
                   fill = category),
               size=1.5,
               pch = 21) + 
    scale_fill_manual(values = vl4,
                      labels = lb4) +
    stat_smooth(method = "lm",
                aes(x = LogArea, y = logS),
                se = FALSE, col = "#299375") +
    xlab(expression(paste("Area (km"^2,")"))) +
    ylab("Species richness") +
    ggtitle('d)') + guides(fill="none")+
    annotate("text", x=2, y=5, size = 5,
             label = g4tex)
    #convert axes to untransformed scale
    ge4 <- x2r(ge4, cont = TRUE, arch = FALSE)
    ge4 <- y2r(ge4)

  dat_cont3 <- filter(dat_cont2, PercEnd > 0.05)
  
  #Adding z and R2 to plot. bquote used to make R2 
  #superscript. Note the is.na() warning can be ignored, just
  #relates to use of annotate() with bquote expression
  ge5lm <- lm(logS~LogArea, data = dat_cont3)
  g5tex <- bquote(
    z == .(round(ge5lm$coefficients[2], 2)) * "," ~ R^2 == .(round(summary(ge5lm)$r.squared, 2)))
  
  ge5 <- ggplot(data = dat_cont3) + 
    geom_point(aes(x = LogArea, y = logS, 
                   fill = category),
               size=1.5,
               pch = 21) + 
    scale_fill_manual(values = vl4,
                      labels = lb4) +
    stat_smooth(method = "lm",
                aes(x = LogArea, y = logS),
                se = FALSE, col = "#299375") +
    xlab(expression(paste("Area (km"^2,")"))) +
    ylab("Species richness") +
    ggtitle('e)') + guides(fill="none") +
    annotate("text", x=2, y=4, size = 5,
             label = g5tex)
  #convert axes to untransformed scale
  ge5 <- x2r(ge5, cont = TRUE, arch = FALSE)
  ge5 <- y2r(ge5)
  
  bottom_row <- gridExtra::arrangeGrob(ge4, ge5, ncol = 2)
  
  # Now arrange everything
 f4a <- gridExtra::grid.arrange(
    ge1, ge2, ge3,            # Top row
    bottom_row,                # Bottom row (spanning 3 columns, with 2 equal-width plots)
    ncol = 3,                  # Overall grid: 3 columns
    layout_matrix = rbind(c(1, 2, 3), c(4, 4, 4)))  # 4 fills entire bottom row
  
  ##Versions with legends included
  ge1b <- ggplot(data = rr) + 
    geom_point(aes(x = end, y = z, 
                   col = Type, size = n),
               alpha = 0.5) +
    scale_color_manual(values = c("#009E73", "#000000")) +
    ylim(yl4) + 
    xlab("Endemism % cut-off") +
    labs(colour="Continents") +
    labs(size = "No. Isl.") +
    scale_size_continuous(limits =  li4, 
                          breaks = br4) +
    ggtitle("a)") 
  
  
  ge4b <- ggplot(data = dat_cont2) + 
    geom_point(aes(x = LogArea, y = logS, 
                   fill = category),
               size=1.5,
               pch = 21) + 
    scale_fill_manual(values = vl4,
                      labels = lb4) +
    stat_smooth(method = "lm",
                aes(x = LogArea, y = logS),
                se = FALSE, col = "#299375") +
    xlab("Log(Area)") + ylab("Log(Richness)") +
    ggtitle('d)') + labs(fill="")
  
  # Now arrange everything
  f4b <- gridExtra::grid.arrange(
    ge1b, ge4b, ncol = 2)  # 4 fills entire bottom row
  
  f4l <- list(f4a, f4b, rr)
  return(f4l)
}

#########################################################
################CORRELOGRAMS##############################
#######################################################

#modLM = modAr; modSpat=spatial_loglog_isar; modMM=mixmod_isar

correl <- function(dat, modLM, modSpat, modMM, title,
                   type = "response"){

  # corr_arc_logS <- ncf::correlog(x=dat$longitude, y=dat$latitude, 
  #                                z=log10(dat$native_count), 
  #                                increment=500, resamp=1,
  #                                latlon = TRUE)
  corr_arc_lm <- ncf::correlog(x=dat$longitude, y=dat$latitude, 
                               z=resid(modLM, type = type), 
                               increment=500, resamp=1,
                               latlon = TRUE)
  corr_arc_lm_spatial <- ncf::correlog(x=dat$longitude, y=dat$latitude, 
                                       z=resid(modSpat$model, 
                                               type = type), 
                                       increment=500, resamp=1,
                                       latlon = TRUE)
  if (!is.null(modMM)){
  corr_arc_lm_mixed <- ncf::correlog(x=dat$longitude, y=dat$latitude, 
                                       z=resid(modMM, type = type), 
                                     increment=500, resamp=1,
                                     latlon = TRUE)
  }

  df.correlago_arc <- data.frame(distance_class = corr_arc_lm$mean.of.class, 
                                 LM = corr_arc_lm$correlation, 
                                 LM_spatial = corr_arc_lm_spatial$correlation)
  
  if (!is.null(modMM)){
    df.correlago_arc$Mixed_effects <- corr_arc_lm_mixed$correlation
    df.correlago_arc <- tidyr::pivot_longer(df.correlago_arc, 
                                            cols = c("LM", 
                                                     "LM_spatial",
                                                     "Mixed_effects"))
    
    df.correlago_arc$name <- factor(df.correlago_arc$name, 
                                    levels = c("LM",
                                               "LM_spatial",
                                               "Mixed_effects"))
  } else {
    df.correlago_arc <- tidyr::pivot_longer(df.correlago_arc, 
                                            cols = c("LM", 
                                                     "LM_spatial"))
    
    df.correlago_arc$name <- factor(df.correlago_arc$name, 
                                    levels = c("LM",
                                               "LM_spatial"))
  }
  

  
  glg <- ggplot(df.correlago_arc, aes(x=distance_class, y=value, 
                               color=name)) + 
    xlab("Distance class") + ylab("Moran's I") +
    geom_line() + scale_color_manual(values =  c("#CB2027",
                                                 "#23d0fc",
                                                 "#265DAB",
                                                 "#5DA5DA",
                                                 "#E69F00"),
                                     name="Models") +
    ggtitle(title) + 
    theme(plot.title = element_text(size=18, hjust = 0,
                                    vjust = 2))
    
  return(glg)
}

#################################################
#####Residuals vs fitted plot#################
#################################################

RFP <- function(mod, title, type = "response"){
  
  if (length(mod) > 2){
  R <- resid(mod, type = type)
  Fm <- fitted(mod)
  } else {
    R <- resid(mod$model, type = type)
    Fm <- fitted(mod$model)
  }
  
  DR <- data.frame(Fm, R)
  RP1 <- ggplot(DR) + 
    geom_point(aes(Fm,R)) + xlab("Fitted") +
    ylab("Residuals") +
    ggtitle(title) + 
    theme(plot.title = element_text(size=18, hjust = 0,
                                    vjust = 2))
  
  RP2 <-  ggplot(DR) + geom_histogram(aes(x=R)) +
    xlab("Residuals") +
    ggtitle(title) + 
    theme(plot.title = element_text(size=18, hjust = 0,
                                    vjust = 2))
  
  RP3 <- list(RP1, RP2)
  return(RP3)
}

RFP2 <- function(modLM, modSpat, modMM, title, 
                 type = "response"){
  tit2a <- paste0(title,"_LM")
  RFPa <- RFP(modLM, tit2a, type)
  tit2b <- paste0(title,"_Spatial_LM")
  RFPb <- RFP(modSpat, tit2b, type)
  if (!is.null(modMM)){
  tit2c <- paste0(title,"_Mixed_Model")
  RFPc <- RFP(modMM, tit2c, type)
  RFPg <- ggpubr::ggarrange(RFPa[[1]],RFPa[[2]],
                                  RFPb[[1]],RFPb[[2]],
                                  RFPc[[1]],RFPc[[2]], 
                                  ncol = 2, nrow = 3)
  } else {
    RFPg <- ggpubr::ggarrange(RFPa[[1]],RFPa[[2]],
                              RFPb[[1]],RFPb[[2]],
                              ncol = 2, nrow = 2)
  }
  return(RFPg)
}

####################################################
#########Convert model output tables into single###########
##########Table for the SI##############################
############################################################

tab_format <- function(rawT, spaT, mixT = NULL){
  
  rawT[,c("Pseudo_R2", "Marginal_R2",	"Conditional_R2")] <- NA
  
  spaT[,c("Adjusted_R2", "Marginal_R2",	"Conditional_R2")] <- NA
  spaT <- relocate(spaT, "Adjusted_R2", .before = "Pseudo_R2")
  spaT <- rename(spaT, Area_t = Area_z,
                 Iso_t = Iso_z)
  
  if (!is.null(mixT)){
    mixT[,c("Adjusted_R2", "Pseudo_R2")] <- NA
    mixT <- relocate(mixT, c("Adjusted_R2", "Pseudo_R2"), 
                     .before = c("Marginal_R2",	"Conditional_R2"))
  }
  
  tfr <- bind_rows(rawT, spaT, mixT)
  if (!is.null(mixT)){
  tfr$Type <- c(rep("LM", 4), rep("LM_spatial", 4),
                rep("Mixed_effects", 4))
  } else {
    tfr$Type <- c(rep("LM", 4), rep("LM_spatial", 4))
    tfr[, c("Marginal_R2",	"Conditional_R2")] <- NULL
  }
  tfr <- relocate(tfr, Type, .before = X)
  colnames(tfr)[2] <- "Model"
  return(tfr)
}

########################################################
##For all islands, fit the extra models for Table 1:
#all islands including continents, then just islands
#with >5% SIEs (both with and without continents)
############################################################
tab1_extra <- function(dat_cont2){
  ge4lm <- lm(logS~LogArea, data = dat_cont2)
  dat_cont3 <- filter(dat_cont2, PercEnd > 0.05)
  ge5lm <- lm(logS~LogArea, data = dat_cont3)
  dat_cont4 <- filter(dat_cont3, 
                      category != "complex_origin")
  ge6lm <- lm(logS~LogArea, data = dat_cont3)
  geL <- list("Cont" = ge4lm, 
              "NoCont_Gt0.05" = ge6lm,
              "Cont_Gt0.05" = ge5lm)
  return(geL)
}
