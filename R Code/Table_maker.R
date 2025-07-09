# #################################################
# ###Main Paper Table 1#############################
# ###################################################

##NEEDS TO BE RUN AFTER THE ABOVE HAS BEEN RUN FOR ALL 
#AND OCEANIC

ResAll <- readRDS("D:\\documents\\Work\\On-going projects\\Global Plant SAR\\Global Plant SAR\\Results\\All\\ResAll_All_islands.RDS")
ResAllEndZer <- readRDS("D:\\documents\\Work\\On-going projects\\Global Plant SAR\\Global Plant SAR\\Results\\All\\ResAllEndZer_All_islands.RDS")
ResOceAll <- readRDS("D:\\documents\\Work\\On-going projects\\Global Plant SAR\\Global Plant SAR\\Results\\Oceanic\\ResAll_Oceanic_islands.RDS")
ResOceEndZer <- readRDS("D:\\documents\\Work\\On-going projects\\Global Plant SAR\\Global Plant SAR\\Results\\Oceanic\\ResAllEndZer_Oceanic_islands.RDS")
ResArchAll <- readRDS("D:\\documents\\Work\\On-going projects\\Global Plant SAR\\Global Plant SAR\\Results\\Archipelago\\ResAll_Archipelagos.RDS")
ResArchEndZer <- readRDS("D:\\documents\\Work\\On-going projects\\Global Plant SAR\\Global Plant SAR\\Results\\Archipelago\\ResAllEndZer_Archipelagos.RDS")

Mod1 <- ResAll[[2]][[2]][[1]][[5]]
Mod2 <- ResAllEndZer[[2]][[2]][[1]][[5]]
Mod3 <- ResOceAll[[2]][[2]][[1]][[5]]
Mod4 <- ResOceEndZer[[2]][[2]][[1]][[5]]
Mod5 <- ResArchAll[[2]][[2]][[1]][[5]]
Mod6 <- ResArchEndZer[[2]][[2]][[1]][[5]]

zs <- c(Mod1$coefficients[2],
        Mod2$coefficients[2],
        Mod3$coefficients[2],
        Mod4$coefficients[2],
        Mod5$coefficients[2],
        Mod6$coefficients[2])

cs <- c(Mod1$coefficients[1],
        Mod2$coefficients[1],
        Mod3$coefficients[1],
        Mod4$coefficients[1],
        Mod5$coefficients[1],
        Mod6$coefficients[1])

rs <- c(summary(Mod1)$r.squared,
        summary(Mod2)$r.squared,
        summary(Mod3)$r.squared,
        summary(Mod4)$r.squared,
        summary(Mod5)$r.squared,
        summary(Mod6)$r.squared)


Ns <- sapply(list(Mod1, Mod2, Mod3, Mod4, Mod5, Mod6),
             function(x){length(x$residuals)})

Tab1 <- data.frame("N" = Ns,
                   "z" = round(zs, 2),
                   "C" = round(cs, 2),
                   "R2" = round(rs, 2))
rownames(Tab1) <- c("All_all", "All_end",
                    "Oce_all", "Oce_end",
                    "Arch_all", "arch_end")

Tab1$IS <- round(c(10^(-cs[1] / zs[1]),
                   10^(-cs[2] / zs[2]),
                   10^(-cs[3] / zs[3]),
                   10^(-cs[4] / zs[4]), NA, NA),2)

write.csv(Tab1, file = "Tab1.csv")


#########Convert model output tables into single###########
##All
rawT <- read.csv("D:\\documents\\Work\\On-going projects\\Global Plant SAR\\Global Plant SAR\\Results\\All\\Raw_mod_output_All_islands.csv")
spaT <- read.csv("D:\\documents\\Work\\On-going projects\\Global Plant SAR\\Global Plant SAR\\Results\\All\\Spatial_mod_output_All_islands.csv")
mixT <- read.csv("D:\\documents\\Work\\On-going projects\\Global Plant SAR\\Global Plant SAR\\Results\\All\\Mixed_mod_output_All_islands.csv")
All_tf <- tab_format(rawT, spaT, mixT)
write.csv(All_tf, file = "Table_S1_Allislands.csv")

##Oceanic
rawT <- read.csv("D:\\documents\\Work\\On-going projects\\Global Plant SAR\\Global Plant SAR\\Results\\Oceanic\\Raw_mod_output_Oceanic_islands.csv")
spaT <- read.csv("D:\\documents\\Work\\On-going projects\\Global Plant SAR\\Global Plant SAR\\Results\\Oceanic\\Spatial_mod_output_Oceanic_islands.csv")
mixT <- read.csv("D:\\documents\\Work\\On-going projects\\Global Plant SAR\\Global Plant SAR\\Results\\Oceanic\\Mixed_mod_output_Oceanic_islands.csv")
Oce_tf <- tab_format(rawT, spaT, mixT)
write.csv(Oce_tf, file = "Table_S1_Oceanic.csv")

##Arch
rawT <- read.csv("D:\\documents\\Work\\On-going projects\\Global Plant SAR\\Global Plant SAR\\Results\\Archipelago\\Raw_mod_output_Archipelagos.csv")
spaT <- read.csv("D:\\documents\\Work\\On-going projects\\Global Plant SAR\\Global Plant SAR\\Results\\Archipelago\\Spatial_mod_output_Archipelagos.csv")
mixT <- NULL
Arch_tf <- tab_format(rawT, spaT, mixT)
write.csv(Arch_tf, file = "Table_S1_Archipelagos.csv")


