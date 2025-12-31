
source("R Code\\ISAR_PLANTS_SOURCE.R")
library(maps)
##########################################################################
#################### GLMM spatial with glmmfields ########################
##########################################################################

options(mc.cores = 4)


##12 knots is arguably quite low for a global dataset, but I have tested increasing to 25 (20,000 iterations)
#and the slope parameter and CIs (which is what we are interested in here) are almost identical: 
#0.73 (0.71-0.75)
m_spatial <- glmmfields::glmmfields(native_count ~ LogArea,
                                    data = datAll,family=nbinom2(link = "log"),
                                    lat = "latitude", lon = "longitude", 
                                    nknots = 12, iter = 2000, chains = 2,
                                    prior_intercept = student_t(3, 0, 10),
                                    prior_beta = student_t(3, 0, 3),
                                    prior_sigma = half_t(3, 0, 3),
                                    prior_gp_theta = half_t(3, 0, 10),
                                    prior_gp_sigma = half_t(3, 0, 3),
                                    seed = 123)

#still warnings
m_spatial2 <- glmmfields::glmmfields(endemic_count ~ LogArea,
                                     data = datAllEndZer,family=nbinom2(link = "log"),
                                     lat = "latitude", lon = "longitude", 
                                     nknots = 30, 
                                     iter = 50000, chains = 4,
                                     prior_intercept = student_t(3, 0, 10),
                                     prior_beta = student_t(3, 0, 3),
                                     prior_sigma = half_t(3, 0, 3),
                                     prior_gp_theta = half_t(3, 0, 10),
                                     prior_gp_sigma = half_t(3, 0, 3),
                                     seed = 123)

 m_list <- list(m_spatial, m_spatial2)
# saveRDS(m_list, file = "spatial_models.rds")
# m_list <- readRDS("D:\\documents\\Work\\On-going projects\\Global Plant SAR\\Global Plant SAR\\Results\\All\\spatial_models.rds")
# m_spatial <- m_list[[1]]; m_spatial2 <- m_list[[2]]

##Info from models
m_list[[1]]
m_list[[2]]
